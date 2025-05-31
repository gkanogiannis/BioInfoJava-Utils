/*
 *
 * BioInfoJava-Utils 
 *
 * Copyright (C) 2021 Anestis Gkanogiannis <anestis@gkanogiannis.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */
package com.gkano.bioinfo.vcf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;

import com.gkano.bioinfo.var.GeneralTools;
import com.gkano.bioinfo.var.Logger;

@SuppressWarnings("FieldMayBeFinal")
public class VCFManager<T> implements Runnable {

    //private static final String POISON_PILL = "__END__";
    private final List<T> POISON_PILL_BATCH = Collections.emptyList();
    private final int batchSize = 1000;

    private final List<String> inputFileNames;
    private int usingThreads;
    private final int maxSizeOfVariantCache;
    private final Function<String, T> variantParser;
    private boolean verbose = false;

    private BlockingQueue<List<T>> variantRawCache;
    private Map<String, int[]> genotypeEncodingCache;

    private List<String> commentData;
    private String headerData;
    private static AtomicInteger currVariantCount;

    private int numVariants;
    private int numSamples;
    private List<String> sampleNames;
    private int ploidy;
    private int maxAlleles;

    private CountDownLatch startSignal;
    private CountDownLatch doneSignal;
    private ExecutorService pool;
    //private Map<Integer, VariantProcessor<T>> variantProcessors;
    private List<CompletableFuture<ProcessorResult>> variantProcessors;

    public static class ProcessorResult {
        int[][] dotProd;
        int[] norm;
        private ProcessorResult(int[][] localDotProd, int[] localNorm) {
            this.dotProd = localDotProd;
            this.norm = localNorm;
        }
        public void merge(ProcessorResult other) {
            for (int i = 0; i < norm.length; i++) {
                norm[i] += other.norm[i];
                for (int j = i; j < norm.length; j++) {
                    dotProd[i][j] += other.dotProd[i][j];
                }
            }
        }
    }

    public VCFManager(
            List<String> inputFileNames,
            int usingThreads,
            Function<String, T> variantParser,
            boolean verbose
    ) {
        this.inputFileNames = inputFileNames;
        this.usingThreads = usingThreads;
        this.variantParser = variantParser;
        this.verbose = verbose;

        this.maxSizeOfVariantCache = 2 * usingThreads;
    }

    public void init() {
        if (inputFileNames == null || inputFileNames.isEmpty()) {
            Logger.error(this, "No VCF input files provided.");
            System.exit(1);
            return;
        }

        variantRawCache = new LinkedBlockingQueue<>(maxSizeOfVariantCache);
        genotypeEncodingCache = new ConcurrentHashMap<>();
        commentData = new ArrayList<>();
        currVariantCount = new AtomicInteger(0);

        int cpus = Runtime.getRuntime().availableProcessors();
        usingThreads = (cpus < usingThreads ? cpus : usingThreads);
        if (verbose) {
            Logger.info(this, "cpus=" + cpus + "\tusing=" + usingThreads);
        }

        startSignal = new CountDownLatch(1);
        doneSignal = new CountDownLatch(1);

        pool = Executors.newFixedThreadPool(usingThreads);

        variantProcessors = new ArrayList<>();
        
        for (int t = 0; t < usingThreads; t++) {
            CompletableFuture<ProcessorResult> variantProcessor = CompletableFuture.supplyAsync(() -> {
                try {
                    startSignal.await();  // ⏸ wait until ploidy is known
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    throw new RuntimeException(e);
                }

                int[][] localDotProd = new int[numSamples][numSamples];
                int[] localNorm = new int[numSamples];

                List<T> batch;
                int[][] variantEncoded;
                int count, step;
                int[] di, dj;
                int norm, dotProd;
                while (true) {
                    try {
                        batch = variantRawCache.take();
                    } catch (InterruptedException e) {
                        Thread.currentThread().interrupt();
                        break;
                    }
                    if (batch.isEmpty()) {
                        break; // poison pill
                    }
                    for (T variant : batch) {
                        count = currVariantCount.incrementAndGet();
                        if(verbose){
                            step = GeneralTools.getAdaptiveVariantStep(count);
                            if (verbose && count % step == 0) {
                                Logger.infoCarret(this, "Variants:\t" + count);
                            }
                        }
                        try {
                            variantEncoded = SNPEncoder.encodeSNPOneHot(
                                (String) variant, ploidy, maxAlleles, genotypeEncodingCache, numSamples);
                        } catch (IllegalArgumentException e) {
                            continue;
                        }
                        for (int i = 0; i < numSamples; i++) {
                            di = variantEncoded[i];
                            norm = 0;
                            for (int k = 0; k < di.length; k++) norm += Integer.bitCount(di[k] & di[k]);
                            localNorm[i] += norm;
                            for (int j = i; j < numSamples; j++) {
                                dj = variantEncoded[j];
                                dotProd = 0;
                                for (int k = 0; k < di.length; k++) dotProd += Integer.bitCount(di[k] & dj[k]);
                                localDotProd[i][j] += dotProd;
                            }
                        }
                    }
                }
                return new ProcessorResult(localDotProd, localNorm);
            }, pool);
            variantProcessors.add(variantProcessor);
        }
    }

    public void awaitFinalization() throws InterruptedException {
        doneSignal.await();
    }

    @Override
    public void run() {
        try {
            Logger.setVerbose(verbose);
            Logger.info(this, "START READ");

            VCFDecoder decoder = new VCFDecoder();
            //byte[][] line : split at tabs, byte[] is a string between tabs
            //Stringh line : a snp line from vcf
            List<T> batch = new ArrayList<>(1000);
            VCFStreamingIterator<String> iterator = new VCFStreamingIterator<>(decoder, verbose, inputFileNames);
            for (String line : iterator) if (line != null) processVariantLine(line, batch);
            // Push the last partial batch if needed
            if (!batch.isEmpty()) variantRawCache.put(batch);
            // putting POISON pills
            for (int i = 0; i < usingThreads; i++) variantRawCache.put(POISON_PILL_BATCH); // use empty list as poison pill
            
            Logger.info(this, "END READ");

            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
            numVariants = currVariantCount.get();

            if (verbose) Logger.info(this, "Processed variants :\t" + numVariants);
            
            doneSignal.countDown();
            
        } catch (InterruptedException e) {
            Logger.error(this, e.getMessage());
            shutdown();
            System.exit(1);
        } finally {
            shutdown();
        }
    }

    private void processVariantLine(String line, List<T> batch) {
        try {
            if (line.startsWith("##")) {
                commentData.add(line);
            } else if (line.startsWith("#")) {
                headerData = line;
                sampleNames = getSampleNamesFromHeader();
                if (sampleNames == null || sampleNames.size() <= 0) {
                    shutdown();
                    System.exit(1);
                    return;
                }
                numSamples = sampleNames.size();
            } else {
                if (ploidy <= 0) {
                    try {
                        int[] ploidy_maxAlleles = SNPEncoder.guessPloidyAndMaxAllele(line);
                        ploidy = ploidy_maxAlleles[0];
                        maxAlleles = ploidy_maxAlleles[1];
                        if (ploidy > 0) {
                            startSignal.countDown();
                        }
                    } catch (IllegalArgumentException e) {
                        Logger.error(this, e.getMessage());
                        shutdown();
                        System.exit(1);
                        return;
                    }
                }
                T parsed = variantParser.apply(line);
                batch.add(parsed);
                if (batch.size() >= batchSize) {
                    variantRawCache.put(new ArrayList<>(batch));
                    batch.clear();
                }
            }
        } catch (InterruptedException e) {
            Logger.error(this, e.getMessage());
        }
    }

    public float[][] reduceDotProdToDistances() {
        int[][] finalDotProd = new int[numSamples][numSamples];
        int[] finalNorm = new int[numSamples];

        for (CompletableFuture<ProcessorResult> vp : variantProcessors) {
            ProcessorResult r = vp.join();
            for (int i = 0; i < numSamples; i++) {
                finalNorm[i] += r.norm[i];
                for (int j = i; j < numSamples; j++) {
                    finalDotProd[i][j] += r.dotProd[i][j];
                }
            }
        }

        float[][] cosineDist = new float[numSamples][numSamples];
        for (int i = 0; i < numSamples; i++) {
            float normI = (float) Math.sqrt(finalNorm[i]);
            for (int j = i; j < numSamples; j++) {
                float normJ = (float) Math.sqrt(finalNorm[j]);
                float dot = finalDotProd[i][j];
                float similarity = ((normI > 0 && normJ > 0) ? (dot / (normI * normJ)) : 0.0f);
                if (j == i && normI == 0 && normJ == 0) {
                    similarity = 1.0f;
                }
                float dist = 1.0f - similarity;
                if(dist<0) dist = 0f;
                cosineDist[i][j] = cosineDist[j][i] = dist;
            }
        }
        return cosineDist;
    }

    public List<String> getSampleNamesFromHeader() {
        if (headerData == null || !headerData.startsWith("#CHROM")) {
            throw new IllegalArgumentException("Line does not start with #CHROM: " + headerData);
        }

        String[] fields = this.headerData.split("\t", -1);
        if (fields.length <= 9) {
            return Collections.emptyList();
        }

        return Arrays.asList(Arrays.copyOfRange(fields, 9, fields.length));
    }

    private void shutdown() {
        pool.shutdown();
    }

    public List<String> getSampleNames() {
        return sampleNames;
    }

    public int getNumVariants() {
        return numVariants;
    }

    public int getNumSamples() {
        return numSamples;
    }

    public int getPloidy() {
        return ploidy;
    }

    public int getMaxAlleles() {
        return maxAlleles;
    }

}
