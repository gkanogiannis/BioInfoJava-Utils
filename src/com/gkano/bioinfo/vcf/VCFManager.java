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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;

import com.gkano.bioinfo.var.Logger;

@SuppressWarnings("FieldMayBeFinal")
public class VCFManager<T> implements Runnable {

    static final String POISON_PILL = "__END__";

    private final List<String> inputFileNames;
    private int usingThreads;
    private final int maxSizeOfVariantCache;
    private final Function<String, T> variantParser;
    private boolean verbose = false;

    private BlockingQueue<T> variantRawCache;
    private Map<String, int[]> genotypeEncodingCache;

    private List<String> commentData;
    private String headerData;
    private AtomicInteger currVariantCount;

    private int numVariants;
    private int numSamples;
    private List<String> sampleNames;
    private int ploidy;
    private int maxAlleles;

    private CountDownLatch startSignal;
    private CountDownLatch doneSignal;
    private ExecutorService pool;
    private Map<Integer, VariantProcessor<T>> variantProcessors;

    public VCFManager(
            List<String> inputFileNames,
            int usingThreads,
            int maxSizeOfVariantCache,
            Function<String, T> variantParser,
            boolean verbose
    ) {
        this.inputFileNames = inputFileNames;
        this.usingThreads = usingThreads;
        this.maxSizeOfVariantCache = maxSizeOfVariantCache;
        this.variantParser = variantParser;
        this.verbose = verbose;
    }

    public void init() {
        if (inputFileNames==null || inputFileNames.isEmpty()) {
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
        doneSignal = new CountDownLatch(usingThreads + 1);

        pool = Executors.newFixedThreadPool(usingThreads);
        variantProcessors = new HashMap<>();
        VariantProcessor.resetCounters();
        for (int i = 0; i < usingThreads; i++) {
            VariantProcessor<T> vp = new VariantProcessor<>(
                    this,
                    startSignal,
                    doneSignal,
                    verbose);
            variantProcessors.put(vp.getId(), vp);
            pool.execute(vp);
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
            VCFStreamingIterator<String> iterator = new VCFStreamingIterator<>(decoder, verbose, inputFileNames);
            for (String line : iterator) {
                if(line != null) processVariantLine(line);
            }

            Logger.info(this, "END READ");
            numVariants = currVariantCount.get();

            // putting POISON
            for (int i = 0; i < usingThreads; i++) {
                putPoison((T) POISON_PILL);
            }

            doneSignal.countDown();
        } catch (InterruptedException e) {
            Logger.error(this, e.getMessage());
            shutdown();
            System.exit(1);
        }
        finally {
            shutdown();
        }
    }

    private void processVariantLine(String line) {
        try {
            if (line.startsWith("##")) {
                commentData.add(line);
            } else if (line.startsWith("#")) {
                headerData = line;
                sampleNames = getSampleNamesFromHeader();
                if(sampleNames==null || sampleNames.size()<=0) {
                    shutdown();
                    System.exit(1);
                    return;
                }
                numSamples = sampleNames.size();
            } else {
                currVariantCount.incrementAndGet();
                if (ploidy <= 0) {
                    try{
                        int[] ploidy_maxAlleles = SNPEncoder.guessPloidyAndMaxAllele(line);
                        ploidy = ploidy_maxAlleles[0];
                        maxAlleles = ploidy_maxAlleles[1];
                        if (ploidy > 0) {
                            startSignal.countDown();
                        }
                    } catch(IllegalArgumentException e){
                        Logger.error(this, e.getMessage());
                        shutdown();
                        System.exit(1);
                        return;
                    }
                }
                T parsed = variantParser.apply(line);
                putVariantRaw(parsed);
            }
        } catch (InterruptedException e) {
            Logger.error(this, e.getMessage());
        }
    }

    private void putVariantRaw(T variant) throws InterruptedException {
        variantRawCache.put(variant); 
    }

    public T getNextVariantRaw() throws InterruptedException {
        return variantRawCache.take();
    }

    public Map<String, int[]> getEncodingCache() {
        return genotypeEncodingCache;
    }

    private void putPoison(T poison) throws InterruptedException{
        variantRawCache.put(poison);
    }

    public float[][] reduceDotProd() {
        int[][] finalDotProd = new int[numSamples][numSamples];
        int[] finalNorm = new int[numSamples];
        
        for (VariantProcessor<T> vp : variantProcessors.values()) {
            int[][] partialDot = vp.getLocalDotProd();
            int[] partialNorm = vp.getLocalNorm();
            for (int i = 0; i < numSamples; i++) {
                finalNorm[i] += partialNorm[i];
                for (int j = i; j < numSamples; j++) {
                    finalDotProd[i][j] += partialDot[i][j];
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

    public void shutdown() {
        try {
            pool.shutdown();
            if (!pool.awaitTermination(10, TimeUnit.MILLISECONDS)) {
                Logger.error(this, "Forcing shutdown...");
                pool.shutdown();
            }
        } catch (InterruptedException ex) {
        }
    }

    @SuppressWarnings("unused")
    private void clear() {
        if (variantRawCache != null) {
            variantRawCache.clear();
        }
    }

    public final List<String> getCommentData() {
        return commentData;
    }

    public final String getHeaderData() {
        return headerData;
    }

    public List<String> getSampleNames() {
        return sampleNames;
    }

    public int getCurrVariantCount() {
        return currVariantCount.get();
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
