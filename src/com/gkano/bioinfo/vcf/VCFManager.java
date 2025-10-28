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

    private int numBootstraps = 0; // Set this from the constructor or a setter

    public static class ProcessorResult {
        int[][][] dotProd; // [replicate][i][j]
        int[][] norm;      // [replicate][i]
        int numReplicates;
        private ProcessorResult(int numReplicates, int numSamples) {
            this.numReplicates = numReplicates;
            this.dotProd = new int[numReplicates][numSamples][numSamples];
            this.norm = new int[numReplicates][numSamples];
        }
        public void merge(ProcessorResult other) {
            for (int r = 0; r < numReplicates; r++) {
                for (int i = 0; i < norm[r].length; i++) {
                    norm[r][i] += other.norm[r][i];
                    for (int j = i; j < norm[r].length; j++) {
                        dotProd[r][i][j] += other.dotProd[r][i][j];
                    }
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
        Logger.setVerbose(verbose);
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
        Logger.info(this, "cpus=" + cpus + "\tusing=" + usingThreads);

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

                int numReplicates = numBootstraps > 0 ? numBootstraps + 1 : 1;
                ProcessorResult result = new ProcessorResult(numReplicates, numSamples);
                java.util.Random rand = new java.util.Random();

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
                        // For each replicate, decide if this SNP is included (with replacement)
                        // Replicate 0 is always the original (no resampling)
                        int[] replicateCounts = new int[numReplicates];
                        replicateCounts[0] = 1; // always include in original
                        for (int r = 1; r < numReplicates; r++) {
                            replicateCounts[r] = poisson1(rand); // 0 or more times
                        }
                        for (int r = 0; r < numReplicates; r++) {
                            if (replicateCounts[r] == 0) continue;
                            for (int i = 0; i < numSamples; i++) {
                                di = variantEncoded[i];
                                norm = 0;
                                for (int k = 0; k < di.length; k++) norm += Integer.bitCount(di[k] & di[k]);
                                result.norm[r][i] += norm * replicateCounts[r];
                                for (int j = i; j < numSamples; j++) {
                                    dj = variantEncoded[j];
                                    dotProd = 0;
                                    for (int k = 0; k < di.length; k++) dotProd += Integer.bitCount(di[k] & dj[k]);
                                    result.dotProd[r][i][j] += dotProd * replicateCounts[r];
                                }
                            }
                        }
                    }
                }
                return result;
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

            Logger.info(this, "Processed variants :\t" + numVariants);
            
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
        try {
            int[][] finalDotProd = new int[numSamples][numSamples];
            int[] finalNorm = new int[numSamples];
    
            for (CompletableFuture<ProcessorResult> vp : variantProcessors) {
                ProcessorResult r = vp.join();
                for (int i = 0; i < numSamples; i++) {
                    finalNorm[i] += r.norm[0][i];
                    for (int j = i; j < numSamples; j++) {
                        finalDotProd[i][j] += r.dotProd[0][i][j];
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
        } catch (Exception e) {
            Logger.error(this, e.getMessage());
            return null;
        }
    }

    // At reduction, for each replicate, sum across threads and convert to distance matrix as before.
    public List<float[][]> reduceDotProdToDistancesBootstraps() {
        try {
            int numReplicates = numBootstraps > 0 ? numBootstraps + 1 : 1;
            int[][][] finalDotProd = new int[numReplicates][numSamples][numSamples];
            int[][] finalNorm = new int[numReplicates][numSamples];
    
            for (CompletableFuture<ProcessorResult> vp : variantProcessors) {
                ProcessorResult r = vp.join();
                for (int rep = 0; rep < numReplicates; rep++) {
                    for (int i = 0; i < numSamples; i++) {
                        finalNorm[rep][i] += r.norm[rep][i];
                        for (int j = i; j < numSamples; j++) {
                            finalDotProd[rep][i][j] += r.dotProd[rep][i][j];
                        }
                    }
                }
            }
    
            List<float[][]> allDistances = new ArrayList<>();
            for (int rep = 0; rep < numReplicates; rep++) {
                float[][] cosineDist = new float[numSamples][numSamples];
                for (int i = 0; i < numSamples; i++) {
                    float normI = (float) Math.sqrt(finalNorm[rep][i]);
                    for (int j = i; j < numSamples; j++) {
                        float normJ = (float) Math.sqrt(finalNorm[rep][j]);
                        float dot = finalDotProd[rep][i][j];
                        float similarity = ((normI > 0 && normJ > 0) ? (dot / (normI * normJ)) : 0.0f);
                        if (j == i && normI == 0 && normJ == 0) {
                            similarity = 1.0f;
                        }
                        float dist = 1.0f - similarity;
                        if(dist<0) dist = 0f;
                        cosineDist[i][j] = cosineDist[j][i] = dist;
                    }
                }
                allDistances.add(cosineDist);
            }
            return allDistances;
        } catch (OutOfMemoryError e) {
            Logger.error(this, "Could not allocate memory for bootstrap distance matrices: use --mem option to increase memory available to JVM.");
            System.exit(256);
            return null;
        } 
        catch (Exception e) {
            Logger.error(this, e.getMessage());
            return null;
        }
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

    private static int poisson1(java.util.Random rand) {
        double L = Math.exp(-1.0);
        int k = 0;
        double p = 1.0;
        do {
            k++;
            p *= rand.nextDouble();
        } while (p > L);
        return k - 1;
    }
}
