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
public class VCFManager implements Runnable {

    //private static final String POISON_PILL = "__END__";
    private final List<String> POISON_PILL_BATCH = Collections.emptyList();
    private final int batchSize = 1000;

    private final List<String> inputFileNames;
    private int usingThreads;
    private final int maxSizeOfVariantCache;
    private final Function<String, String> variantParser;
    private boolean verbose = false;

    private BlockingQueue<List<String>> variantRawCache;
    private Map<String, int[]> genotypeEncodingCache;

    private List<String> commentData;
    private String headerData;
    private AtomicInteger currVariantCount;
    private AtomicInteger skippedVariantCount;

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
    private List<CompletableFuture<EmbeddingProcessorResult>> embeddingProcessors;

    private int numBootstraps = 0; // Set this from the constructor or a setter

    // Embedding mode fields
    private boolean embeddingMode = false;
    private Map<String, double[]> variantEmbeddings;
    private int embeddingDim = -1;
    private VariantKeyExtractor keyExtractor;

    /**
     * Result container for standard genotype-based processing.
     */
    public static class ProcessorResult {
        long[][][] dotProd; // [replicate][i][j]
        long[][] norm;      // [replicate][i]
        int numReplicates;
        private ProcessorResult(int numReplicates, int numSamples) {
            this.numReplicates = numReplicates;
            this.dotProd = new long[numReplicates][numSamples][numSamples];
            this.norm = new long[numReplicates][numSamples];
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

    /**
     * Result container for embedding-based processing.
     * Stores accumulated sample embeddings in embedding space.
     */
    public static class EmbeddingProcessorResult {
        double[][][] sampleEmbeddings; // [replicate][sample][embeddingDim]
        int numReplicates;
        int embeddingDim;

        private EmbeddingProcessorResult(int numReplicates, int numSamples, int embeddingDim) {
            this.numReplicates = numReplicates;
            this.embeddingDim = embeddingDim;
            this.sampleEmbeddings = new double[numReplicates][numSamples][embeddingDim];
        }

        public void merge(EmbeddingProcessorResult other) {
            for (int r = 0; r < numReplicates; r++) {
                for (int i = 0; i < sampleEmbeddings[r].length; i++) {
                    for (int d = 0; d < embeddingDim; d++) {
                        sampleEmbeddings[r][i][d] += other.sampleEmbeddings[r][i][d];
                    }
                }
            }
        }
    }

    public VCFManager(
            List<String> inputFileNames,
            int usingThreads,
            Function<String, String> variantParser,
            boolean verbose
    ) {
        this.inputFileNames = inputFileNames;
        this.usingThreads = usingThreads;
        this.variantParser = variantParser;
        this.verbose = verbose;
        this.maxSizeOfVariantCache = 4 * usingThreads;
    }

    /**
     * Enable embedding mode with pre-loaded embeddings.
     *
     * @param embeddings Map from variant key (CHROM:POS:REF:ALT) to embedding vector
     * @param keyFormat  Format for extracting variant keys from VCF lines
     */
    public void setEmbeddings(Map<String, double[]> embeddings, VariantKeyExtractor.KeyFormat keyFormat) {
        if (embeddings == null || embeddings.isEmpty()) {
            this.embeddingMode = false;
            this.variantEmbeddings = null;
            this.embeddingDim = -1;
            this.keyExtractor = null;
            return;
        }
        this.embeddingMode = true;
        this.variantEmbeddings = embeddings;
        this.embeddingDim = VariantEmbeddingLoader.getEmbeddingDimension(embeddings);
        this.keyExtractor = new VariantKeyExtractor(keyFormat);
        Logger.info(this, "Embedding mode enabled: " + embeddings.size() +
                    " embeddings with dimension " + embeddingDim);
    }

    /**
     * Check if embedding mode is enabled.
     */
    public boolean isEmbeddingMode() {
        return embeddingMode;
    }

    /**
     * Get embedding dimension (only valid in embedding mode).
     */
    public int getEmbeddingDim() {
        return embeddingDim;
    }

    public void init() {
        Logger.setVerbose(verbose);
        if (inputFileNames == null || inputFileNames.isEmpty()) {
            Logger.error(this, "No VCF input files provided.");
            throw new IllegalArgumentException("No VCF input files provided.");
        }

        variantRawCache = new LinkedBlockingQueue<>(maxSizeOfVariantCache);
        genotypeEncodingCache = new ConcurrentHashMap<>();
        commentData = new ArrayList<>();
        currVariantCount = new AtomicInteger(0);
        skippedVariantCount = new AtomicInteger(0);

        int cpus = Runtime.getRuntime().availableProcessors();
        usingThreads = (cpus < usingThreads ? cpus : usingThreads);
        Logger.info(this, "cpus=" + cpus + "\tusing=" + usingThreads);

        startSignal = new CountDownLatch(1);
        doneSignal = new CountDownLatch(1);

        pool = Executors.newFixedThreadPool(usingThreads);

        if (embeddingMode) {
            initEmbeddingProcessors();
        } else {
            initStandardProcessors();
        }
    }

    /**
     * Initialize worker threads for standard genotype-based processing.
     */
    private void initStandardProcessors() {
        variantProcessors = new ArrayList<>();

        for (int t = 0; t < usingThreads; t++) {
            CompletableFuture<ProcessorResult> variantProcessor = CompletableFuture.supplyAsync(() -> {
                try {
                    // Wait until ploidy is known
                    if (!startSignal.await(5, TimeUnit.MINUTES)) {
                        throw new IllegalStateException("Timeout waiting for ploidy/maxAlleles initialization");
                    }
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    throw new RuntimeException(e);
                }

                int numReplicates = numBootstraps > 0 ? numBootstraps + 1 : 1;
                ProcessorResult result = new ProcessorResult(numReplicates, numSamples);
                java.util.Random rand = new java.util.Random();

                List<String> batch;
                int[][] variantEncoded;
                int count, step;
                int[] di, dj;
                long norm, dotProd;
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
                    for (String variant : batch) {
                        count = currVariantCount.incrementAndGet();
                        if(verbose){
                            step = GeneralTools.getAdaptiveVariantStep(count);
                            if (count % step == 0) {
                                Logger.infoCarret(this, "Variants:\t" + count);
                            }
                        }
                        try {
                            variantEncoded = SNPEncoder.encodeSNPOneHot(
                                variant, ploidy, maxAlleles, genotypeEncodingCache, numSamples);
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
                                for (int k = 0; k < di.length; k++) {
                                    norm += Integer.bitCount(di[k]);
                                }
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

    /**
     * Initialize worker threads for embedding-based processing.
     * Uses kernel trick: K[v,w] = e_v · e_w decomposed as h_i = Σ_v contrib_v × e_v
     */
    private void initEmbeddingProcessors() {
        embeddingProcessors = new ArrayList<>();

        for (int t = 0; t < usingThreads; t++) {
            CompletableFuture<EmbeddingProcessorResult> embeddingProcessor = CompletableFuture.supplyAsync(() -> {
                try {
                    // Wait until ploidy is known
                    if (!startSignal.await(5, TimeUnit.MINUTES)) {
                        throw new IllegalStateException("Timeout waiting for ploidy/maxAlleles initialization");
                    }
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    throw new RuntimeException(e);
                }

                int numReplicates = numBootstraps > 0 ? numBootstraps + 1 : 1;
                EmbeddingProcessorResult result = new EmbeddingProcessorResult(numReplicates, numSamples, embeddingDim);
                java.util.Random rand = new java.util.Random();

                List<String> batch;
                int count, step;

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

                    for (String variant : batch) {
                        count = currVariantCount.incrementAndGet();
                        if (verbose) {
                            step = GeneralTools.getAdaptiveVariantStep(count);
                            if (count % step == 0) {
                                Logger.infoCarret(this, "Variants:\t" + count);
                            }
                        }

                        // Extract variant key and lookup embedding
                        String variantKey = keyExtractor.extractKey(variant);
                        if (variantKey == null) {
                            skippedVariantCount.incrementAndGet();
                            continue;
                        }

                        double[] embedding = variantEmbeddings.get(variantKey);
                        if (embedding == null) {
                            // Skip variants without embeddings
                            skippedVariantCount.incrementAndGet();
                            continue;
                        }

                        // Compute allele dosage for each sample
                        // Dosage = count of alternate alleles (0, 1, 2, ... or -1 for missing)
                        int[] alleleDosage;
                        try {
                            alleleDosage = SNPEncoder.computeAlleleDosage(variant, numSamples);
                        } catch (IllegalArgumentException e) {
                            skippedVariantCount.incrementAndGet();
                            continue;
                        }

                        // Bootstrap resampling
                        int[] replicateCounts = new int[numReplicates];
                        replicateCounts[0] = 1; // always include in original
                        for (int r = 1; r < numReplicates; r++) {
                            replicateCounts[r] = poisson1(rand);
                        }

                        // Accumulate embeddings for each sample
                        // h_i[d] += dosage × replicateCount × embedding[d]
                        // where dosage = count of alternate alleles (0, 1, 2, ...)
                        for (int r = 0; r < numReplicates; r++) {
                            if (replicateCounts[r] == 0) continue;

                            for (int i = 0; i < numSamples; i++) {
                                int dosage = alleleDosage[i];

                                if (dosage <= 0) continue; // Missing or reference-only genotype

                                // Weight = dosage × replicateCount
                                double weight = dosage * replicateCounts[r];

                                // Accumulate weighted embedding
                                for (int d = 0; d < embeddingDim; d++) {
                                    result.sampleEmbeddings[r][i][d] += weight * embedding[d];
                                }
                            }
                        }
                    }
                }
                return result;
            }, pool);
            embeddingProcessors.add(embeddingProcessor);
        }
    }

    public void awaitFinalization() throws InterruptedException {
        doneSignal.await();
    }

    @Override
    public void run() {
        try {
            Logger.info(this, "START READ" + (embeddingMode ? " (embedding mode)" : ""));

            VCFDecoder decoder = new VCFDecoder();
            List<String> batch = new ArrayList<>(2500);
            try (VCFStreamingIterator iterator = new VCFStreamingIterator(decoder, verbose, inputFileNames)) {
                for (String line : iterator) {
                    if (line != null) processVariantLine(line, batch);
                }
            }
            // Push the last partial batch if needed
            if (!batch.isEmpty()) variantRawCache.put(batch);
            // If ploidy was never inferred (no data lines), release workers to exit cleanly
            if (ploidy <= 0) {
                startSignal.countDown();
            }
            // putting POISON pills
            for (int i = 0; i < usingThreads; i++) variantRawCache.put(POISON_PILL_BATCH);

            Logger.info(this, "END READ");

            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
            numVariants = currVariantCount.get();

            Logger.info(this, "Processed variants :\t" + numVariants);
            if (embeddingMode) {
                int skipped = skippedVariantCount.get();
                int used = numVariants - skipped;
                Logger.info(this, "Variants with embeddings:\t" + used +
                           " (skipped " + skipped + " without embeddings)");
            }

            doneSignal.countDown();

        } catch (Exception e) {
            Logger.error(this, e.getMessage());
            shutdown();
            doneSignal.countDown();
            Thread.currentThread().interrupt();
        } finally {
            shutdown();
        }
    }

    private void processVariantLine(String line, List<String> batch) throws Exception {
        try {
            if (line.startsWith("##")) {
                commentData.add(line);
            } else if (line.startsWith("#")) {
                headerData = line;
                sampleNames = getSampleNamesFromHeader();
                if (sampleNames == null || sampleNames.size() <= 0) {
                    throw new IllegalStateException("No samples detected from #CHROM header.");
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
                        throw new IllegalStateException("Failed to infer ploidy / maxAlleles: " + e.getMessage(), e);
                    }
                }
                String parsed = variantParser.apply(line);
                batch.add(parsed);
                if (batch.size() >= batchSize) {
                    variantRawCache.put(new ArrayList<>(batch));
                    batch.clear();
                }
            }
        } catch (IllegalStateException | InterruptedException | IllegalArgumentException e) {
            throw  new RuntimeException(e);
        }
    }

    /**
     * Reduce standard dot products to distance matrix.
     */
    public double[][] reduceDotProdToDistances() {
        if (embeddingMode) {
            return reduceEmbeddingsToDistances();
        }
        try {
            long[][] finalDotProd = new long[numSamples][numSamples];
            long[] finalNorm = new long[numSamples];

            for (CompletableFuture<ProcessorResult> vp : variantProcessors) {
                ProcessorResult r = vp.join();
                for (int i = 0; i < numSamples; i++) {
                    finalNorm[i] += r.norm[0][i];
                    for (int j = i; j < numSamples; j++) {
                        finalDotProd[i][j] += r.dotProd[0][i][j];
                    }
                }
            }

            double[][] cosineDist = new double[numSamples][numSamples];
            for (int i = 0; i < numSamples; i++) {
                double normI = Math.sqrt(finalNorm[i]);
                for (int j = i; j < numSamples; j++) {
                    double normJ = Math.sqrt(finalNorm[j]);
                    double dot = finalDotProd[i][j];
                    double similarity = ((normI > 0 && normJ > 0) ? (dot / (normI * normJ)) : 0.0f);
                    if (j == i && normI == 0 && normJ == 0) {
                        similarity = 1.0;
                    }
                    double dist = 1.0 - similarity;
                    if(dist<0) dist = 0.0;
                    cosineDist[i][j] = cosineDist[j][i] = dist;
                }
            }
            return cosineDist;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Reduce embedding accumulations to distance matrix.
     * Computes cosine distance in embedding space.
     */
    public double[][] reduceEmbeddingsToDistances() {
        try {
            // Merge sample embeddings from all worker threads
            double[][] finalEmbeddings = new double[numSamples][embeddingDim];

            for (CompletableFuture<EmbeddingProcessorResult> ep : embeddingProcessors) {
                EmbeddingProcessorResult r = ep.join();
                for (int i = 0; i < numSamples; i++) {
                    for (int d = 0; d < embeddingDim; d++) {
                        finalEmbeddings[i][d] += r.sampleEmbeddings[0][i][d];
                    }
                }
            }

            // Compute pairwise cosine distances
            return computeCosineDistancesFromEmbeddings(finalEmbeddings);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Reduce with bootstrap replicates for standard mode.
     */
    public List<double[][]> reduceDotProdToDistancesBootstraps() {
        if (embeddingMode) {
            return reduceEmbeddingsToDistancesBootstraps();
        }
        try {
            int numReplicates = numBootstraps > 0 ? numBootstraps + 1 : 1;
            long[][][] finalDotProd = new long[numReplicates][numSamples][numSamples];
            long[][] finalNorm = new long[numReplicates][numSamples];

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

            List<double[][]> allDistances = new ArrayList<>();
            for (int rep = 0; rep < numReplicates; rep++) {
                double[][] cosineDist = new double[numSamples][numSamples];
                for (int i = 0; i < numSamples; i++) {
                    double normI = Math.sqrt(finalNorm[rep][i]);
                    for (int j = i; j < numSamples; j++) {
                        double normJ = Math.sqrt(finalNorm[rep][j]);
                        double dot = finalDotProd[rep][i][j];
                        double similarity = ((normI > 0 && normJ > 0) ? (dot / (normI * normJ)) : 0.0);
                        if (j == i && normI == 0 && normJ == 0) {
                            similarity = 1.0;
                        }
                        double dist = 1.0 - similarity;
                        if(dist<0) dist = 0.0;
                        cosineDist[i][j] = cosineDist[j][i] = dist;
                    }
                }
                allDistances.add(cosineDist);
            }
            return allDistances;
        } catch (OutOfMemoryError e) {
            Logger.error(this, "Could not allocate memory for bootstrap distance matrices: use --mem option to increase memory available to JVM.");
            throw new RuntimeException("Could not allocate memory for bootstrap distance matrices: " + e.getMessage(), e);
        }
        catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Reduce with bootstrap replicates for embedding mode.
     */
    public List<double[][]> reduceEmbeddingsToDistancesBootstraps() {
        try {
            int numReplicates = numBootstraps > 0 ? numBootstraps + 1 : 1;

            // Merge sample embeddings from all worker threads for all replicates
            double[][][] finalEmbeddings = new double[numReplicates][numSamples][embeddingDim];

            for (CompletableFuture<EmbeddingProcessorResult> ep : embeddingProcessors) {
                EmbeddingProcessorResult r = ep.join();
                for (int rep = 0; rep < numReplicates; rep++) {
                    for (int i = 0; i < numSamples; i++) {
                        for (int d = 0; d < embeddingDim; d++) {
                            finalEmbeddings[rep][i][d] += r.sampleEmbeddings[rep][i][d];
                        }
                    }
                }
            }

            // Compute pairwise cosine distances for each replicate
            List<double[][]> allDistances = new ArrayList<>();
            for (int rep = 0; rep < numReplicates; rep++) {
                double[][] cosineDist = computeCosineDistancesFromEmbeddings(finalEmbeddings[rep]);
                allDistances.add(cosineDist);
            }

            return allDistances;
        } catch (OutOfMemoryError e) {
            Logger.error(this, "Could not allocate memory for bootstrap distance matrices: use --mem option to increase memory available to JVM.");
            throw new RuntimeException("Could not allocate memory for bootstrap distance matrices: " + e.getMessage(), e);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Compute pairwise cosine distances from sample embeddings.
     *
     * @param embeddings Sample embeddings [numSamples][embeddingDim]
     * @return Distance matrix [numSamples][numSamples]
     */
    private double[][] computeCosineDistancesFromEmbeddings(double[][] embeddings) {
        int n = embeddings.length;
        double[][] cosineDist = new double[n][n];

        // Pre-compute norms
        double[] norms = new double[n];
        for (int i = 0; i < n; i++) {
            double sumSq = 0.0;
            for (int d = 0; d < embeddings[i].length; d++) {
                sumSq += embeddings[i][d] * embeddings[i][d];
            }
            norms[i] = Math.sqrt(sumSq);
        }

        // Compute pairwise cosine distances
        for (int i = 0; i < n; i++) {
            double normI = norms[i];
            for (int j = i; j < n; j++) {
                double normJ = norms[j];

                // Compute dot product
                double dot = 0.0;
                for (int d = 0; d < embeddings[i].length; d++) {
                    dot += embeddings[i][d] * embeddings[j][d];
                }

                double similarity;
                if (normI > 0 && normJ > 0) {
                    similarity = dot / (normI * normJ);
                } else if (i == j) {
                    similarity = 1.0; // Same sample with zero embedding
                } else {
                    similarity = 0.0;
                }

                double dist = 1.0 - similarity;
                if (dist < 0) dist = 0.0;
                if (dist > 2.0) dist = 2.0; // Cosine distance bounds

                cosineDist[i][j] = cosineDist[j][i] = dist;
            }
        }

        return cosineDist;
    }

    public List<String> getSampleNamesFromHeader() {
        if (headerData == null || !headerData.startsWith("#CHROM")) {
            throw new IllegalArgumentException("VCF header missing #CHROM line: " + headerData);
        }

        final String[] fields = headerData.split("\t", -1);
        if (fields.length <= 9) {
            return Collections.emptyList();
        }

        return Arrays.asList(Arrays.copyOfRange(fields, 9, fields.length));
    }

    private void shutdown() {
        pool.shutdownNow();
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

    public void setNumBootstraps(int replicates) {
        if (replicates < 0) throw new IllegalArgumentException("bootstrap replicates must be >= 0");
        this.numBootstraps = replicates;
    }

    /**
     * Get number of variants skipped (only meaningful in embedding mode).
     */
    public int getSkippedVariants() {
        return skippedVariantCount != null ? skippedVariantCount.get() : 0;
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
