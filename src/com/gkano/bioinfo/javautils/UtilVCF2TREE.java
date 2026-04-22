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
package com.gkano.bioinfo.javautils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.gkano.bioinfo.tree.Clade;
import com.gkano.bioinfo.tree.HierarchicalCluster;
import com.gkano.bioinfo.var.GeneralTools;
import com.gkano.bioinfo.var.Logger;
import com.gkano.bioinfo.vcf.SNPEncoder;
import com.gkano.bioinfo.vcf.VCFManager;
import com.gkano.bioinfo.vcf.VariantEmbeddingLoader;
import com.gkano.bioinfo.vcf.VariantKeyExtractor;
import com.gkano.bioinfo.vcf.WindowPolicy;
import com.gkano.bioinfo.vcf.WindowedDistanceWriter;

@SuppressWarnings("FieldMayBeFinal")
@Parameters(commandDescription = "VCF2TREE")
public class UtilVCF2TREE {

    public UtilVCF2TREE() {
    }

    public static String getUtilName() {
        return "VCF2TREE";
    }

    @SuppressWarnings("unused")
    @Parameter(names = "--help", help = true)
    private boolean help;

    @Parameter(names = {"-v", "--verbose"})
    private boolean verbose = false;

    @Parameter(description = "<positional input files>")
    private List<String> positionalInputFiles = new ArrayList<>();

    @Parameter(names = {"-i", "--input"}, description = "VCF input file(s)", variableArity = true)
    private List<String> namedInputFiles = new ArrayList<>();

    @Parameter(names = {"-o", "--output"}, description = "Tree output file")
    private String outputFile;

    @Parameter(names = {"--numberOfThreads", "-t"})
    private int numOfThreads = 1;

    @Parameter(names = {"--bootstrap", "-b"}, description = "Number of bootstrap replicates")
    private int numBootstraps = 0;

    // Embedding options
    @Parameter(names = {"-e", "--embeddings"},
               description = "Variant embeddings file for embedding-based distance calculation")
    private String embeddingsFile;

    @Parameter(names = {"--embeddings-format"},
               description = "Embeddings file format: TSV or HUGGINGFACE (auto-detected if not specified)")
    private String embeddingsFormat;

    @Parameter(names = {"--variant-key"},
               description = "Variant key format for embedding lookup: CHROM_POS, CHROM_POS_REF_ALT, or VCF_ID")
    private String variantKeyFormat = "CHROM_POS_REF_ALT";

    @Parameter(names = {"--window-bp"},
               description = "Emit one tree per genomic window of N base pairs (mutually exclusive with --window-variants)")
    private Integer windowBp;

    @Parameter(names = {"--window-variants"},
               description = "Emit one tree per N consecutive variants (mutually exclusive with --window-bp)")
    private Integer windowVariants;

    @Parameter(names = {"--step"},
               description = "Window step (defaults to window size, i.e. tiled). Sliding windows not yet implemented.")
    private Integer windowStep;

    @Parameter(names = {"--min-variants"},
               description = "Minimum number of variants required to emit a window (default 1)")
    private int windowMinVariants = 1;

    //@SuppressWarnings("unused")
    //@Parameter(names = { "--ignoreMissing", "-m" })
    //private boolean ignoreMissing = false;

    //@SuppressWarnings("unused")
    //@Parameter(names={"--onlyHets", "-h"})
    //private boolean onlyHets = false;

    //@SuppressWarnings("unused")
    //@Parameter(names={"--ignoreHets", "-g"})
    //private boolean ignoreHets = false;

    public void go() {
        try (PrintStream ops = GeneralTools.getPrintStreamOrExit(outputFile, this)) {
            VCFManager vcfm = new VCFManager(
                    Stream.concat(positionalInputFiles.stream(), namedInputFiles.stream()).collect(Collectors.toList()),
                    numOfThreads,
                    SNPEncoder.StringToStringParser,
                    verbose);

            // Set number of bootstraps
            vcfm.setNumBootstraps(numBootstraps);

            WindowPolicy windowPolicy = buildWindowPolicy();
            if (windowPolicy != null) {
                if (numBootstraps > 0) {
                    Logger.error(this, "bootstrap (-b) is not yet supported with windowed output");
                    return;
                }
                if (embeddingsFile != null && !embeddingsFile.isEmpty()) {
                    Logger.error(this, "embeddings (-e) are not yet supported with windowed output");
                    return;
                }
                final boolean v = verbose;
                WindowedDistanceWriter writer = (chrom, start, end, n, names, d) -> {
                    ops.println("# window chrom=" + chrom + " start=" + start + " end=" + end
                            + " nvariants=" + n + " nsamples=" + names.length);
                    try {
                        ops.println((String) new HierarchicalCluster(v)
                                .hclusteringTree(names, d, null)[0]);
                    } catch (Exception e) {
                        Logger.error(this, "Failed to build per-window tree: " + e.getMessage());
                    }
                };
                vcfm.setWindowing(windowPolicy, writer);
                vcfm.init();
                new Thread(vcfm).start();
                vcfm.awaitFinalization();
                ops.flush();
                return;
            }

            // Load and set embeddings if provided
            if (embeddingsFile != null && !embeddingsFile.isEmpty()) {
                Logger.info(this, "Loading embeddings from: " + embeddingsFile);
                Map<String, double[]> embeddings;
                if (embeddingsFormat != null && !embeddingsFormat.isEmpty()) {
                    VariantEmbeddingLoader.EmbeddingFormat format =
                        VariantEmbeddingLoader.EmbeddingFormat.valueOf(embeddingsFormat.toUpperCase());
                    embeddings = VariantEmbeddingLoader.loadEmbeddings(embeddingsFile, format);
                } else {
                    embeddings = VariantEmbeddingLoader.loadEmbeddings(embeddingsFile);
                }
                VariantKeyExtractor.KeyFormat keyFormat = VariantKeyExtractor.parseFormat(variantKeyFormat);
                vcfm.setEmbeddings(embeddings, keyFormat);
            }

            vcfm.init();
            new Thread(vcfm).start();
            vcfm.awaitFinalization();

            List<String> sampleNames = vcfm.getSampleNames();
            if (sampleNames.size() < 2) {
                Logger.error(this, "At least two samples are required to build a tree.");
                return;
            }

            if (numBootstraps > 0) {
                List<double[][]> allDistances = vcfm.reduceDotProdToDistancesBootstraps();
                HierarchicalCluster hc = new HierarchicalCluster(verbose);
                Logger.info(this, allDistances.size() + " distance matrices computed (1 original + " + (allDistances.size() - 1) + " bootstraps).");

                if (vcfm.isEmbeddingMode()) {
                    int used = vcfm.getNumVariants() - vcfm.getSkippedVariants();
                    Logger.info(this, "Used " + used + " variants with embeddings (skipped " + vcfm.getSkippedVariants() + ")");
                }

                // Original tree
                Object[] originalTreeAndRoot = hc.hclusteringTree(sampleNames.toArray(String[]::new), allDistances.get(0), null);
                String originalTree = (String) originalTreeAndRoot[0];
                Clade originalRoot = (Clade) originalTreeAndRoot[1];

                // Collect bootstrap trees
                List<Clade> bootstrapTrees = new ArrayList<>();
                for (int b = 1; b < allDistances.size(); b++) {
                    Object[] bootTreeAndRoot = hc.hclusteringTree(sampleNames.toArray(String[]::new), allDistances.get(b), null);
                    String bootTree = (String) bootTreeAndRoot[0];
                    Clade bootRoot = (Clade) bootTreeAndRoot[1];
                    bootstrapTrees.add(bootRoot);
                }

                // Summarize bootstrap support and annotate tree
                Clade.annotateTreeWithBootstrap(originalRoot, bootstrapTrees, sampleNames.size(), numBootstraps);
                ops.println(Clade.cladeToString(originalRoot));
            } else {
                double[][] distances = vcfm.reduceDotProdToDistances();

                if (vcfm.isEmbeddingMode()) {
                    int used = vcfm.getNumVariants() - vcfm.getSkippedVariants();
                    Logger.info(this, "Used " + used + " variants with embeddings (skipped " + vcfm.getSkippedVariants() + ")");
                }

                HierarchicalCluster hc = new HierarchicalCluster(verbose);
                String treeString = (String) hc.hclusteringTree(sampleNames.toArray(String[]::new), distances, null)[0];
                ops.println(treeString);
            }
            ops.close();
        } catch (Exception e) {
            Logger.error(this, e.getMessage());
        }
    }

    private WindowPolicy buildWindowPolicy() {
        if (windowBp != null && windowVariants != null) {
            throw new IllegalArgumentException("--window-bp and --window-variants are mutually exclusive");
        }
        if (windowBp == null && windowVariants == null) {
            return null;
        }
        WindowPolicy.Mode mode = (windowBp != null) ? WindowPolicy.Mode.BP : WindowPolicy.Mode.VARIANTS;
        int size = (windowBp != null) ? windowBp : windowVariants;
        int step = (windowStep != null) ? windowStep : size;
        return new WindowPolicy(mode, size, step, windowMinVariants);
    }
}
