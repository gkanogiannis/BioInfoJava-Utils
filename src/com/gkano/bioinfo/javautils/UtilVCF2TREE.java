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

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.gkano.bioinfo.tree.HierarchicalCluster;
import com.gkano.bioinfo.var.Logger;
import com.gkano.bioinfo.vcf.SNPEncoder;
import com.gkano.bioinfo.vcf.VCFManager;
import com.gkano.bioinfo.vcf.VariantManager;
import com.gkano.bioinfo.vcf.VariantProcessor;

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

	@Parameter(names = { "-v", "--verbose"})
	private boolean verbose = false;

	@Parameter(description = "VCF positional input files")
    private List<String> positionalInputFiles = new ArrayList<>();

	@Parameter(names = { "-i", "--input" }, description = "VCF input file(s)", variableArity = true)
    private List<String> namedInputFiles = new ArrayList<>();
	
	@Parameter(names = { "-o", "--output" }, description = "Tree output file")
    private String outputFile;

	@Parameter(names = { "--numberOfThreads", "-t" })
	private int numOfThreads = 1;
	
	@SuppressWarnings("unused")
	@Parameter(names = { "--ignoreMissing", "-m" })
	private boolean ignoreMissing = false;
	
	@SuppressWarnings("unused")
	@Parameter(names={"--onlyHets", "-h"})
	private boolean onlyHets = false;
	
	@SuppressWarnings("unused")
	@Parameter(names={"--ignoreHets", "-g"})
	private boolean ignoreHets = false;

	public void go() {
		try {
			//Output PrintStream
			PrintStream ops = System.out;
			if(outputFile != null) {
				try {
					ops = new PrintStream(outputFile);
				} 
				catch (FileNotFoundException e) {
					System.err.println("Error: Cannot write to " + outputFile);
                	return;
				}
			}

			// Merge all VCF inputs into one list
            List<String> inputFileNames = new ArrayList<>();
            inputFileNames.addAll(positionalInputFiles);
            inputFileNames.addAll(namedInputFiles);

			if (inputFileNames.isEmpty()) {
                System.err.println("Error: No VCF input files provided.");
                return;
            }

			int cpus = Runtime.getRuntime().availableProcessors();
			int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
			if(verbose) System.err.println("cpus=" + cpus);
			if(verbose) System.err.println("using=" + usingThreads);

			CountDownLatch startSignal = new CountDownLatch(1);
			CountDownLatch doneSignal = new CountDownLatch(usingThreads + 1);

			//if(verbose) System.err.println(GeneralTools.time() + " START ");

			ExecutorService pool = Executors.newFixedThreadPool(usingThreads + 1);

			Map<Integer, VariantProcessor<String>> variantProcessors = new HashMap<>();
			VariantManager<String> vm = new VariantManager<>(10*usingThreads);
			
			VCFManager<String> vcfm = new VCFManager<>(vm, SNPEncoder.StringToStringParser, inputFileNames, startSignal, doneSignal, verbose);
            pool.execute(vcfm);

			VariantProcessor.resetCounters();
			// Starting threads
			for (int i = 0; i < usingThreads; i++) {
				VariantProcessor<String> vp = new VariantProcessor<>(vm, vcfm, startSignal, doneSignal, verbose);
				variantProcessors.put(vp.getId(), vp);
				pool.execute(vp);
			}

			doneSignal.await();
			pool.shutdown();
			
			List<String> sampleNames = SNPEncoder.getSampleNamesFromHeader(vcfm.getHeaderData());
			
			if(verbose) System.err.printf("\rProcessed variants : \t%8d\n", VariantProcessor.getVariantCount().get());

			// Calculate distances
            //double[][] distances = vm.calculateCosineDistances();
			float[][] distances = vm.reduce(variantProcessors.values());

			//HCluster tree
			HierarchicalCluster hc = new HierarchicalCluster(verbose);
			String treeString = hc.hclusteringTree(sampleNames.toArray(new String[0]), distances, null);
			ops.println(treeString);
			ops.flush();
			ops.close();
		} 
		catch (InterruptedException e) {
			Logger.error(this, e.getMessage());
		}
	}
	
}
