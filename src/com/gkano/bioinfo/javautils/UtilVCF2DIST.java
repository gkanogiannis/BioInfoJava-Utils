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
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.gkano.bioinfo.var.GeneralTools;
import com.gkano.bioinfo.var.Logger;
import com.gkano.bioinfo.vcf.CalculateDistancesCOSINE;
import com.gkano.bioinfo.vcf.VCFManager;
import com.gkano.bioinfo.vcf.VariantManager;
import com.gkano.bioinfo.vcf.VariantProcessor;

@SuppressWarnings("FieldMayBeFinal")
@Parameters(commandDescription = "VCF2DIST")
public class UtilVCF2DIST {

	public UtilVCF2DIST() {
	}

	public static String getUtilName() {
		return "VCF2DIST";
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

	@Parameter(names = { "-o", "--output" }, description = "Distance output file")
    private String outputFile;

	@Parameter(names = { "--numberOfThreads", "-t" })
	private int numOfThreads = 1;
	
	@Parameter(names = { "--ignoreMissing", "-m" })
	private boolean ignoreMissing = false;
	
	@Parameter(names={"--onlyHets", "-h"})
	private boolean onlyHets = false;
	
	@Parameter(names={"--ignoreHets", "-g"})
	private boolean ignoreHets = false;

	@Parameter(names={"--useMappedBuffer"}, description="Use MappedByteBuffer for reading input files. Not compatible with piped input.")
	private boolean useMappedBuffer = false;

	public void go() {
		try {
			//Output PrintStream
			PrintStream ops = System.out;
			if(outputFile != null) {
				try {
					ops = new PrintStream(outputFile);
				} 
				catch (FileNotFoundException e) {
					Logger.error(this, "Cannot write to " + outputFile);
                	return;
				}
			}

			// Merge all VCF inputs into one list
            List<String> inputFileNames = new ArrayList<>();
            inputFileNames.addAll(positionalInputFiles);
            inputFileNames.addAll(namedInputFiles);

			if (inputFileNames.isEmpty()) {
				Logger.error(this, "No VCF input files provided.");
                return;
            }

			int cpus = Runtime.getRuntime().availableProcessors();
			int usingThreads = (cpus < (numOfThreads+0) ? cpus : (numOfThreads+0));
			if(verbose) System.err.println("cpus=" + cpus);
			if(verbose) System.err.println("using=" + usingThreads);

			CountDownLatch startSignal = new CountDownLatch(1);
			CountDownLatch doneSignal = new CountDownLatch(usingThreads + 1);

			//if(verbose) System.err.println(GeneralTools.time() + " START ");

			ExecutorService pool = Executors.newFixedThreadPool(usingThreads + 1);

			//Map<Integer, VariantProcessor> variantProcessors = new HashMap<>();
			VariantManager vm = new VariantManager();
			VCFManager vcfm = new VCFManager(vm, inputFileNames, startSignal, doneSignal, useMappedBuffer, verbose);
			pool.execute(vcfm);

			VariantProcessor.resetCounters();
			// Starting threads
			for (int i = 0; i < usingThreads; i++) {
				VariantProcessor vp = new VariantProcessor(vm, vcfm, startSignal, doneSignal, verbose);
				//variantProcessors.put(vp.getId(), vp);
				pool.execute(vp);
			}

			doneSignal.await();
			pool.shutdown();
			
			vm.populateSampleVariant();

			List<String> sampleNames = new ArrayList<>();
			for (int i = 9; i < vcfm.getHeaderData().length; i++) {
				sampleNames.add(new String(vcfm.getHeaderData()[i]));
			}

			if(verbose) System.err.printf("\rProcessed variants : \t%8d\n", VariantProcessor.getVariantCount().get());

			// Calculate distances
			//double[][] distances = SampleVariantTools.calculateDistances(samplesToVariantsData, sampleNames, vcfm, ignoreHets, onlyHets, false);
			CalculateDistancesCOSINE.resetCounters();
			CalculateDistancesCOSINE fj = new CalculateDistancesCOSINE(verbose);
			double[][] distances = fj.calculateDistances(usingThreads, sampleNames, vm, vcfm, ignoreHets, onlyHets, ignoreMissing);
			
			// Print data
			@SuppressWarnings("unused")
			int sampleCounter = 0;
			ops.println(vm.getNumSamples()+"\t"+vm.getNumVariants());
			for(int i=0; i<vm.getNumSamples(); i++) {
				String sampleName1 = sampleNames.get(i);
				ops.print(sampleName1);
				for (int j=0;j<vm.getNumSamples();j++) {
					@SuppressWarnings("unused")
					String sampleName2 = sampleNames.get(j);
					ops.print("\t"+GeneralTools.decimalFormat.format(distances[i][j]));
				}
				ops.println("");
				sampleCounter++;
				//if(++sampleCounter % 10 == 0 && verbose) System.err.println(GeneralTools.time()+" Samples Processed : \t"+sampleCounter);
			}
			ops.flush();
			ops.close();
			//if(verbose) System.err.println(GeneralTools.time()+" Samples Processed : \t"+sampleCounter);

		} 
		catch (InterruptedException e) {
			Logger.error(this, e.getMessage());
		}
	}
}
