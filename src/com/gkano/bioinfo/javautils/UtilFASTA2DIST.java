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
import java.util.concurrent.ConcurrentHashMap;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.gkano.bioinfo.fasta2.DistanceCalculator;
import com.gkano.bioinfo.fasta2.FastaManager;
import com.gkano.bioinfo.fasta2.SequenceD2;
import com.gkano.bioinfo.var.GeneralTools;
import com.gkano.bioinfo.var.Logger;

@SuppressWarnings("FieldMayBeFinal")
@Parameters(commandDescription = "FASTA2DIST")
public class UtilFASTA2DIST {

	public UtilFASTA2DIST() {
	}

	public static String getUtilName() {
		return "FASTA2DIST";
	}

	@SuppressWarnings("unused")
	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(names = { "-v", "--verbose"})
	private boolean verbose = false;

	@Parameter(description = "<positional input files>")
    private List<String> positionalInputFiles = new ArrayList<>();

	@Parameter(names = { "-i", "--input" }, description = "FASTA input file(s)", variableArity = true)
    private List<String> namedInputFiles = new ArrayList<>();

	@Parameter(names = { "-o", "--output" }, description = "Distance output file")
    private String outputFile;
	
	@Parameter(names = {"--isFastq","-q"}, description = "Input is FASTQ (default: false)")
	private boolean isFastq = false;

	@Parameter(names = {"--kmerSize","-k"}, description = "Kmer size (default: 4)")
	private Integer k = 4;

	@Parameter(names = { "--normalize", "-n" }, description = "Normalize probabilities (default: false)")
	private boolean normalize = false;

	@Parameter(names = { "--numberOfThreads", "-t" }, description = "Number of threads (default: 1)")
	private int numOfThreads = 1;

	public void go() {
		try (PrintStream ops = GeneralTools.getPrintStreamOrExit(outputFile, this)) {

			// Merge all FASTA inputs into one list
            List<String> inputFileNames = new ArrayList<>();
            inputFileNames.addAll(positionalInputFiles);
            inputFileNames.addAll(namedInputFiles);

			// Create manager with integrated thread control
			var manager = new FastaManager.Builder(inputFileNames)
				.withProcessorThreads(numOfThreads)
				.withKmerSize(k)
				.withNormalization(normalize)
				.keepQualities(isFastq)
				.withVerbose(verbose)
				.build();

			// Initialize - this starts all threads internally
			manager.init();

			// Wait for completion
			manager.awaitCompletion();

			// Get results
			ConcurrentHashMap<Integer, SequenceD2> seqVectors = manager.getResults();
			List<String> seqNames = manager.getSequenceNames();
			List<Integer> seqIds = manager.getSequenceIds();

			// Calculate distances. Pass seqIds so that row/col order of the
			// returned matrix matches seqNames (both in parse order).
			var calc = new DistanceCalculator(numOfThreads);
			double[][] distances = calc.computeD2Distances(seqVectors, seqIds);
			
			// Print data
			int seqCounter = 0;
			ops.println(seqNames.size());
			for(int i=0; i<seqNames.size(); i++) {
				String seqName1 = seqNames.get(i);
				ops.print(seqName1);
				for (int j=0;j<seqNames.size();j++) {
					String seqName2 = seqNames.get(j);
					ops.print("\t"+GeneralTools.decimalFormat.format(distances[i][j]));
				}
				ops.println("");
				seqCounter++;
				//if(++sampleCounter % 10 == 0 && verbose) System.err.println(GeneralTools.time()+" Samples Processed : \t"+sampleCounter);
			}
			ops.flush();
			ops.close();
			//if(verbose) System.err.println(GeneralTools.time()+" Samples Processed : \t"+sampleCounter);
		} 
		catch (Exception e) {
			Logger.error(this, e.getMessage());
		}
	}

}
