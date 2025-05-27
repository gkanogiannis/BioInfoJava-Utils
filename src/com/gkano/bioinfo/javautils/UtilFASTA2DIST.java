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
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.gkano.bioinfo.fasta.CalculateDistancesD2;
import com.gkano.bioinfo.fasta.FastaManager;
import com.gkano.bioinfo.fasta.SequenceD2;
import com.gkano.bioinfo.fasta.SequenceProcessor;
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

	@Parameter(description = "FASTA positional input files")
    private List<String> positionalInputFiles = new ArrayList<>();

	@Parameter(names = { "-i", "--input" }, description = "FASTA input file(s)", variableArity = true)
    private List<String> namedInputFiles = new ArrayList<>();

	@Parameter(names = { "-o", "--output" }, description = "Distance output file")
    private String outputFile;
	
	@Parameter(names = {"--isFastq","-q"}, description = "Input is FASTQ (default: false)")
	private boolean isFastq = false;

	@Parameter(names = {"--kmerSize","-k"}, description = "Kmer size (default: 4)")
	private Integer k = 4;

	@Parameter(names = { "--normalize", "-n" })
	private boolean normalize = false;
	
	@Parameter(names = { "--numberOfThreads", "-t" })
	private int numOfThreads = 1;

	@Parameter(names={"--useMappedBuffer"}, description="Use MappedByteBuffer for reading input files. Not compatible with piped input.")
	private boolean useMappedBuffer = false;

	@SuppressWarnings("unused")
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

			// Merge all FASTA inputs into one list
            List<String> inputFileNames = new ArrayList<>();
            inputFileNames.addAll(positionalInputFiles);
            inputFileNames.addAll(namedInputFiles);

			if (inputFileNames.isEmpty()) {
                System.err.println("Error: No FASTA input files provided.");
                return;
            }

			int cpus = Runtime.getRuntime().availableProcessors();
			int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
			if(verbose) System.err.println("cpus=" + cpus);
			if(verbose) System.err.println("using=" + usingThreads);

			ConcurrentHashMap<Integer, SequenceD2> seqVectors = new ConcurrentHashMap<>();

			CountDownLatch startSignal = new CountDownLatch(1);
			CountDownLatch doneSignal = new CountDownLatch(usingThreads + 1);

			//if(verbose) System.err.println(GeneralTools.time() + " START ");

			ExecutorService pool = Executors.newFixedThreadPool(usingThreads + 1);

			//Map<Integer, SequenceProcessor> sequenceProcessors = new HashMap<>();
			FastaManager frm = new FastaManager(isFastq, false, inputFileNames, startSignal, doneSignal, useMappedBuffer);
			pool.execute(frm);

			SequenceProcessor.resetCounters();
			// Starting threads
			for (int i = 0; i < usingThreads; i++) {
				SequenceProcessor sp = new SequenceProcessor(seqVectors, frm, k, normalize, startSignal, doneSignal, verbose);
				//sequenceProcessors.put(sp.getId(), sp);
				pool.execute(sp);
			}

			doneSignal.await();
			pool.shutdown();
				
			List<String> seqNames = frm.getSequenceNames();
			List<Integer> seqIds = frm.getSequenceIds();
			
			// Calculate distances
			CalculateDistancesD2 fj = new CalculateDistancesD2(verbose);
			double[][] distances = fj.calculateDistances(usingThreads, seqVectors, seqNames, seqIds, frm);
			
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
		catch (InterruptedException e) {
			Logger.error(this, e.getMessage());
		}
	}

}
