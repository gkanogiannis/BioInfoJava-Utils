/*
 *
 * BioInfoJava-Utils ciat.agrobio.javautils.UtilVCF2DIST
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
package ciat.agrobio.javautils;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import ciat.agrobio.core.CalculateDistancesCOSINE;
import ciat.agrobio.core.VariantProcessor;
import ciat.agrobio.io.VCFManager;

@Parameters(commandDescription = "VCF2DIST")
public class UtilVCF2DIST {

	private static UtilVCF2DIST instance = new UtilVCF2DIST();

	private UtilVCF2DIST() {
	}

	public static UtilVCF2DIST getInstance() {
		return instance;
	}

	public static String getUtilName() {
		return "VCF2DIST";
	}

	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(description = "Input_File", required = true)
	private String inputFileName;

	@Parameter(names = { "--numberOfThreads", "-t" })
	private int numOfThreads = 1;
	
	@Parameter(names = { "--ignoremissing", "-m" })
	private boolean ignoremissing = false;
	
	@Parameter(names={"--onlyHets", "-h"})
	private boolean onlyHets = false;
	
	@Parameter(names={"--ignoreHets", "-g"})
	private boolean ignoreHets = false;

	public void go() {
		try {
			int cpus = Runtime.getRuntime().availableProcessors();
			int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
			//System.err.println("cpus=" + cpus);
			//System.err.println("using=" + usingThreads);

			ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>> samplesToVariantsData = new ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>>();

			CountDownLatch startSignal = new CountDownLatch(1);
			CountDownLatch doneSignal = new CountDownLatch(usingThreads + 1);

			//System.err.println(GeneralTools.time() + " START ");

			ExecutorService pool = Executors.newFixedThreadPool(usingThreads + 1);

			Map<Integer, VariantProcessor> variantProcessors = new HashMap<Integer, VariantProcessor>();
			VCFManager vcfm = new VCFManager(inputFileName, startSignal, doneSignal, samplesToVariantsData);
			pool.execute(vcfm);

			VariantProcessor.resetCounters();
			// Starting threads
			for (int i = 0; i < usingThreads; i++) {
				VariantProcessor vp = new VariantProcessor(samplesToVariantsData, vcfm, startSignal, doneSignal);
				variantProcessors.put(vp.getId(), vp);
				pool.execute(vp);
			}

			doneSignal.await();
			pool.shutdown();

			List<String> sampleNames = new ArrayList<String>();
			for (int i = 9; i < vcfm.getHeaderData().size(); i++) {
				sampleNames.add(vcfm.getHeaderData().get(i));
			}

			//System.err.printf("\rProcessed variants : \t%8d\n", VariantProcessor.getVariantCount().get());

			// Calculate distances
			//double[][] distances = SampleVariantTools.calculateDistances(samplesToVariantsData, sampleNames, vcfm, ignoreHets, onlyHets, false);
			CalculateDistancesCOSINE fj = new CalculateDistancesCOSINE();
			double[][] distances = fj.calculateDistances(usingThreads, samplesToVariantsData, sampleNames, vcfm, ignoreHets, onlyHets, ignoremissing);
			
			// Print data
			int sampleCounter = 0;
			DecimalFormat df = new DecimalFormat("#.############"); 
			System.out.println(sampleNames.size()+"\t"+vcfm.getVariantIds().size());
			for(int i=0; i<sampleNames.size(); i++) {
				String sampleName1 = sampleNames.get(i);
				System.out.print(sampleName1);
				for (int j=0;j<sampleNames.size();j++) {
					String sampleName2 = sampleNames.get(j);
					System.out.print("\t"+df.format(distances[i][j]));
				}
				System.out.println("");
				sampleCounter++;
				//if(++sampleCounter % 10 == 0) System.err.println(GeneralTools.time()+" Samples Processed : \t"+sampleCounter);
			}
			//System.err.println(GeneralTools.time()+" Samples Processed : \t"+sampleCounter);
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
}
