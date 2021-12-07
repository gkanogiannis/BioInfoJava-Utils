/*
 *
 * BioInfoJava-Utils ciat.agrobio.javautils.UtilVCF2TREE
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
import ciat.agrobio.hcluster.HierarchicalCluster;
import ciat.agrobio.io.VCFManager;

@Parameters(commandDescription = "VCF2TREE")
public class UtilVCF2TREE {

	private static UtilVCF2TREE instance = new UtilVCF2TREE();

	private UtilVCF2TREE() {
	}

	public static UtilVCF2TREE getInstance() {
		return instance;
	}

	public static String getUtilName() {
		return "VCF2TREE";
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
			
			//HCluster
			//String treeString = hclusteringTree(individualNames.toArray(new String[individualNames.size()]), distances);
			String treeString = HierarchicalCluster.hclusteringTree(sampleNames.toArray(new String[sampleNames.size()]), distances);
			System.out.println(treeString);
			
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
