/*
 *
 * BioInfoJava-Utils ciat.agrobio.javautils.UtilVCFRemoveClones
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

import java.io.FileWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import ciat.agrobio.core.SampleVariantTools;
import ciat.agrobio.core.CalculateDistancesCOSINE;
import ciat.agrobio.core.GeneralTools;
import ciat.agrobio.core.Variant;
import ciat.agrobio.core.VariantProcessor;
import ciat.agrobio.hcluster.HierarchicalCluster;
import ciat.agrobio.io.VCFManager;

@Parameters(commandDescription = "VCFRemoveClones")
public class UtilVCFRemoveClones {

	private static UtilVCFRemoveClones instance = new UtilVCFRemoveClones();

	private UtilVCFRemoveClones() {
	}

	public static UtilVCFRemoveClones getInstance() {
		return instance;
	}

	public static String getUtilName() {
		return "VCFRemoveClones";
	}

	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(description = "Input_File", required = true)
	private String inputFileName;

	@Parameter(names = { "--cutHeight", "-c" })
	private Double cutHeight = 0.075;
	
	@Parameter(names = { "--numberOfThreads", "-t" })
	private int numOfThreads = 1;
	
	@Parameter(names={"--onlyHets", "-h"})
	private boolean onlyHets = false;
	
	@Parameter(names={"--ignoreHets", "-g"})
	private boolean ignoreHets = false;
	
	@Parameter(names={"--output", "-o"}, required=true)
	private String output = null;
	
	private static VCFManager vcfm;
	
	@Parameter(names={"--extra", "-e"})
	private boolean extra = false;

	public void go() {
		try {
			int cpus = Runtime.getRuntime().availableProcessors();
			int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
			System.err.println("cpus=" + cpus);
			System.err.println("using=" + usingThreads);

			ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>> samplesToVariantsData = new ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>>();

			CountDownLatch startSignal = new CountDownLatch(1);
			CountDownLatch doneSignal = new CountDownLatch(usingThreads + 1);

			System.err.println(GeneralTools.time() + " START ");

			ExecutorService pool = Executors.newFixedThreadPool(usingThreads + 1);

			Map<Integer, VariantProcessor> variantProcessors = new HashMap<Integer, VariantProcessor>();
			vcfm = new VCFManager(inputFileName, startSignal, doneSignal, samplesToVariantsData);
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
			
			System.err.printf("\rProcessed variants : \t%8d\n", VariantProcessor.getVariantCount().get());
		
			// Calculate distances
			//double[][] distances = SampleVariantTools.calculateDistances(samplesToVariantsData, sampleNames, vcfm, ignoreHets, onlyHets, false);
			CalculateDistancesCOSINE fj = new CalculateDistancesCOSINE();
			double[][] distances = fj.calculateDistances(usingThreads, samplesToVariantsData, sampleNames, vcfm, ignoreHets, onlyHets, false);
			
			// Calculate missing data
			//HashMap<String, Double> missingData = calculateMissingData(samplesToVariantsData, sampleNames);
			HashMap<String, Double> missingData = SampleVariantTools.calculateMissingData(samplesToVariantsData, sampleNames, vcfm);
			
			//HCluster
			//TreeMap<Integer, TreeSet<String>> clusters = hclusteringClusters(sampleNames.toArray(new String[sampleNames.size()]), distances, null, cutHeight);
			TreeMap<Integer, TreeSet<String>> clusters = HierarchicalCluster.hclusteringClusters(sampleNames.toArray(new String[sampleNames.size()]), distances, null, cutHeight, extra);
			
			List<String> sampleNamesKeep = new ArrayList<String>();
			//Remove from each cluster(except 0) all the samples except the one with the least missing data
			for(Entry<Integer, TreeSet<String>> entry : clusters.entrySet()){
				int clusterId = entry.getKey();
				if(clusterId==0) {
					TreeSet<String> cluster = entry.getValue();
					for(String name : cluster){
						sampleNamesKeep.add(name);
					}
					continue;
				}
				PriorityQueue<ClusterItem> clusterItems = new PriorityQueue<ClusterItem>();
				TreeSet<String> cluster = entry.getValue();
				System.err.println("Cluster\t"+clusterId+"="+cluster.size());
				for(String name : cluster){
					System.err.println("\t"+name+"\t"+missingData.get(name));
					clusterItems.add(new ClusterItem(name, missingData.get(name)));
				}
				//pop the head (least missing data)
				ClusterItem least = clusterItems.poll();
				System.err.println("Keep\t"+least.getName()+"\n");
				sampleNamesKeep.add(least.getName());
			}
			
			//Print kept samples
			System.out.println("Kept samples : "+sampleNamesKeep.size());
			for(String sampleName : sampleNamesKeep){
				System.out.println(sampleName);
			}
			
			//Write final VCF
			Writer writer = new FileWriter(output);
			//System.out.println("##fileformat=VCFv4.2");
			//System.out.println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
			writer.write(String.join("\n", vcfm.getCommentData())+"\n");
			writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
			for(String sampleName : sampleNames){
				if(sampleNamesKeep.contains(sampleName))
					writer.write("\t"+sampleName);
			}
			writer.write("\n");
			Iterator<Integer> iterator = vcfm.getVariantIds().iterator();
			while (iterator.hasNext()) {
				Variant variant = vcfm.getStaticVariants().get(iterator.next());
				if(variant==null)
					continue;
				writer.write(variant.getData()[0]+"\t"+
								 variant.getData()[1]+"\t"+
								 variant.getData()[2]+"\t"+
								 variant.getData()[3]+"\t"+
								 variant.getData()[4]+"\t"+
								 variant.getData()[5]+"\t"+
								 variant.getData()[6]+"\t"+
								 variant.getData()[7]+"\t"+
								 variant.getData()[8]);
				for(int i=0; i<sampleNames.size(); i++) {
					String sampleName = sampleNames.get(i);
					if(!sampleNamesKeep.contains(sampleName))
						continue;
					Integer[] GTCode = (Integer[])samplesToVariantsData.get(sampleName).get(variant.getVariantId());
					if(GTCode==null || GTCode.length!=3)
						writer.write("\t"+"./.");
					//else if(ignoreHets && (GTCode[1]==1))
						//writer.write("\t"+"./.");
					else
						writer.write("\t"+variant.getData()[i+9]);
				}
				writer.write("\n");
			}
			writer.close();
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private final class ClusterItem implements Comparable<ClusterItem>{
		private String name;
		private double missingData;
		
		public ClusterItem(String name, double missingData) {
			super();
			this.name = name;
			this.missingData = missingData;
		}

		public String getName() {
			return name;
		}

		@SuppressWarnings("unused")
		public void setName(String name) {
			this.name = name;
		}

		public double getMissingData() {
			return missingData;
		}

		@SuppressWarnings("unused")
		public void setMissingData(double missingData) {
			this.missingData = missingData;
		}

		@Override
		public int compareTo(ClusterItem o) {
			if(missingData<o.getMissingData())
				return -1;
			else if(missingData>o.getMissingData())
				return 1;
			else
				return 0;
		}
		
		
	}
}
