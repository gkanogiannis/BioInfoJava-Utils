/*
 *
 * BioInfoJava-Utils ciat.agrobio.javautils.UtilVCFKeepVariants
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

import java.io.File;
import java.io.FileWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import ciat.agrobio.core.GeneralTools;
import ciat.agrobio.core.Variant;
import ciat.agrobio.core.VariantProcessor;
import ciat.agrobio.io.VCFManager;

@Parameters(commandDescription = "VCFKeepVariants")
public class UtilVCFKeepVariants {

	private static UtilVCFKeepVariants instance = new UtilVCFKeepVariants();

	private UtilVCFKeepVariants() {
	}

	public static UtilVCFKeepVariants getInstance() {
		return instance;
	}

	public static String getUtilName() {
		return "VCFKeepVariants";
	}

	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(description = "Input_File", required = true)
	private String inputFileName;
	
	@Parameter(names={"--variantList", "-v"}, description = "Variant_List_File", required = true)
	private String variantListFileName;
	
	@Parameter(names = { "--numberOfThreads", "-t" })
	private int numOfThreads = 1;

	@Parameter(names={"--ignoreHets", "-g"})
	private boolean ignoreHets = false;
	
	@Parameter(names={"--output", "-o"}, required=true)
	private String output = null;
	
	private static VCFManager vcfm;

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
		
			//Read variantList to keep
			List<String> variantVCFIdKeep = new ArrayList<String>();
			Scanner scanner = new Scanner(new File(variantListFileName));
			while(scanner.hasNextLine()) {
				variantVCFIdKeep.add(scanner.nextLine().trim().split("[\\s,=]+")[0].trim());
			}
			scanner.close();
			
			//Write final VCF keep only the variants in variantVCFIdKeep
			Writer writer = new FileWriter(output);
			//System.out.println("##fileformat=VCFv4.2");
			//System.out.println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
			writer.write(String.join("\n", vcfm.getCommentData())+"\n");
			writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
			for(String sampleName : sampleNames){
				writer.write("\t"+sampleName);
			}
			writer.write("\n");
			Iterator<Integer> iterator = vcfm.getVariantIds().iterator();
			while (iterator.hasNext()) {
				Variant variant = vcfm.getStaticVariants().get(iterator.next());
				if(variant==null)
					continue;
				else if(variant.getVariantVCFId().trim().equals("\\.") && !variantVCFIdKeep.contains(String.valueOf(variant.getVariantId())) )
					continue;
				else if(!variantVCFIdKeep.contains(variant.getVariantVCFId().trim()))
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
					Integer[] GTCode = (Integer[])samplesToVariantsData.get(sampleName).get(variant.getVariantId());
					if(GTCode==null || GTCode.length!=3)
						writer.write("\t"+"./.");
					else if(ignoreHets && (GTCode[1]==1))
						writer.write("\t"+"./.");
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
	
}
