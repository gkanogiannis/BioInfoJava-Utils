/*
 *
 * BioInfoJava-Utils ciat.agrobio.core.SampleVariantTools
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
package ciat.agrobio.core;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;

import ciat.agrobio.io.VCFManager;

public class SampleVariantTools {

	public static double[][] calculateDistances(ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>> samplesToVariantsData, List<String> sampleNames, VCFManager vcfm, boolean ignoreHets, boolean onlyHets, boolean euclidean) {
		try {
			DecimalFormat df = new DecimalFormat("#.############"); 
			int numOfSamples = sampleNames.size();
			double[][] distances = new double[numOfSamples][numOfSamples]; 
			int sampleCounter = 0;
			//System.err.println(numOfSamples);
			for (int i=0;i<numOfSamples;i++) {
				String sampleName1 = sampleNames.get(i);
				//System.err.print(sampleName1);
				for (int j=0;j<numOfSamples;j++) {
					String sampleName2 = sampleNames.get(j);
					double dot = 0.0;
					double eucl = 0.0;
					double norm1 = 0.0;
					double norm2 = 0.0;
					Iterator<Integer> iterator = vcfm.getVariantIds().iterator();
					while (iterator.hasNext()) {
						Integer variantId = iterator.next();
						Integer[] GTCode1 = (Integer[])samplesToVariantsData.get(sampleName1).get(variantId);
						Integer[] GTCode2 = (Integer[])samplesToVariantsData.get(sampleName2).get(variantId);
						
						if((GTCode1==null && GTCode2==null) || (GTCode1.length==0 && GTCode2.length==0))
							continue;
						else if(ignoreHets && ((GTCode1!=null&&GTCode1.length==3&&GTCode1[1]==1) || (GTCode2!=null&&GTCode2.length==3&&GTCode2[1]==1)))
							continue;
						else if(onlyHets && ((GTCode1!=null&&GTCode1.length==3&&GTCode1[1]!=1) && (GTCode2!=null&&GTCode2.length==3&&GTCode2[1]!=1)))
							continue;
						else if((GTCode1==null || GTCode1.length==0) && (GTCode2!=null&&GTCode2.length==3))
							norm2 += ( (GTCode2[0]*GTCode2[0]) + (GTCode2[1]*GTCode2[1]) + (GTCode2[2]*GTCode2[2]) );
						else if((GTCode2==null || GTCode2.length==0) && (GTCode1!=null&&GTCode1.length==3))
							norm1 += ( (GTCode1[0]*GTCode1[0]) + (GTCode1[1]*GTCode1[1]) + (GTCode1[2]*GTCode1[2]) );
						else {
							dot += ( (GTCode1[0]*GTCode2[0]) + (GTCode1[1]*GTCode2[1]) + (GTCode1[2]*GTCode2[2]) );
							eucl += ( Math.pow(GTCode1[0]-GTCode2[0], 2.0) + Math.pow(GTCode1[1]-GTCode2[1], 2.0) + Math.pow(GTCode1[2]-GTCode2[2], 2.0) );
							norm1 += ( (GTCode1[0]*GTCode1[0]) + (GTCode1[1]*GTCode1[1]) + (GTCode1[2]*GTCode1[2]) );
							norm2 += ( (GTCode2[0]*GTCode2[0]) + (GTCode2[1]*GTCode2[1]) + (GTCode2[2]*GTCode2[2]) );
						}
					}
					if(euclidean) {
						eucl = Math.sqrt(eucl);
						distances[i][j] = eucl;
						//System.out.print("\t"+df.format(eucl));
					}
					else {
						double cosine = dot / (Math.sqrt(norm1)*Math.sqrt(norm2));
						if(Double.isInfinite(cosine) || Double.isNaN(cosine)) cosine = 0.0;
						
						double distance = 1.0 - cosine;
						if(distance<=0.0) distance = 0.0;
						
						//double distance = -1.0*Math.log(cosine);
						//if(distance<=0.0) distance = 0.0;
						
						distances[i][j] = Double.parseDouble(df.format(distance));
						
						//System.out.print("\t"+df.format(distance));
					}
				}
				if(++sampleCounter % 50 == 0) System.err.println(GeneralTools.time()+" Samples Processed : \t"+sampleCounter);
			}
			System.err.println(GeneralTools.time()+" Samples Processed : \t"+sampleCounter);
			
			return distances;
		} 
		catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	public static HashMap<String, Double> calculateMissingData(ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>> samplesToVariantsData, List<String> sampleNames, VCFManager vcfm) {
		try {
			HashMap<String, Double> missingData = new HashMap<String, Double>();
			DecimalFormat df = new DecimalFormat("#.############");
			for (int i=0;i<sampleNames.size();i++) {
				String sampleName = sampleNames.get(i);
				double currentMissingData = 0.0;
				Iterator<Integer> iterator = vcfm.getVariantIds().iterator();
				while (iterator.hasNext()) {
					Integer variantId = iterator.next();
					Integer[] GTCode = (Integer[])samplesToVariantsData.get(sampleName).get(variantId);
					if(GTCode==null || GTCode.length==0)
						currentMissingData++;
				}
				missingData.put(sampleName, currentMissingData/(double)vcfm.getVariantIds().size());
				System.err.println(sampleName+"\t"+df.format(missingData.get(sampleName)));
			}
			
			return missingData;
		} 
		catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
}
