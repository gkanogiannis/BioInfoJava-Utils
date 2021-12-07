/*
 *
 * BioInfoJava-Utils ciat.agrobio.javautils.UtilVCF2CSV
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

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import ciat.agrobio.core.GeneralTools;
import ciat.agrobio.core.GenotypeEncoder;

@Parameters(commandDescription = "VCF2CSV")
public class UtilVCF2CSV {
	
	private static UtilVCF2CSV instance = new UtilVCF2CSV();
	
	private UtilVCF2CSV() {}
	
	public static UtilVCF2CSV getInstance() {return instance;}
	
	public static String getUtilName() {
		return "VCF2CSV";
	}
	
	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(description = "Input_File", required = true)
	private List<String> inputFileNames = new ArrayList<String>();
	
	public void go() {
		try {
			InputStream fis = Files.newInputStream(Paths.get(inputFileNames.get(0)));
			BufferedReader br;
			if(inputFileNames.get(0).endsWith(".gz")) {
				br = new BufferedReader(new InputStreamReader(new GZIPInputStream(fis),"UTF-8"));
			}
			else {
				br = new BufferedReader(new InputStreamReader(fis, "UTF-8"));
			}
			
			int varCounter = 0;
			String line; 
			List<String> headerData = null;
			while((line=br.readLine()) != null) {
				if(line.startsWith("##")) continue;
				//Print header Individual names
				else if (line.startsWith("#")) {
					//System.err.println(line);
					headerData = Arrays.asList(line.split("\\t"));
					int numOfInd = headerData.size() - 9;
					System.err.println("Num of Ind = "+numOfInd);
					System.out.print("Variant");
					for(int i=9; i<headerData.size(); i++) {
						System.out.print(","+headerData.get(i));
					}
					System.out.println("");
					continue;
				}
				else {
					//Print Genotype Code for each Individual. 
					List<String> variantData = Arrays.asList(line.split("\\t"));
					String varName = variantData.get(0)+"_"+variantData.get(1);
					String format = variantData.get(8);
					int indexGT = Arrays.asList(format.split(":")).indexOf("GT");
					
					System.out.print(varName+"=Hom_Ref");
					for(int i=9; i<headerData.size(); i++) {
						String GT = variantData.get(i).split(":")[indexGT];
						String GTCodeString = "?";
					   	Integer[] GTCode = GenotypeEncoder.encodeGTComplex(GT);
					   	if(GTCode!=null) {
					   		GTCodeString = String.valueOf(GTCode[0]);
					   	}
						System.out.print("," + GTCodeString);
					}
					System.out.print(varName+"=Het");
					for(int i=9; i<headerData.size(); i++) {
						String GT = variantData.get(i).split(":")[indexGT];
						String GTCodeString = "?";
					   	Integer[] GTCode = GenotypeEncoder.encodeGTComplex(GT);
					   	if(GTCode!=null) {
					   		GTCodeString = String.valueOf(GTCode[1]);
					   	}
						System.out.print("," + GTCodeString);
					}
					System.out.print(varName+"=Hom_Alt");
					for(int i=9; i<headerData.size(); i++) {
						String GT = variantData.get(i).split(":")[indexGT];
						String GTCodeString = "?";
					   	Integer[] GTCode = GenotypeEncoder.encodeGTComplex(GT);
					   	if(GTCode!=null) {
					   		GTCodeString = String.valueOf(GTCode[2]);
					   	}
						System.out.print("," + GTCodeString);
					}
					System.out.println("");
					if(++varCounter % 1000 == 0) System.err.println(GeneralTools.time()+" Variants Processed : \t"+varCounter);
				}
				
			}
			br.close();
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
}
