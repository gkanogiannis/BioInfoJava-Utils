/*
 *
 * BioInfoJava-Utils ciat.agrobio.javautils.UtilVCFsIntersection
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
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import ciat.agrobio.core.GeneralTools;
import ciat.agrobio.core.Variant;

@Parameters(commandDescription = "VCFsIntersection")
public class UtilVCFsIntersection {
	
	private static UtilVCFsIntersection instance = new UtilVCFsIntersection();
	
	private UtilVCFsIntersection() {}
	
	public static UtilVCFsIntersection getInstance() {return instance;}
	
	public static String getUtilName() {
		return "VCFsIntersection";
	}
	
	@Parameter(names = "--help", help = true)
	private boolean help;
	
	@Parameter(names={"--inputVCF1", "-i1"}, required=true)
	private String inputFileName1;
	
	@Parameter(names={"--inputVCF2", "-i2"}, required=true)
	private String inputFileName2;
	
	@Parameter(names={"--intersect", "-x"})
	private boolean interesct = false; // then join
	
	@Parameter(names={"--noINFOfield", "-n"})
	private boolean noInfoField = false; // dont print the INFO field (just "INFO")
	
	public void go() {
		try {
			InputStream fis1 = Files.newInputStream(Paths.get(inputFileName1));
			BufferedReader br1;
			if(inputFileName1.endsWith(".gz")) {
				br1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(fis1),"UTF-8"));
			}
			else {
				br1 = new BufferedReader(new InputStreamReader(fis1, "UTF-8"));
			}
			InputStream fis2 = Files.newInputStream(Paths.get(inputFileName2));
			BufferedReader br2;
			if(inputFileName2.endsWith(".gz")) {
				br2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(fis2),"UTF-8"));
			}
			else {
				br2 = new BufferedReader(new InputStreamReader(fis2, "UTF-8"));
			}
			
			//Comment1
			String line1; 
			List<String> headerData1 = null;
			while((line1=br1.readLine()) != null) {
				if(line1.startsWith("##")) {
					System.out.println(line1);
					continue;
				}
				//Header1
				else if (line1.startsWith("#")) {
					headerData1 = Arrays.asList(line1.split("\\s+"));
					break;
				}
			}
			
			//Comment2
			String line2; 
			List<String> headerData2 = null;
			while((line2=br2.readLine()) != null) {
				if(line2.startsWith("##")) {
					System.out.println(line2);
					continue;
				}
				//Header1
				else if (line2.startsWith("#")) {
					headerData2 = Arrays.asList(line2.split("\\s+"));
					break;
				}
			}
			
			//New Header
			System.out.print(headerData1.get(0));
			//Header for vcf1
			for(int i=1; i<headerData1.size(); i++) {
				System.out.print("\t" + headerData1.get(i));
			}
			//Header for vcf2
			for(int i=9; i<headerData2.size(); i++) {
				System.out.print("\t" + headerData2.get(i));
			}
			System.out.println("");
			
			//New Variants
			int varCounter = 0;
			if(interesct) { //Intersect
				line1=br1.readLine(); varCounter++;
				line2=br2.readLine(); varCounter++;
				while(line1!=null && line2!=null) {
					Variant var1 = new Variant(0, line1.split("\\t"));
					Variant var2 = new Variant(0, line2.split("\\t"));
					int compare = var1.compareTo(var2);
					if(compare==0) {
						//Print first 7 fields
						for(int i=0; i<=6; i++) {
							System.out.print(var1.getData()[i]+"\t");
						}
						//Print INFO field
						if(noInfoField) {
							System.out.print("INFO=NO");
						}
						else {
							System.out.print("VCF=1;"+var1.getData()[7]);
							System.out.print(";VCF=2;"+var2.getData()[7]);
						}
						//Print FORMAT field
						System.out.print("\t"+var1.getData()[8]);
						//Print genotypes for vcf1
						for(int i=9; i<var1.getData().length; i++) {
							System.out.print("\t"+var1.getData()[i]);
						}
						//Print genotypes for vcf2
						for(int i=9; i<var2.getData().length; i++) {
							System.out.print("\t"+var2.getData()[i]);
						}
						
						System.out.println("");
						line1=br1.readLine(); varCounter++;
						line2=br2.readLine(); varCounter++;
					}
					else if(compare<0) {
						line1=br1.readLine(); varCounter++;
					}
					else {
						line2=br2.readLine(); varCounter++;
					}
						
					if(varCounter % 10000 == 0) System.err.println(GeneralTools.time()+" Variants Processed : \t"+varCounter);
				}
			}
			else { //Join
				
			}
			
			br1.close();
			br2.close();
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
}
