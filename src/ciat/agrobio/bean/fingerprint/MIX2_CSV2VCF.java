/*
 *
 * BioInfoJava-Utils ciat.agrobio.bean.fingerprint.MIX2_CSV2VCF
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
package ciat.agrobio.bean.fingerprint;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.Scanner;

public class MIX2_CSV2VCF {

	public static void main(String[] args) {
		try {
			//Load markers data from first input csv file
			//NCBI_ss,Chr,Pos,Ref,Alt
			ArrayList<String[]> markers_data = new ArrayList<String[]>();
			Scanner scanner = new Scanner(new FileInputStream(args[0]));
			//Skip header line
			scanner.nextLine();
			while(scanner.hasNextLine()) {
				markers_data.add(scanner.nextLine().trim().split(","));
			}
			scanner.close();
			
			//Load samples names from second input txt file
			ArrayList<String> samples_names = new ArrayList<String>();
			scanner = new Scanner(new FileInputStream(args[1]));
			while(scanner.hasNextLine()) {
				String[] line = scanner.nextLine().trim().split("\t");
				samples_names.add(line[0].trim());
			}
			scanner.close();
			
			//Load genotype data from third input csv file. 
			//marker X sample
			ArrayList<String[]> genotype_data = new ArrayList<String[]>();
			scanner = new Scanner(new FileInputStream(args[2]));
			while(scanner.hasNextLine()) {
				String[] line = scanner.nextLine().trim().split(",");
				genotype_data.add(line);
			}
			scanner.close();
			
			//Write the vcf file. Converting 0,1,2 to 0/0,1/0,1/1 genotypes
			System.out.println("##fileformat=VCFv4.2");
			System.out.println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
			System.out.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
			for(String sampleName : samples_names){
				System.out.print("\t"+sampleName);
			}
			System.out.println("");
			
			for(int i=0; i<markers_data.size(); i++) {
				String[] marker_data = markers_data.get(i);
				
				String id = marker_data[0];
				if(id.isEmpty())
					id = ".";
				
				String refbase = marker_data[3];
				String altbase = marker_data[4];
				
				System.out.print(marker_data[1]+"\t"+
								 marker_data[2]+"\t"+
								 id+"\t"+
								 refbase+"\t"+
								 altbase+"\t"+
								 "."+"\t"+
								 "."+"\t"+
								 "NO_INFO"+"\t"+
								 "GT");
				for(int j=0; j<samples_names.size(); j++) {
					String genotype_original = genotype_data.get(i)[j];
					String genotype = "./.";
					if(genotype_original.equalsIgnoreCase("?") || genotype_original.isEmpty() || genotype_original.equalsIgnoreCase("NA")) {
						genotype = "./.";
					}
					else {
						if(genotype_original.equalsIgnoreCase("0"))
							genotype = "0/0";
						else if(genotype_original.equalsIgnoreCase("1"))
							genotype = "0/1";
						else if(genotype_original.equalsIgnoreCase("2"))
							genotype = "1/1";
					}
					System.out.print("\t"+genotype);
				}
				System.out.println("");
			}	
		} 
		catch (Exception e) {
			e.printStackTrace();
		}

	}
	
	@SuppressWarnings("unused")
	private static String RC(String in) {
		if(in.equalsIgnoreCase("A"))
			return "T";
		else if(in.equalsIgnoreCase("T"))
			return "A";
		else if(in.equalsIgnoreCase("C"))
			return "G";
		else if(in.equalsIgnoreCase("G"))
			return "C";
		else return "";

	}

}
