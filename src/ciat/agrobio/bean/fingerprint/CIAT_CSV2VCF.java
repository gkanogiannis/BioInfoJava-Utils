/*
 *
 * BioInfoJava-Utils ciat.agrobio.bean.fingerprint.CIAT_CSV2VCF
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
import java.util.Arrays;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class CIAT_CSV2VCF {

	public static void main(String[] args) {
		try {
			//Load markers data from first input csv file
			//NCBIAssayID(ss#),markerID-shortID,chr,position,Ref
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
			
			//Load genotype data from third input csv file. First field is sample name and it is skipped.
			//sample X marker
			ArrayList<String[]> genotype_data = new ArrayList<String[]>();
			scanner = new Scanner(new FileInputStream(args[2]));
			while(scanner.hasNextLine()) {
				String[] line = scanner.nextLine().trim().split(",");
				genotype_data.add(Arrays.copyOfRange(line, 1, line.length));
			}
			scanner.close();
			
			//Write the vcf file. Converting A:A etc to 0/0,1/0,1/1 genotypes
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
					id = marker_data[1];
				if(id.isEmpty())
					id = ".";
				
				String refalt;
				String refbase;
				String altbase;
				try {
					Matcher m = Pattern.compile("_[ATCG]_[ATCG]_").matcher(marker_data[1]);
					if(m.find())
						refalt = m.group().replaceAll("_", "");
					else
						refalt = ".";
					refbase = refalt.substring(0, 1);
					altbase = refalt.substring(1, 2);
				} 
				catch (Exception e) {
					refbase = marker_data[4];
					altbase = ".";
				}
				
				System.out.print(marker_data[2]+"\t"+
								 marker_data[3].replaceAll(",","").replaceAll("\\.","")+"\t"+
								 id+"\t"+
								 refbase+"\t"+
								 altbase+"\t"+
								 "."+"\t"+
								 "."+"\t"+
								 "NO_INFO"+"\t"+
								 "GT");
				for(int j=0; j<samples_names.size(); j++) {
					String genotype_original = genotype_data.get(j)[i];
					String genotype = "./.";
					if(genotype_original.equalsIgnoreCase("?") || genotype_original.isEmpty() || !genotype_original.contains(":")) {
						genotype = "./.";
					}
					else {
						String al1 = genotype_original.split(":")[0].trim();
						String al2 = genotype_original.split(":")[1].trim();
						if(!al1.equalsIgnoreCase(al2)) {
							//heterozygote;
							if(al1.equalsIgnoreCase(refbase))
								genotype = "0/1";
							else if(al2.equalsIgnoreCase(refbase))
								genotype = "1/0";
							else genotype = "0/1";
						}
						else {
							//homozygote
							if(al1.equalsIgnoreCase(refbase) || al1.equalsIgnoreCase(RC(refbase)))
								genotype = "0/0";
							else if (al1.equalsIgnoreCase(altbase) || al1.equalsIgnoreCase(RC(altbase)))
								genotype = "1/1";
							else
								genotype = "0/0";
						}
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
