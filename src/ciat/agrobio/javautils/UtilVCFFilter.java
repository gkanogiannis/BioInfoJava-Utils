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
import ciat.agrobio.core.Variant;

@Parameters(commandDescription = "VCFFilter")
public class UtilVCFFilter {
	
	private static UtilVCFFilter instance = new UtilVCFFilter();
	
	private UtilVCFFilter() {}
	
	public static UtilVCFFilter getInstance() {return instance;}
	
	public static String getUtilName() {
		return "VCFFilter";
	}
	
	@Parameter(names = "--help", help = true)
	private boolean help;
	
	@Parameter(description = "Input_File", required = true)
	private List<String> inputFileNames = new ArrayList<String>();
	
	@Parameter(names={"--heterozygosityFilter","-h"})
	private double heterFilter = 0.10;
	
	@Parameter(names={"--missingdataFilter","-m"})
	private double misDataMaxFilter = 0.80;
	
	@SuppressWarnings("unused")
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
			
			//Comment
			String line; 
			List<String> headerData = null;
			int varCounter = 0;
			int varCounterRemain = 0;
			double avrgHet = 0.0;
			double avrgMid = 0.0;
			while((line=br.readLine()) != null) {
				if(line.startsWith("##")) {
					System.out.println(line);
					continue;
				}
				//Header
				else if (line.startsWith("#")) {
					headerData = Arrays.asList(line.split("\\s+"));
					System.out.println(line);
					continue;
				}
				//Variants
				else {
					Variant var = new Variant(++varCounter,line.split("\\s+"));
					//Heter filter
					double heterozygosity = var.calculateHeterozygosity();
					if(heterozygosity<heterFilter) {
						//Missing data filter
						double misData = var.calculateMissingData();
						if(misData<misDataMaxFilter) {
							++varCounterRemain;
							avrgHet += heterozygosity;
							avrgMid += misData;
							System.out.println(var.toString());
						}
						else {
							//System.err.println(var.getName()+"\tMiD="+var.calculateMissingData());
						}
					}
					else {
						//System.err.println(var.getName()+"\tHet="+var.calculateHeterozygosity());
					}
					if(varCounter % 10000 == 0) System.err.println(GeneralTools.time()+" Variants Processed : \t"+varCounter);
				}
			}
			System.err.println("Variants Filtered      = "+(varCounter-varCounterRemain));
			System.err.println("Variants Remain        = "+(varCounterRemain));
			System.err.println("Average Heterozygosity = "+(avrgHet/varCounterRemain));
			System.err.println("Average Missing Data   = "+(avrgMid/varCounterRemain));
			
			br.close();
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
}

