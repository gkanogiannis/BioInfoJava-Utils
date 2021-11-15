package ciat.agrobio.syngenta2017;

import java.io.BufferedReader;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import ciat.agrobio.core.GeneralTools;

public class HxG2DIST {

	private HxG2DIST(String inputFileName, String separator) {
		this.inputFileName = inputFileName;
		this.separator = separator;
	}

	private String inputFileName;
	private String separator;
	
	public void go() {
		try {
			BufferedReader br = new BufferedReader(new FileReader(inputFileName));
			String line = null;
			
			HashMap<String, HashMap<String, Object>> variantsToIndividualsData = new HashMap<String, HashMap<String, Object>>();
			
			List<String> variantNames = new ArrayList<String>();
			line = br.readLine();
			for (String i : line.split(separator)) {
				if(!i.trim().startsWith("G")) continue;
				variantNames.add(i.trim());
				variantsToIndividualsData.put(i.trim(), new HashMap<String, Object>());
			}
			
			List<String> individualNames = new ArrayList<String>();
			while ((line=br.readLine())!=null) {
				String[] data = line.split(separator);
				String individualName = data[0].trim();
				individualNames.add(individualName);
				for (int i=0; i<variantNames.size(); i++) {
					variantsToIndividualsData.get(variantNames.get(i)).put(individualName, data[i+1]);
				}
			}
			br.close();
			

			// Print data
			int sampleCounter = 0;
			DecimalFormat df = new DecimalFormat("#.############"); 
			System.out.println(individualNames.size());
			for (String individualName1 : individualNames) {
				System.out.print(individualName1);
				for (String individualName2 : individualNames) {
					double similarity = 0.0;
					double norm1 = 0.0;
					double norm2 = 0.0;
					Iterator<String> iterator = variantNames.iterator();
					while (iterator.hasNext()) {
						try {
							String variantName = iterator.next();
							Double GTCode1 = Double.parseDouble((String)variantsToIndividualsData.get(variantName).get(individualName1));
							Double GTCode2 = Double.parseDouble((String)variantsToIndividualsData.get(variantName).get(individualName2));
							if(GTCode1 == 0)
								continue;
							if(GTCode2 == 0)
								continue;
							similarity = similarity + (GTCode1*GTCode2);
							norm1 = norm1 + (GTCode1*GTCode1);
							norm2 = norm2 + (GTCode2*GTCode2);
						} 
						catch (Exception e) {
							continue;
						}
					}
					similarity = similarity / (Math.sqrt(norm1)*Math.sqrt(norm2));
					if(Double.isInfinite(similarity) || Double.isNaN(similarity)) similarity = 0.0;
					double distance = 1.0 - similarity;
					if(distance<0.0) distance = 0.0;
					
					System.out.print("\t"+df.format(distance));
				}
				System.out.println("");
				if(++sampleCounter % 10 == 0) System.err.println(GeneralTools.time()+" Samples Processed : \t"+sampleCounter);
			}
			System.err.println(GeneralTools.time()+" Samples Processed : \t"+sampleCounter);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		HxG2DIST tool = new HxG2DIST(args[0], args[1]);
		tool.go();
	}
}
