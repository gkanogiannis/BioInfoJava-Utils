package ciat.agrobio.syngenta2017;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

public class HxG2ARFF {

	private HxG2ARFF(String inputFileName, String separator) {
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
			

			//Print attributes (variants)
			System.out.println("@relation unnamed\n");
			System.out.println("@attribute\t" + "Hybrid" + "\tstring");
			Iterator<String> iteratorVariantNames = variantNames.iterator();
		    while(iteratorVariantNames.hasNext()) {
		    	String variantName = iteratorVariantNames.next();
		    	System.out.println("@attribute\t" + variantName + "\tnumeric");
		    }
		    //Print data
			System.out.println("\n@data");
			for(String individualName : individualNames) {
				System.out.print(individualName);
				iteratorVariantNames = variantNames.iterator();
			    while(iteratorVariantNames.hasNext()) {
			    	String variantName = iteratorVariantNames.next();
			    	if(((String)variantsToIndividualsData.get(variantName).get(individualName)).equals("NA")) {
			    		System.out.print(","+"?");
			    	}
			    	else {
			    		System.out.print(","+(String)variantsToIndividualsData.get(variantName).get(individualName));
			    	}
			    }
				System.out.println("");
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		HxG2ARFF tool = new HxG2ARFF(args[0], args[1]);
		tool.go();
	}
}
