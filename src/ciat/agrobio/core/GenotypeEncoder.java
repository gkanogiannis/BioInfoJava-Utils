package ciat.agrobio.core;

import java.util.HashMap;
import java.util.Map;

public class GenotypeEncoder {
	
	private static Map<String, Integer[]> mapComplex;
	
	
	static {
		mapComplex = new HashMap<String, Integer[]>();
		mapComplex.put("0/0", new Integer[]{1, 0, 0});
		mapComplex.put("1/1", new Integer[]{0, 0, 1});
		mapComplex.put("0/1", new Integer[]{0, 1, 0});
		mapComplex.put("1/0", new Integer[]{0, 1, 0});
		mapComplex.put("./.", new Integer[]{});
		
		mapComplex.put("0|0", new Integer[]{1, 0, 0});
		mapComplex.put("1|1", new Integer[]{0, 0, 1});
		mapComplex.put("0|1", new Integer[]{0, 1, 0});
		mapComplex.put("1|0", new Integer[]{0, 1, 0});
		mapComplex.put(".|.", new Integer[]{});
	}
	
	/*
	static {
		mapComplex = new HashMap<String, Integer[]>();
		mapComplex.put("0/0", new Integer[]{1, 0, 1});
		mapComplex.put("1/1", new Integer[]{-1, 0, -1});
		mapComplex.put("0/1", new Integer[]{0, 1, 0});
		mapComplex.put("1/0", new Integer[]{0, 1, 0});
		mapComplex.put("./.", new Integer[]{});
	}
	*/
	
	public static String genotypeToAlleles(String GT, String Ref, String Alt, boolean mode2) {
		if(mode2) {
			if(GT.equals("0/0"))
				return Ref+"/"+Ref;
			else if(GT.equals("1/1"))
				return Alt+"/"+Alt;
			else if(GT.equals("0/1") || GT.equals("1/0"))
				return Ref+"/"+Alt;
			else
				return "N/N";
		}
		else {
			if(GT.equals("0/0"))
				return Ref;
			else if(GT.equals("1/1"))
				return Alt;
			else if(GT.equals("0/1") || GT.equals("1/0"))
				return Ref+"/"+Alt;
			else
				return "N";
		}
	}
	
	public static Integer[] encodeGTComplex(String GT) {
		Integer[] ret =  mapComplex.get(GT);
		if(ret!=null)
			return ret;
		else
			return new Integer[]{};
	}
	
	
	public static String encodeGTComplexReverse(Integer[] GT) {
		try {
			if(GT==null || GT.length!=3)
				return "./.";
			else {
				if(GT[0]==1 && GT[1]==0 && GT[2]==0) return "0/0";
				else if(GT[0]==0 && GT[1]==0 && GT[2]==1) return "1/1";
				else if(GT[0]==0 && GT[1]==1 && GT[2]==0) return "0/1";
				else return "./.";
			}
		}
		catch (Exception e) {
			return "./.";
		}
	}
	
	/*
	public static String encodeGTComplexReverse(Integer[] GT) {
		try {
			if(GT==null || GT.length!=3)
				return "./.";
			else {
				if(GT[0]==1 && GT[1]==0 && GT[2]==1) return "0/0";
				else if(GT[0]==-1 && GT[1]==0 && GT[2]==-1) return "1/1";
				else if(GT[0]==0 && GT[1]==1 && GT[2]==0) return "0/1";
				else return "./.";
			}
		}
		catch (Exception e) {
			return "./.";
		}
	}
	*/
	
	/*
	public static Double encodeGTSimple(String GT) {
		try {
			Integer i1 = Integer.parseInt(GT.split("/")[0]);
			Integer i2 = Integer.parseInt(GT.split("/")[1]);
			return (i1+i2)/1.0;
		} 
		catch (Exception e) {
			return -1.0;
		}
	}
	
	public static String encodeGTSimpleReverse(Double GT) {
		try {
			if(GT==0.0)
				return "0/0";
			else if(GT==2.0)
				return "1/1";
			else if(GT==1.0)
				return "0/1";
			else
				return "./.";
		} 
		catch (Exception e) {
			return "./.";
		}
	}
	*/
}
