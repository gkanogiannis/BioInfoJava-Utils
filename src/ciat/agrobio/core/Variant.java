/*
 *
 * BioInfoJava-Utils ciat.agrobio.core.Variant
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

import java.util.Comparator;
import java.util.stream.IntStream;

public class Variant implements Comparable<Variant>{
	
	private int variantId;
	private String[] data;
	private static AlphanumericComparator alnumComparator;
	
	static {
		alnumComparator = new AlphanumericComparator();
	}
	
	public Variant(int variantId, String[] data) {
		this.variantId = variantId;
		this.data = data;
	}

	public int getVariantId() {
		return variantId;
	}
	
	public String getVariantVCFId() {
		return data[2];
	}
	
	public String[] getData() {
		return data;
	}
	
	public String getVariantName() {
		return data[0]+"\t|\t"+data[1]+"\t|\t"+data[3]+"\t|\t"+data[4];
	}
	
	public String getInfo() {
		return data[7];
	}
	
	public String getFormat() {
		return data[8];
	}
	
	public double calculateHeterozygosity() {
		double het = 0.0;
		for(int i=9; i<data.length; i++) {
			String GT = data[i].split(":")[0];
			if(GT.equalsIgnoreCase("1/0") || GT.equalsIgnoreCase("0/1")) {
				het += 1.0;
			}
		}
		return het/(data.length-9);
	}
	
	public double calculateMissingData() {
		double misdata = 0.0;
		for(int i=9; i<data.length; i++) {
			String GT = data[i].split(":")[0];
			if(GT.equalsIgnoreCase("./.")) {
				misdata += 1.0;
			}
		}
		return misdata/(data.length-9);
	}
	
	public double calculateCommonGenotypes(int includeIndex, int...excludeIndexes) {
		double common = 0.0;
		for(int i=9; i<data.length; i++) {
			final int ii = i;
			if(IntStream.of(excludeIndexes).anyMatch(x -> x == ii)) continue;
			if(includeIndex == i) continue;
			String GT1 = data[i].split(":")[0];
			String GT2 = data[includeIndex].split(":")[0];
			if(GT1.split("/")[0].equalsIgnoreCase(GT2.split("/")[0]) || GT1.split("/")[0].equalsIgnoreCase(GT2.split("/")[1])) {
				common += 0.5;
			}
			if(GT1.split("/")[1].equalsIgnoreCase(GT2.split("/")[0]) || GT1.split("/")[1].equalsIgnoreCase(GT2.split("/")[1])) {
				common += 0.5;
			}
		}
		return common/(data.length-9-1-excludeIndexes.length);
	}

	@Override
	public int compareTo(Variant other) {
		String chr1 = this.getVariantName().split("\t\\|\t")[0];
    	String chr2 = other.getVariantName().split("\t\\|\t")[0];
    	Integer pos1 = Integer.valueOf(this.getVariantName().split("\t\\|\t")[1]);
    	Integer pos2 = Integer.valueOf(other.getVariantName().split("\t\\|\t")[1]);
    	String ref1 = this.getVariantName().split("\t\\|\t")[2];
    	String ref2 = other.getVariantName().split("\t\\|\t")[2];
    	String alt1 = this.getVariantName().split("\t\\|\t")[3];
    	String alt2 = other.getVariantName().split("\t\\|\t")[3];
    	
    	int cComp = alnumComparator.compare(chr1, chr2);
    	//int cComp = chr1.compareTo(chr2);
        int pComp = pos1.compareTo(pos2);
        int rComp = ref1.compareTo(ref2);
        int aComp = alt1.compareTo(alt2);
        
        if (cComp != 0) {
           return cComp;
        } 
        else if(pComp != 0){
           return pComp;
        }
        else if(rComp != 0){
            return rComp;
        }
        else {
            return aComp;
        }
	}
	
	@Override
	public String toString() {
		return String.join("\t", data);
	}
		
	static class AlphanumericComparator implements Comparator<String> {
	    public int compare(String firstString, String secondString) {
	        if (secondString == null || firstString == null) {
	            return 0;
	        }
	 
	        int lengthFirstStr = firstString.length();
	        int lengthSecondStr = secondString.length();
	 
	        int index1 = 0;
	        int index2 = 0;
	 
	        while (index1 < lengthFirstStr && index2 < lengthSecondStr) {
	            char ch1 = firstString.charAt(index1);
	            char ch2 = secondString.charAt(index2);
	 
	            char[] space1 = new char[lengthFirstStr];
	            char[] space2 = new char[lengthSecondStr];
	 
	            int loc1 = 0;
	            int loc2 = 0;
	 
	            do {
	                space1[loc1++] = ch1;
	                index1++;
	 
	                if (index1 < lengthFirstStr) {
	                    ch1 = firstString.charAt(index1);
	                } else {
	                    break;
	                }
	            } while (Character.isDigit(ch1) == Character.isDigit(space1[0]));
	 
	            do {
	                space2[loc2++] = ch2;
	                index2++;
	 
	                if (index2 < lengthSecondStr) {
	                    ch2 = secondString.charAt(index2);
	                } else {
	                    break;
	                }
	            } while (Character.isDigit(ch2) == Character.isDigit(space2[0]));
	 
	            String str1 = new String(space1);
	            String str2 = new String(space2);
	 
	            int result;
	 
	            if (Character.isDigit(space1[0]) && Character.isDigit(space2[0])) {
	                Integer firstNumberToCompare = new Integer(Integer
	                        .parseInt(str1.trim()));
	                Integer secondNumberToCompare = new Integer(Integer
	                        .parseInt(str2.trim()));
	                result = firstNumberToCompare.compareTo(secondNumberToCompare);
	            } else {
	                result = str1.compareTo(str2);
	            }
	 
	            if (result != 0) {
	                return result;
	            }
	        }
	        return lengthFirstStr - lengthSecondStr;
	    }
	}
	
	/*
	static class VariantNameComparator implements Comparator<String> {
	    public int compare(String varName1, String varName2){
	    	String chr1 = varName1.split("\t\\|\t")[0];
	    	String chr2 = varName2.split("\t\\|\t")[0];
	    	Integer pos1 = Integer.valueOf(varName1.split("\t\\|\t")[1]);
	    	Integer pos2 = Integer.valueOf(varName2.split("\t\\|\t")[1]);
	    	String ref1 = varName1.split("\t\\|\t")[2];
	    	String ref2 = varName2.split("\t\\|\t")[2];
	    	String alt1 = varName1.split("\t\\|\t")[3];
	    	String alt2 = varName2.split("\t\\|\t")[3];
	    	
    		int cComp = alnumComparator.compare(chr1, chr2);
	    	//int cComp = chr1.compareTo(chr2);
	        int pComp = pos1.compareTo(pos2);
	        int rComp = ref1.compareTo(ref2);
	        int aComp = alt1.compareTo(alt2);
	        
	        if (cComp != 0) {
	           return cComp;
	        } 
	        else if(pComp != 0){
	           return pComp;
	        }
	        else if(rComp != 0){
	            return rComp;
	        }
	        else {
	            return aComp;
	        }
	    }
	}
	*/
}
