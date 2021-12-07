/*
 *
 * BioInfoJava-Utils ciat.agrobio.hcluster.ClusterPoisson
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
package ciat.agrobio.hcluster;

public class ClusterPoisson {	
	private double genomeAbundance;
	private double genomeLength;

	private int lowLimit;
	private int highLimit;
	
	@SuppressWarnings("unused")
	private ClusterPoisson(){
	}
	
	/*
	private double scaleValue(double value){
		return ((newMax-newMin) / (max-min)) * (value-min) + newMin;
	}	
	*/
	
	public ClusterPoisson(double genomeAbundance, double genomeLength) {
		this.genomeAbundance = genomeAbundance;
		this.genomeLength = genomeLength;
	}

	public double getProbability(int count, ClusterPoisson[] poissons, int myIndex){
		double probability = genomeLength;//1.0;//genomeSize;
		double sum = 0.0;
		for(int i=0; i<poissons.length; i++){
			ClusterPoisson other = poissons[i];
			if(i==myIndex){
				sum += genomeLength;
				continue;
			}
			double lnSum = lnPoissonProbabilitySum(other.genomeAbundance, genomeAbundance, count, other.genomeLength);
			double tempSum = Math.exp(lnSum);
			if(Double.isInfinite(tempSum)){
				return 0.0;
			}
			sum += tempSum;
			
			//sum += Math.pow(other.genomeAbundance / genomeAbundance, count) * Math.exp(genomeAbundance - other.genomeAbundance) * other.genomeSize; 
		}
		
		return probability/sum;
	}
	
	public static double lnPoissonProbabilitySum(double a1, double a2, int count, double size){
		return (double)count * (Math.log(a1)-Math.log(a2)) + (a2-a1) + Math.log(size);
	}
	
	public double getGenomeAbundance() {
		return genomeAbundance;
	}

	public double getGenomeLength() {
		return genomeLength;
	}

	public int getLowLimit() {
		return lowLimit;
	}

	public void setLowLimit(int lowLimit) {
		this.lowLimit = lowLimit;
	}

	public int getHighLimit() {
		return highLimit;
	}

	public void setHighLimit(int highLimit) {
		this.highLimit = highLimit;
	}

}
