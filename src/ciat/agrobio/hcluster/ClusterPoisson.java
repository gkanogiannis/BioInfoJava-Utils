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
