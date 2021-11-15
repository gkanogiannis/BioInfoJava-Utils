package ciat.agrobio.core;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;
import java.util.concurrent.RecursiveTask;
import java.util.concurrent.atomic.AtomicInteger;

import ciat.agrobio.io.VCFManager;

public class CalculateDistancesCOSINE {

	public CalculateDistancesCOSINE() {
	}
	
	public double[][] calculateDistances(int numOfThreads, ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>> samplesToVariantsData, List<String> sampleNames, VCFManager vcfm, boolean ignoreHets, boolean onlyHets, boolean ignoremissing) {
		try {
			int numOfSamples = sampleNames.size();
			double[][] distances = new double[numOfSamples][numOfSamples]; 
			CalculateDistancesTask task = new CalculateDistancesTask(samplesToVariantsData, sampleNames, vcfm, ignoreHets, onlyHets, ignoremissing, distances);
			ForkJoinPool pool = new ForkJoinPool(numOfThreads);
			pool.execute(task);
			
			return task.get();
		} 
		catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
}

@SuppressWarnings("serial")
class CalculateDistancesTask extends RecursiveTask<double[][]> {
	private ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>> samplesToVariantsData;
	private List<String> sampleNames;
	private VCFManager vcfm;
	private boolean ignoreHets, onlyHets, ignoremissing;
	
	private double[][] distances;
	
	public CalculateDistancesTask(ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>> samplesToVariantsData, 
			   					  List<String> sampleNames, VCFManager vcfm, boolean ignoreHets, boolean onlyHets, boolean ignoremissing, double[][] distances) {
		
		this.samplesToVariantsData = samplesToVariantsData;
		this.sampleNames = sampleNames; this.vcfm = vcfm; this.ignoreHets = ignoreHets;
		this.onlyHets = onlyHets; this.ignoremissing = ignoremissing; this.distances = distances;
	}

	@SuppressWarnings("rawtypes")
	protected double[][] compute() {
		List<ForkJoinTask> children = new ArrayList<ForkJoinTask>();
		for ( int row=sampleNames.size()-1; row>=0; row--) {
			children.add(new CalculateDistancesChildTask(samplesToVariantsData, sampleNames, vcfm, ignoreHets, onlyHets, ignoremissing, distances, row));
		}
		invokeAll(children);
		return distances;
	}
	
}

@SuppressWarnings("serial")
class CalculateDistancesChildTask extends RecursiveAction {
	private static AtomicInteger sampleCounter = new AtomicInteger(0);
	private static AtomicInteger taskCount = new AtomicInteger(0);
	private final int id = taskCount.getAndIncrement();
	
	private ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>> samplesToVariantsData;
	private List<String> sampleNames;
	private VCFManager vcfm;
	private boolean ignoreHets, onlyHets, ignoremissing;
	
	private double[][] distances;
	private int row;
	
	public CalculateDistancesChildTask(ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>> samplesToVariantsData, List<String> sampleNames, 
									   VCFManager vcfm, boolean ignoreHets, boolean onlyHets, boolean ignoremissing, double[][] distances, int row) {
		this.samplesToVariantsData = samplesToVariantsData;
		this.sampleNames = sampleNames; this.vcfm = vcfm; this.ignoreHets = ignoreHets;
		this.onlyHets = onlyHets; this.ignoremissing = ignoremissing; this.distances = distances; this.row = row;
	}

	@Override
	protected void compute() {
		DecimalFormat df = new DecimalFormat("#.############"); 
		try {
			String sampleName1 = sampleNames.get(row);
			//System.err.print(sampleName1);
			for (int column=sampleNames.size()-1;column>=row;column--) {
				String sampleName2 = sampleNames.get(column);
				double dot = 0.0;
				double notmissing = 0.0;
				double norm1 = 0.0;
				double norm2 = 0.0;
				Iterator<Integer> iterator = vcfm.getVariantIds().iterator();
				while (iterator.hasNext()) {
					Integer variantId = iterator.next();
					Integer[] GTCode1 = (Integer[])samplesToVariantsData.get(sampleName1).get(variantId);
					Integer[] GTCode2 = (Integer[])samplesToVariantsData.get(sampleName2).get(variantId);
					
					if((GTCode1==null && GTCode2==null) || (GTCode1.length==0 && GTCode2.length==0))
						continue;
					else if(ignoreHets && ((GTCode1!=null&&GTCode1.length==3&&GTCode1[1]==1) || (GTCode2!=null&&GTCode2.length==3&&GTCode2[1]==1)))
						continue;
					else if(onlyHets && ((GTCode1!=null&&GTCode1.length==3&&GTCode1[1]!=1) && (GTCode2!=null&&GTCode2.length==3&&GTCode2[1]!=1)))
						continue;
					else if((GTCode1==null || GTCode1.length==0) && (GTCode2!=null&&GTCode2.length==3))
						norm2 += ( (GTCode2[0]*GTCode2[0]) + (GTCode2[1]*GTCode2[1]) + (GTCode2[2]*GTCode2[2]) );
					else if((GTCode2==null || GTCode2.length==0) && (GTCode1!=null&&GTCode1.length==3))
						norm1 += ( (GTCode1[0]*GTCode1[0]) + (GTCode1[1]*GTCode1[1]) + (GTCode1[2]*GTCode1[2]) );
					else {
						dot += ( (GTCode1[0]*GTCode2[0]) + (GTCode1[1]*GTCode2[1]) + (GTCode1[2]*GTCode2[2]) );
						norm1 += ( (GTCode1[0]*GTCode1[0]) + (GTCode1[1]*GTCode1[1]) + (GTCode1[2]*GTCode1[2]) );
						norm2 += ( (GTCode2[0]*GTCode2[0]) + (GTCode2[1]*GTCode2[1]) + (GTCode2[2]*GTCode2[2]) );
						notmissing++;
					}
				}
				if(ignoremissing) {	
					double cosine = dot / notmissing;
					if(Double.isInfinite(cosine) || Double.isNaN(cosine)) cosine = 0.0;
					
					double distance = 1.0 - cosine;
					if(distance<=0.0) distance = 0.0;
					
					distances[row][column] = Double.parseDouble(df.format(distance));
					distances[column][row] = distances[row][column];
				}
				else {
					double cosine = dot / (Math.sqrt(norm1)*Math.sqrt(norm2));
					if(Double.isInfinite(cosine) || Double.isNaN(cosine)) cosine = 0.0;
					
					double distance = 1.0 - cosine;
					if(distance<=0.0) distance = 0.0;
					
					//double distance = -1.0*Math.log(cosine);
					//if(distance<=0.0) distance = 0.0;
					
					distances[row][column] = Double.parseDouble(df.format(distance));
					distances[column][row] = distances[row][column];
				}
			}
			int count = sampleCounter.incrementAndGet();
			if(count % 50 == 0) System.err.println(GeneralTools.time()+" CalculateDistancesChildTask("+id+"): "+"\t"+count);
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
