package ciat.agrobio.javautils;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import ciat.agrobio.core.GeneralTools;
import ciat.agrobio.core.Variant;
import ciat.agrobio.core.VariantProcessor;
import ciat.agrobio.io.VCFManager;

@Parameters(commandDescription = "VCF2ARFF")
public class UtilVCF2ARFF {
	
private static UtilVCF2ARFF instance = new UtilVCF2ARFF();
	
	private UtilVCF2ARFF() {}
	
	public static UtilVCF2ARFF getInstance() {return instance;}
	
	public static String getUtilName() {
		return "VCF2ARFF";
	}
	
	@Parameter(names = "--help", help = true)
	private boolean help;
	
	@Parameter(description = "Input_File", required = true)
	private String inputFileName;
	
	@Parameter(names={"--numberOfThreads", "-t"})
	private int numOfThreads = 1;
	
	@Parameter(names={"--ignoreHets", "-g"})
	private boolean ignoreHets = false;
	
	private static VCFManager vcfm;
	
	public void go() {
		try {
			int cpus = Runtime.getRuntime().availableProcessors();
			int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
			System.err.println("cpus=" + cpus);
			System.err.println("using=" + usingThreads);

			ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>> samplesToVariantsData = new ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>>();

			CountDownLatch startSignal = new CountDownLatch(1);
			CountDownLatch doneSignal = new CountDownLatch(usingThreads + 1);

			System.err.println(GeneralTools.time() + " START ");

			ExecutorService pool = Executors.newFixedThreadPool(usingThreads + 1);

			Map<Integer, VariantProcessor> variantProcessors = new HashMap<Integer, VariantProcessor>();
			vcfm = new VCFManager(inputFileName, startSignal, doneSignal, samplesToVariantsData);
			pool.execute(vcfm);

			VariantProcessor.resetCounters();
			// Starting threads
			for (int i = 0; i < usingThreads; i++) {
				VariantProcessor vp = new VariantProcessor(samplesToVariantsData, vcfm, startSignal, doneSignal);
				variantProcessors.put(vp.getId(), vp);
				pool.execute(vp);
			}

			doneSignal.await();
			pool.shutdown();

			List<String> sampleNames = new ArrayList<String>();
			for (int i = 9; i < vcfm.getHeaderData().size(); i++) {
				sampleNames.add(vcfm.getHeaderData().get(i));
			}
			
			System.err.printf("\rProcessed variants : \t%8d\n", VariantProcessor.getVariantCount().get());
			
			//Print ARFF file
			//Print attributes (variants)
			System.out.println("@relation\t"+(new File(inputFileName).getName())+"\n");
			Iterator<Integer> variantIdIterator = vcfm.getVariantIds().iterator();
			while (variantIdIterator.hasNext()) {
				Variant variant = vcfm.getStaticVariants().get(variantIdIterator.next());
				if(variant==null)
					continue;
		    	System.out.println("@attribute\t" + variant.getVariantVCFId() + "=Hom_Ref" + "\tnumeric");
			    System.out.println("@attribute\t" + variant.getVariantVCFId() + "=Het" + "\tnumeric");
			    System.out.println("@attribute\t" + variant.getVariantVCFId() + "=Hom_Alt" + "\tnumeric");
		    }
		    System.out.println("@attribute\tclass\t"+"{"+String.join(",", sampleNames)+"}");
		    
			//Print data
			System.out.println("\n@data");
			for(String sampleName : sampleNames) {
				variantIdIterator = vcfm.getVariantIds().iterator();
			    while(variantIdIterator.hasNext()) {
			    	Integer variantId = variantIdIterator.next();
			    	String GTCodeString;
			    	Integer[] GTCode = (Integer[])samplesToVariantsData.get(sampleName).get(variantId);
			    	if(GTCode==null || GTCode.length==0)
			    		GTCodeString = "?,?,?";
			    	else if(ignoreHets && (GTCode[1]==1))
			    		GTCodeString = "?,?,?";
			    	else
			    		GTCodeString = String.valueOf(GTCode[0])+","+String.valueOf(GTCode[1])+","+String.valueOf(GTCode[2]);
			    	System.out.print(GTCodeString+",");
			    }
				System.out.println(sampleName);
			}
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
}
