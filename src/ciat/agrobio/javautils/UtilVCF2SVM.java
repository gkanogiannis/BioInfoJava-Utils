package ciat.agrobio.javautils;

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
import ciat.agrobio.core.VariantProcessor;
import ciat.agrobio.io.VCFManager;

@Parameters(commandDescription = "VCF2SVM")
public class UtilVCF2SVM {
	
	private static UtilVCF2SVM instance = new UtilVCF2SVM();
	
	private UtilVCF2SVM() {}
	
	public static UtilVCF2SVM getInstance() {return instance;}
	
	public static String getUtilName() {
		return "VCF2SVM";
	}
	
	@Parameter(names = "--help", help = true)
	private boolean help;
	
	@Parameter(description = "Input_File", required = true)
	private String inputFileName;
	
	@Parameter(names={"--numberOfThreads", "-t"})
	private int numOfThreads = 1;
	
	public void go() {
		try {
			int cpus = Runtime.getRuntime().availableProcessors();
		    int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
		    System.err.println("cpus="+cpus);
		    System.err.println("using="+usingThreads);
		    
		    ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>> individualsToVariantsData = new ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>>();
		    
		    CountDownLatch startSignal = new CountDownLatch(1);
		    CountDownLatch doneSignal = new CountDownLatch(usingThreads+1);
		    
		    System.err.println(GeneralTools.time()+" START ");
		    
			ExecutorService pool = Executors.newFixedThreadPool(usingThreads+1);
			
			Map<Integer, VariantProcessor> variantProcessors = new  HashMap<Integer, VariantProcessor>();
			VCFManager vcfm = new VCFManager(inputFileName, startSignal, doneSignal, individualsToVariantsData);
			pool.execute(vcfm);
			
			VariantProcessor.resetCounters();
			//Starting threads
			for(int i=0; i<usingThreads; i++){
				VariantProcessor vp = new VariantProcessor(individualsToVariantsData, vcfm, startSignal, doneSignal);
				variantProcessors.put(vp.getId(), vp);
				pool.execute(vp);
			}
			
			doneSignal.await();
			pool.shutdown();
				
			List<String> individualNames = new ArrayList<String>();
			for(int i=9; i<vcfm.getHeaderData().size(); i++) {
				individualNames.add(vcfm.getHeaderData().get(i));
			}
			
			System.err.printf("\rProcessed variants : \t%8d\n",VariantProcessor.getVariantCount().get());
			
			//Print data
			for(String individualName : individualNames) {
				int varCounter = 0;
				System.out.print(individualName);
				Iterator<Integer> iterator = vcfm.getVariantIds().iterator();
			    while(iterator.hasNext()) {
			    	Integer variantId = iterator.next();
			    	String GTCodeString;
			    	Integer[] GTCode = (Integer[])individualsToVariantsData.get(individualName).get(variantId);
			    	if(GTCode!=null) {
				    	GTCodeString = String.valueOf(++varCounter)+":"+String.valueOf(GTCode[0])+"\t"+
				    				   String.valueOf(++varCounter)+":"+String.valueOf(GTCode[1])+"\t"+
				    				   String.valueOf(++varCounter)+":"+String.valueOf(GTCode[2]);
				    	System.out.print("\t"+GTCodeString);
				    }
				    else {
				    	++varCounter; ++varCounter; ++varCounter;
				    }
			    }
				System.out.println("");
			}
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
}
