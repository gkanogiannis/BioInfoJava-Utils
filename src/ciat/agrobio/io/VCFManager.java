package ciat.agrobio.io;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.collections4.list.SetUniqueList;

import ciat.agrobio.core.GeneralTools;
import ciat.agrobio.core.Variant;

public class VCFManager implements Runnable{
	private BlockingQueue<Variant> variants = null;
	private ConcurrentHashMap<Integer,Variant> staticVariants = null;
	private List<Integer> variantIds = null;
	private List<String> variantNames = null;
	
	private List<String> commentData;
	private List<String> headerData;
	
	private AtomicInteger currVariantId = new AtomicInteger(0);
	
	private String inputFileName = null;
	private boolean done = false;
	private CountDownLatch startSignal = null;
	private CountDownLatch doneSignal = null;
	
	private ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>> individualsToVariantsData;
	
	public VCFManager(String inputFileName, CountDownLatch startSignal, CountDownLatch doneSignal, ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>> individualsToVariantsData) {
		this.inputFileName = inputFileName;
		this.startSignal = startSignal;
		this.doneSignal = doneSignal;
		this.individualsToVariantsData = individualsToVariantsData;
		
		this.variants = new LinkedBlockingQueue<Variant>();
		this.staticVariants = new ConcurrentHashMap<Integer,Variant>();
		this.variantIds = SetUniqueList.setUniqueList(new ArrayList<Integer>());
		this.variantNames = SetUniqueList.setUniqueList(new ArrayList<String>());
		
		this.commentData = new ArrayList<>();
	}
	
	@SuppressWarnings("unused")
	private void clear(){
		if(variants!=null){
			variants.clear();
		}
		variants = new LinkedBlockingQueue<Variant>();
		
		if(staticVariants!=null){
			staticVariants.clear();
		}
		staticVariants = new ConcurrentHashMap<Integer,Variant>();
		
		if(variantIds!=null){
			variantIds.clear();
		}
		variantIds = SetUniqueList.setUniqueList(new ArrayList<Integer>());
		
		if(variantNames!=null){
			variantNames.clear();
		}
		variantNames = SetUniqueList.setUniqueList(new ArrayList<String>());
	}
	
	public final List<String> getCommentData() {
		return commentData;
	}
	
	public final List<String> getHeaderData() {
		return headerData;
	}
	
	public final ConcurrentHashMap<Integer,Variant> getStaticVariants(){
		return staticVariants;
	}
	
	public final List<Integer> getVariantIds() {
		return variantIds;
	}
	
	public final List<String> getVariantNames() {
		return variantNames;
	}
	
	public boolean hasMore() {
		if(!done){
			return true;
		}
		else{
			return !this.variants.isEmpty();
		}
	}
	
	private boolean putVariant(Variant variant) {
		try {
			variants.put(variant);
			staticVariants.put(variant.getVariantId(), variant);
			variantIds.add(variant.getVariantId());
			variantNames.add(variant.getVariantName());
			return true;
		} 
		catch (Exception e) {
			variants.remove(variant);
			staticVariants.remove(variant.getVariantId());
			variantIds.remove(variant.getVariantId());
			variantNames.remove(variant.getVariantName());
			return false;
		}
	}
	
	public Variant getNextVariant() {
		try {
			return variants.poll(100, TimeUnit.MILLISECONDS);
		} 
		catch (Exception e) {
			return null;
		}
	}
	
	public void run() {
		try{
			done = false;
			startSignal.countDown();
			
		    System.err.println(GeneralTools.time()+" VCFManager: START READ");
		    
		    List<File> inputFiles = new ArrayList<File>();
		    inputFiles.add(new File(inputFileName).getCanonicalFile());
		    
		    //Check if files exist
		    for(File f : inputFiles){
		    	if(!f.exists() || !f.canRead()){
		    		System.err.println("\tERROR : File "+f+"\n\tdoes not exist ot cannot be read. Exiting.");
		    		System.exit(1);
		    	}
		    }
		    
		    VCFDecoder decoder = new VCFDecoder();
		    VCFIterator<byte[]> iterator = VCFIterator.create(decoder, inputFiles);
		    for (List<byte[]> chunk : iterator) {
		    	System.err.println(GeneralTools.time()+" VCFManager: Chunk with "+chunk.size()+" lines.");
		    	for(byte[] line : chunk) {		
		    		if(line[0]=='#' && line[1]=='#') {
		    			commentData.add((new String(line)).trim());
		    			continue;
		    		}
		    		else if(line[0]=='#') {
		    			headerData = Arrays.asList((new String(line)).trim().split("\\t"));
		    			for(int i=9; i<headerData.size(); i++) {
		    				individualsToVariantsData.put(headerData.get(i), new ConcurrentHashMap<Integer, Object>());
		    			}
		    		}
		    		else {
		    			Variant variant = new Variant(currVariantId.incrementAndGet(), (new String(line)).trim().split("\\t"));
		    			if(!putVariant(variant)) {
		    				currVariantId.decrementAndGet();
		    			}
		    		}	
		    	}
		    	chunk.clear();
		    	chunk = null;
		    }
		    
			System.err.println(GeneralTools.time()+" VCFManager: END READ");
			done = true;
			doneSignal.countDown();
		}
		catch(Exception e){
			e.printStackTrace();
			System.exit(0);
		}
	}

}
