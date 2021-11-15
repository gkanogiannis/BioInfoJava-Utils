package ciat.agrobio.javautils;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.Scanner;
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

@Parameters(commandDescription = "VCFClassifyFromTxtModel")
public class UtilVCFClassifyFromTxtModel {

	private static UtilVCFClassifyFromTxtModel instance = new UtilVCFClassifyFromTxtModel();

	private UtilVCFClassifyFromTxtModel() {
	}

	public static UtilVCFClassifyFromTxtModel getInstance() {
		return instance;
	}

	public static String getUtilName() {
		return "VCFClassifyFromTxtModel";
	}

	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(description = "Input_File", required = true)
	private String inputFileName;
	
	@Parameter(names={"--libLinearTxtModel", "-m"}, description = "LibLinear_Txt_Model_File", required = true)
	private String libLinearTxtModelFileName;
	
	@Parameter(names = { "--numberOfThreads", "-t" })
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
		
			//Read liblinear models
			Map<String, LibLinearModel> models = readLibLinearModels(libLinearTxtModelFileName);
			
			
			//Get classifications
			for(String sampleName : sampleNames) {
				Map<String, Double> sampleX = new HashMap<String, Double>();
				Iterator<Integer> variantIdIterator = vcfm.getVariantIds().iterator();
				while (variantIdIterator.hasNext()) {
					Variant variant = vcfm.getStaticVariants().get(variantIdIterator.next());
					if(variant==null)
						continue;
					Integer[] GTCode = (Integer[])samplesToVariantsData.get(sampleName).get(variant.getVariantId());
					if(GTCode==null || GTCode.length!=3)
						continue;
					else if(ignoreHets && (GTCode[1]==1))
						continue;
					else {
						sampleX.put(variant.getVariantVCFId()+"=Hom_Ref", (double)GTCode[0]);
						sampleX.put(variant.getVariantVCFId()+"=Het", (double)GTCode[1]);
						sampleX.put(variant.getVariantVCFId()+"=Hom_Alt", (double)GTCode[2]);
					}
				}
				System.out.print(sampleName+"\t:\t");
				PriorityQueue<MyEntry> queue = new PriorityQueue<MyEntry>();
				Iterator<Entry<String, LibLinearModel>> modelIterator = models.entrySet().iterator();
				while(modelIterator.hasNext()) {
					Map.Entry<String, LibLinearModel> modelEntry = (Map.Entry<String, LibLinearModel>)modelIterator.next();
					double reply = modelEntry.getValue().reply(sampleX);
					queue.add(new MyEntry(reply, modelEntry.getKey()));
					//if(reply>0.0)
						//System.out.print(modelEntry.getKey()+"("+reply+")"+",");
				}
				//System.out.print(models.get(sampleName).reply(sampleX));
				MyEntry max = queue.poll();
				System.out.print(max.toString());
				boolean yes = sampleName.equalsIgnoreCase(max.getName());
				System.out.print("\t"+yes);
				System.out.println();
			}	
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private Map<String, LibLinearModel> readLibLinearModels(String libLinearModelsFileName) {
		try {
			HashMap<String, LibLinearModel> models = new HashMap<String, LibLinearModel>();
			Scanner scanner = new Scanner(new File(libLinearModelsFileName));
			String modelName=null;
			Map<String, Double> variantCoefMap=null;
			while(scanner.hasNextLine()) {
				String line = scanner.nextLine().trim();
				if(line.isEmpty()) {
					continue;
				}
				else if(line.contains("Model") && modelName==null) {
					modelName = line.split("[\\s]+")[line.split("[\\s]+").length-1];
					variantCoefMap = new HashMap<String, Double>();
				}
				else if(line.contains("Model") && modelName!=null) {
					models.put(modelName, new LibLinearModel(modelName, variantCoefMap));
					modelName = line.split("[\\s]+")[line.split("[\\s]+").length-1];
					variantCoefMap = new HashMap<String, Double>();
				}
				else {
					String[] data = line.replaceAll("[\\s]", "").split("\\*");
					variantCoefMap.put(data[1].trim(), Double.valueOf(data[0].trim()));
				}
			}
			models.put(modelName, new LibLinearModel(modelName, variantCoefMap));
			scanner.close();
			
			return models;
		} 
		catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	private final class LibLinearModel {
		private String modelName;
		Map<String, Double> variantCoefMap;
		
		public LibLinearModel(String modelName, Map<String, Double> variantCoefMap) {
			super();
			this.modelName = modelName;
			this.variantCoefMap = variantCoefMap;
		}
		
		public double reply(Map<String, Double> sampleX) {
			double retValue = 0.0;
			for(Entry<String, Double> entryX : sampleX.entrySet()) {
				if(variantCoefMap.containsKey(entryX.getKey())) {
					retValue += variantCoefMap.get(entryX.getKey()) * entryX.getValue();
				}
			}
			retValue += variantCoefMap.get("1.0") * 1.0;
			return retValue;
		}
		
		public String toString() {
			return modelName;
		}
	}
	
	private final class MyEntry implements Comparable<MyEntry>{
		private double value;
		private String name;
		
		public MyEntry(double value, String name) {
			super();
			this.value = value;
			this.name = name;
		}

		public double getValue() {
			return value;
		}

		@SuppressWarnings("unused")
		public void setValue(double value) {
			this.value = value;
		}

		public String getName() {
			return name;
		}

		@SuppressWarnings("unused")
		public void setName(String name) {
			this.name = name;
		}

		@Override
		public int compareTo(MyEntry o) {
			if(value<o.getValue())
				return 1;
			else if(value>o.getValue())
				return -1;
			else
			return 0;
		}
		
		public String toString() {
			return name+"("+value+")";
		}
		
	}
	
}
