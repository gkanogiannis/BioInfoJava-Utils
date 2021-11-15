package ciat.agrobio.javautils;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import ciat.agrobio.core.CalculateDistancesD2;
import ciat.agrobio.core.SequenceProcessor;
import ciat.agrobio.io.FastaManager;
import ciat.agrobio.io.SequenceD2;

@Parameters(commandDescription = "FASTA2DIST")
public class UtilFASTA2DIST {

	private static UtilFASTA2DIST instance = new UtilFASTA2DIST();

	private UtilFASTA2DIST() {
	}

	public static UtilFASTA2DIST getInstance() {
		return instance;
	}

	public static String getUtilName() {
		return "FASTA2DIST";
	}

	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(names = {"--inputFile","-i"}, description = "Input_File", required = true)
	private List<String> inputFileNames;
	
	@Parameter(names = {"--kmerSize","-k"}, description = "Kmer size", required = true)
	private Integer k;

	@Parameter(names = { "--normalize", "-n" })
	private boolean normalize = false;
	
	@Parameter(names = { "--numberOfThreads", "-t" })
	private int numOfThreads = 1;

	@SuppressWarnings("unused")
	public void go() {
		try {
			int cpus = Runtime.getRuntime().availableProcessors();
			int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
			System.err.println("cpus=" + cpus);
			System.err.println("using=" + usingThreads);

			ConcurrentHashMap<Integer, SequenceD2> seqVectors = new ConcurrentHashMap<Integer, SequenceD2>();

			CountDownLatch startSignal = new CountDownLatch(1);
			CountDownLatch doneSignal = new CountDownLatch(usingThreads + 1);

			//System.err.println(GeneralTools.time() + " START ");

			ExecutorService pool = Executors.newFixedThreadPool(usingThreads + 1);

			Map<Integer, SequenceProcessor> sequenceProcessors = new HashMap<Integer, SequenceProcessor>();
			FastaManager frm = new FastaManager(false, inputFileNames, startSignal, doneSignal);
			pool.execute(frm);

			SequenceProcessor.resetCounters();
			// Starting threads
			for (int i = 0; i < usingThreads; i++) {
				SequenceProcessor sp = new SequenceProcessor(seqVectors, frm, k, normalize, startSignal, doneSignal);
				sequenceProcessors.put(sp.getId(), sp);
				pool.execute(sp);
			}

			doneSignal.await();
			pool.shutdown();
				
			List<String> seqNames = frm.getSequenceNames();
			List<Integer> seqIds = frm.getSequenceIds();
			
			// Calculate distances
			CalculateDistancesD2 fj = new CalculateDistancesD2();
			double[][] distances = fj.calculateDistances(usingThreads, seqVectors, seqNames, seqIds, frm);
			
			// Print data
			int seqCounter = 0;
			DecimalFormat df = new DecimalFormat("#.############"); 
			System.out.println(seqNames.size());
			for(int i=0; i<seqNames.size(); i++) {
				String seqName1 = seqNames.get(i);
				System.out.print(seqName1);
				for (int j=0;j<seqNames.size();j++) {
					String seqName2 = seqNames.get(j);
					System.out.print("\t"+df.format(distances[i][j]));
				}
				System.out.println("");
				seqCounter++;
				//if(++sampleCounter % 10 == 0) System.err.println(GeneralTools.time()+" Samples Processed : \t"+sampleCounter);
			}
			//System.err.println(GeneralTools.time()+" Samples Processed : \t"+sampleCounter);
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}

}
