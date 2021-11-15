package ciat.agrobio.javautils;

import java.util.TreeMap;
import java.util.TreeSet;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import ciat.agrobio.core.GeneralTools;
import ciat.agrobio.hcluster.HierarchicalCluster;

@Parameters(commandDescription = "DIST2Clusters")
public class UtilDIST2Clusters {

	private static UtilDIST2Clusters instance = new UtilDIST2Clusters();

	private UtilDIST2Clusters() {
	}

	public static UtilDIST2Clusters getInstance() {
		return instance;
	}

	public static String getUtilName() {
		return "DIST2Clusters";
	}

	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(description = "Input_File", required = true)
	private String inputFileName;
	
	@Parameter(names = { "--numberOfThreads", "-t" })
	private int numOfThreads = 1;
	
	@Parameter(names = { "--cutHeight", "-c" })
	private Double cutHeight = null;
	
	@Parameter(names = { "--minClusterSize", "-m" })
	private Integer minClusterSize = 1;
	
	@Parameter(names={"--extra", "-e"})
	private boolean extra = false;

	@SuppressWarnings("unused")
	public void go() {
		try {
			int cpus = Runtime.getRuntime().availableProcessors();
			int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
			System.err.println("cpus=" + cpus);
			System.err.println("using=" + usingThreads);

			//Read distances matrix and sample names
			//Object[] data = readDistancesSamples(inputFileName);
			Object[] data = GeneralTools.readDistancesSamples(inputFileName);
			
			//HCluster
			//TreeMap<Integer, TreeSet<String>> clusters = hclusteringClusters((String[])data[1], (double[][])data[0], minClusterSize, cutHeight);
			TreeMap<Integer, TreeSet<String>> clusters = HierarchicalCluster.hclusteringClusters((String[])data[1], (double[][])data[0], minClusterSize, cutHeight, extra);
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
