package ciat.agrobio.javautils;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import ciat.agrobio.core.GeneralTools;
import ciat.agrobio.hcluster.HierarchicalCluster;

@Parameters(commandDescription = "DIST2TREE")
public class UtilDIST2TREE {

	private static UtilDIST2TREE instance = new UtilDIST2TREE();

	private UtilDIST2TREE() {
	}

	public static UtilDIST2TREE getInstance() {
		return instance;
	}

	public static String getUtilName() {
		return "DIST2TREE";
	}

	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(description = "Input_File", required = true)
	private String inputFileName;
	
	@Parameter(names = { "--numberOfThreads", "-t" })
	private int numOfThreads = 1;

	@SuppressWarnings("unused")
	public void go() {
		try {
			int cpus = Runtime.getRuntime().availableProcessors();
			int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
			System.err.println("cpus=" + cpus);
			System.err.println("using=" + usingThreads);

			//Read distances matrix and sample names
			Object[] data = GeneralTools.readDistancesSamples(inputFileName);
			
			//HCluster
			//String treeString = hclusteringTree((String[])data[1], (double[][])data[0]);
			String treeString = HierarchicalCluster.hclusteringTree((String[])data[1], (double[][])data[0]);
			System.out.println(treeString);
			 
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
