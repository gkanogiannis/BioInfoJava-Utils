/*
 *
 * BioInfoJava-Utils ciat.agrobio.javautils.UtilDIST2Clusters
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
package ciat.agrobio.javautils;

import java.io.FileInputStream;
import java.io.PrintStream;
import java.util.TreeMap;
import java.util.TreeSet;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import ciat.agrobio.core.GeneralTools;
import ciat.agrobio.hcluster.HierarchicalCluster;

@Parameters(commandDescription = "DIST2Clusters")
public class UtilDIST2Clusters {

	//private final GeneralTools gTools = GeneralTools.getInstance();
	
	public UtilDIST2Clusters() {
	}

	public static String getUtilName() {
		return "DIST2Clusters";
	}

	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(names = { "-v", "--verbose"})
	private boolean verbose = false;

	@Parameter(description = "Positional Input File")
	private String positionalInputFile;

	@Parameter(names = { "-i", "--input" }, description = "Input file (overrides positional)")
    private String namedInputFile;

	@Parameter(names = { "-o", "--output" }, description = "Tree output file")
    private String outputFile;
	
	//@Parameter(names = { "--numberOfThreads", "-t" })
	//private int numOfThreads = 1;
	
	@Parameter(names = { "--cutHeight", "-c" })
	private Double cutHeight = null;
	
	@Parameter(names = { "--minClusterSize", "-m" })
	private Integer minClusterSize = 1;
	
	@Parameter(names={"--extra", "-e"})
	private boolean extra = false;

	@SuppressWarnings("unused")
	public void go() {
		try {
			//Output PrintStream
			PrintStream ops = System.out;
			if(outputFile != null) {
				try {
					ops = new PrintStream(outputFile);
				} 
				catch (Exception e) {
					System.err.println("Error: Cannot write to " + outputFile);
                return;
				}
			}

			// Select Input File
            String inputFileName = namedInputFile!=null?namedInputFile:positionalInputFile;

			if (inputFileName == null) {
                System.err.println("Error: No input file provided.");
                return;
            }

			//int cpus = Runtime.getRuntime().availableProcessors();
			//int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
			//System.err.println("cpus=" + cpus);
			//System.err.println("using=" + usingThreads);

			//Read distances matrix and sample names
			Object[] data;
			if ("-".equals(inputFileName)) {
				// Read data from stdin
				data = GeneralTools.readDistancesSamples(System.in);
			} else {
				// Read from file
				FileInputStream fis = new FileInputStream(inputFileName);
				data = GeneralTools.readDistancesSamples(fis);
				fis.close();
			}
			
			//HCluster cluster
			HierarchicalCluster hc = new HierarchicalCluster(true);
			TreeMap<Integer, TreeSet<String>> clusters = hc.hclusteringClusters((String[])data[1], (double[][])data[0], minClusterSize, cutHeight, extra, ops);

			ops.flush();
			ops.close();
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
