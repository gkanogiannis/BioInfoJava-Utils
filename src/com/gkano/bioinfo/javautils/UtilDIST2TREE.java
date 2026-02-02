/*
 *
 * BioInfoJava-Utils 
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
package com.gkano.bioinfo.javautils;

import java.io.FileInputStream;
import java.io.PrintStream;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.gkano.bioinfo.tree.HierarchicalCluster;
import com.gkano.bioinfo.var.GeneralTools;
import com.gkano.bioinfo.var.Logger;

@SuppressWarnings("FieldMayBeFinal")
@Parameters(commandDescription = "DIST2TREE")
public class UtilDIST2TREE {

	//private final GeneralTools gTools = GeneralTools.getInstance();

	public  UtilDIST2TREE() {
	}

	public static String getUtilName() {
		return "DIST2TREE";
	}

	@SuppressWarnings("unused")
	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(names = { "-v", "--verbose"})
	private boolean verbose = false;

	@Parameter(description = "<positional input file>")
	private String positionalInputFile;

	@Parameter(names = { "-i", "--input" }, description = "Input file (overrides positional)")
    private String namedInputFile;

	@Parameter(names = { "-o", "--output" }, description = "Tree output file")
    private String outputFile;
	
	//@Parameter(names = { "--numberOfThreads", "-t" })
	//private int numOfThreads = 1;

	public void go() {
		try (PrintStream ops = GeneralTools.getPrintStreamOrExit(outputFile, this)) {

			// Select Input File
            String inputFileName = namedInputFile!=null?namedInputFile:positionalInputFile;

			if (inputFileName == null) {
				Logger.error(this, "No input file provided.");
				return;
            }

			//int cpus = Runtime.getRuntime().availableProcessors();
			//int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
			//if(verbose) System.err.println("cpus=" + cpus);
			//if(verbose) System.err.println("using=" + usingThreads);

			//Read distances matrix and sample names
			Object[] data;
			if ("-".equals(inputFileName)) {
				// Read data from stdin
				data = GeneralTools.readDistancesSamples(System.in);
			} else {
            	try ( // Read from file
                	FileInputStream fis = new FileInputStream(inputFileName)) {
                		data = GeneralTools.readDistancesSamples(fis);
                    }
			}
			
			//HCluster tree
			HierarchicalCluster hc = new HierarchicalCluster(verbose);
			String treeString = (String) hc.hclusteringTree((String[])data[1], (double[][])data[0], null)[0];
			ops.println(treeString);
			ops.close();
		} 
		catch (Exception e) {
			Logger.error(this, e.getMessage());
		}
	}
	
}
