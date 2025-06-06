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

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.gkano.bioinfo.tree.HierarchicalCluster;
import com.gkano.bioinfo.var.GeneralTools;
import com.gkano.bioinfo.var.Logger;
import com.gkano.bioinfo.vcf.SNPEncoder;
import com.gkano.bioinfo.vcf.VCFManager;

@SuppressWarnings("FieldMayBeFinal")
@Parameters(commandDescription = "VCF2TREE")
public class UtilVCF2TREE {

	public UtilVCF2TREE() {
	}

	public static String getUtilName() {
		return "VCF2TREE";
	}

	@SuppressWarnings("unused")
	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(names = { "-v", "--verbose"})
	private boolean verbose = false;

	@Parameter(description = "<positional input files>")
    private List<String> positionalInputFiles = new ArrayList<>();

	@Parameter(names = { "-i", "--input" }, description = "VCF input file(s)", variableArity = true)
    private List<String> namedInputFiles = new ArrayList<>();
	
	@Parameter(names = { "-o", "--output" }, description = "Tree output file")
    private String outputFile;

	@Parameter(names = { "--numberOfThreads", "-t" })
	private int numOfThreads = 1;
	
	//@SuppressWarnings("unused")
	//@Parameter(names = { "--ignoreMissing", "-m" })
	//private boolean ignoreMissing = false;
	
	//@SuppressWarnings("unused")
	//@Parameter(names={"--onlyHets", "-h"})
	//private boolean onlyHets = false;
	
	//@SuppressWarnings("unused")
	//@Parameter(names={"--ignoreHets", "-g"})
	//private boolean ignoreHets = false;

	public void go() {
		try (PrintStream ops = GeneralTools.getPrintStreamOrExit(outputFile, this)) {
			VCFManager<String> vcfm = new VCFManager<>(
                    Stream.concat(positionalInputFiles.stream(), namedInputFiles.stream()).collect(Collectors.toList()),
                    numOfThreads,
                    SNPEncoder.StringToStringParser,
                    verbose);
            vcfm.init();
			new Thread(vcfm).start();
            vcfm.awaitFinalization();

            // Calculate distances
            float[][] distances = vcfm.reduceDotProdToDistances();
			List<String> sampleNames = vcfm.getSampleNames();

			//HCluster tree
			HierarchicalCluster hc = new HierarchicalCluster(verbose);
			String treeString = hc.hclusteringTree(sampleNames.toArray(String[]::new), distances, null);
			ops.println(treeString);
			ops.close();
		} 
		catch (InterruptedException e) {
			Logger.error(this, e.getMessage());
			System.exit(1);
		}
	}
	
}
