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

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.gkano.bioinfo.var.GeneralTools;
import com.gkano.bioinfo.var.Logger;

@Parameters(commandDescription = "FASTA2TREE")
public class UtilFASTA2TREE {

    public UtilFASTA2TREE() {
	}

	public static String getUtilName() {
		return "FASTA2DIST";
	}

	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(names = { "-v", "--verbose"})
	private boolean verbose = false;

	@Parameter(description = "<positional input files>")
    private List<String> positionalInputFiles = new ArrayList<>();

	@Parameter(names = { "-i", "--input" }, description = "FASTA input file(s)", variableArity = true)
    private List<String> namedInputFiles = new ArrayList<>();

	@Parameter(names = { "-o", "--output" }, description = "Tree output file")
    private String outputFile;
	
	@Parameter(names = {"--isFastq","-q"}, description = "Input is FASTQ (default: false)")
	private boolean isFastq = false;

	@Parameter(names = {"--kmerSize","-k"}, description = "Kmer size (default: 4)")
	private Integer k = 4;

	@Parameter(names = { "--normalize", "-n" })
	private boolean normalize = false;
	
	@Parameter(names = { "--numberOfThreads", "-t" })
	private int numOfThreads = 1;

	public void go() {
		try (PrintStream ops = GeneralTools.getPrintStreamOrExit(outputFile, this)) {
			
			
		} 
		catch (Exception e) {
			Logger.error(this, e.getMessage());
		}
	}
}
