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
package com.gkano.bioinfo.tree;

import java.io.IOException;
import java.io.StringWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;

import com.gkano.bioinfo.var.Logger;

public class Clade extends DefaultMutableTreeNode {
	private static final long serialVersionUID = -8392803888823319176L;
	private double branchLength;
	private int bootstrapSupport = -1;

    public Clade() {
        this(null, 1);
    }

    public Clade(String nodeLabel, double branchLength) {
    	super(nodeLabel);
    	this.branchLength = branchLength;    
    }

    public double getBranchLength() {
        return branchLength;
    }

    public void setBranchLength(double distance) {
    	if(distance < 0)
    		throw new IllegalArgumentException("Branch length cannot be negative.");
        branchLength = distance;
    }
    
    public void setBootstrapSupport(int value) { this.bootstrapSupport = value; }
    public int getBootstrapSupport() { return this.bootstrapSupport; }

	public static String cladeToString(Clade node) {
		PhylipWriter writer = new PhylipWriter();
        StringWriter sw = new StringWriter();
        writer.setOutput(sw);
        try {
            writer.write(new DefaultTreeModel(node));
        } catch (IOException ex) {
            Logger.error(Clade.class, ex.getMessage());
        }
        return sw.toString().replace("\n", "");
	}

    // Extract clades from a Clade tree
	public static Set<Set<String>> extractCladesFromTree(Clade node, int totalLeaves) {
	    Set<Set<String>> clades = new HashSet<>();
	    collectClades(node, clades, totalLeaves);
	    return clades;
	}

	private static Set<String> collectClades(Clade node, Set<Set<String>> clades, int totalLeaves) {
	    Set<String> leaves = new HashSet<>();
	    if (node.isLeaf()) {
	        leaves.add(node.toString());
	    } else {
	        for (int i = 0; i < node.getChildCount(); i++) {
	            Clade child = (Clade) node.getChildAt(i);
	            leaves.addAll(collectClades(child, clades, totalLeaves));
	        }
	        // Only add non-trivial clades (not root, not singletons)
	        if (leaves.size() > 1 && leaves.size() < totalLeaves) {
	            clades.add(new HashSet<>(leaves));
	        }
	    }
	    return leaves;
	}

	public static void annotateTreeWithBootstrap(Clade originalTree, List<Clade> bootstrapTrees, int totalLeaves, int numBootstraps) {
	    Set<Set<String>> originalClades = Clade.extractCladesFromTree(originalTree, totalLeaves);

	    // Count clade frequencies in bootstraps
	    Map<Set<String>, Integer> cladeCounts = new HashMap<>();
	    for (Set<String> clade : originalClades) {
	        cladeCounts.put(clade, 0);
	    }
	    for (Clade bootTree : bootstrapTrees) {
	        Set<Set<String>> bootClades = Clade.extractCladesFromTree(bootTree, totalLeaves);
	        for (Set<String> clade : originalClades) {
	            if (bootClades.contains(clade)) {
	                cladeCounts.put(clade, cladeCounts.get(clade) + 1);
	            }
	        }
	    }

	    // Annotate the original tree
	    annotateCladeWithBootstrap(originalTree, cladeCounts, totalLeaves, numBootstraps);
	}

	public static Set<String> annotateCladeWithBootstrap(Clade node, Map<Set<String>, Integer> cladeCounts, int totalLeaves, int numBootstraps) {
	    Set<String> leaves = new HashSet<>();
	    if (node.isLeaf()) {
	        leaves.add(node.toString());
	    } else {
	        for (int i = 0; i < node.getChildCount(); i++) {
	            Clade child = (Clade) node.getChildAt(i);
	            leaves.addAll(annotateCladeWithBootstrap(child, cladeCounts, totalLeaves, numBootstraps));
	        }
	        // Only annotate non-trivial clades (not root, not singletons)
	        if (leaves.size() > 1 && leaves.size() < totalLeaves) {
	            for (Set<String> clade : cladeCounts.keySet()) {
	                if (clade.equals(leaves)) {
	                    int support = (int) Math.round(100.0 * cladeCounts.get(clade) / numBootstraps);
	                    node.setBootstrapSupport(support);
	                    break;
	                }
	            }
	        }
	    }
	    return leaves;
	}

}
