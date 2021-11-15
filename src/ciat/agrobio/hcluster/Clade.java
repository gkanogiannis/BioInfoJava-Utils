package ciat.agrobio.hcluster;

import javax.swing.tree.DefaultMutableTreeNode;

public class Clade extends DefaultMutableTreeNode {
	private static final long serialVersionUID = -8392803888823319176L;
	private double branchLength;

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

}
