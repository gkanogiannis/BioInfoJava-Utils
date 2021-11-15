package ciat.agrobio.hcluster;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Iterator;
import java.util.List;

import javax.swing.tree.DefaultMutableTreeNode;

public class TreeDataModel {
	@SuppressWarnings("rawtypes")
	protected List data = new ArrayList();
    protected DefaultMutableTreeNode root;

    public DefaultMutableTreeNode getRoot() {
        return root;
    }

    @SuppressWarnings({ "rawtypes", "unchecked" })
	public void setRoot(DefaultMutableTreeNode root) { //, boolean isBreadthFirst)
        this.root = root;
        data.clear();

        Enumeration en = root.breadthFirstEnumeration();
        //Enumeration en = isBreadthFirst ? root.breadthFirstEnumeration()
        //                   : root.depthFirstEnumeration();
        while(en.hasMoreElements()) {
            data.add(en.nextElement());
        }
    }

    @SuppressWarnings("unchecked")
	public void add(Object item) {
        data.add(item);
    }
    
    public int size() {
        return data.size();
    }

    public Object get(int i) {
        return data.get(i);
    }

    public int indexOf(Object obj) {
        return data.indexOf(obj);
    }

    @SuppressWarnings("rawtypes")
	public Iterator iterator() {
        return data.iterator();
    }
}
