/*
 *
 * BioInfoJava-Utils ciat.agrobio.bean.fingerprint.WekaFeatureSelection
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
package ciat.agrobio.bean.fingerprint;

import weka.attributeSelection.AttributeSelection;
import weka.attributeSelection.ClassifierSubsetEval;
import weka.attributeSelection.GreedyStepwise;

import weka.classifiers.functions.LibLINEAR;
import weka.core.Instances;
import weka.core.converters.ConverterUtils.DataSource;

public class WekaFeatureSelection {

	private String inputFileName;
	private int slots;
	private boolean back = false;
	
	private WekaFeatureSelection(String inputFileName, int slots, boolean back) {
		this.inputFileName = inputFileName;
		this.slots = slots;
		this.back = back;
	}
	
	public void go() {
		try {
			DataSource source = new DataSource(inputFileName);
			Instances instances = source.getDataSet();
			if (instances.classIndex() == -1) instances.setClassIndex(instances.numAttributes()-1);
			System.out.println("Num instances \t"+instances.numInstances());
			System.out.println("Num attributes \t"+instances.numAttributes());
			System.out.println("Class index \t"+instances.classIndex());
			
			AttributeSelection attsel = new AttributeSelection();
			LibLINEAR classifier = new LibLINEAR();
			classifier.setNumDecimalPlaces(6);
			ClassifierSubsetEval eval = new ClassifierSubsetEval();
			GreedyStepwise search = new GreedyStepwise();
			
			eval.setIRClassValue("class");
			eval.setClassifier(classifier);
			search.setDebuggingOutput(true);
			search.setNumExecutionSlots(slots);
			search.setSearchBackwards(back);
			attsel.setEvaluator(eval);
			attsel.setSearch(search);
			attsel.SelectAttributes(instances);
		
			
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		WekaFeatureSelection tool = new WekaFeatureSelection(args[0], Integer.parseInt(args[1]), Boolean.parseBoolean(args[2]));
		tool.go();
	}
}
