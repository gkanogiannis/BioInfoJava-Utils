/*
 *
 * BioInfoJava-Utils ciat.agrobio.bean.fingerprint.WekaBuildClassifier
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

import weka.classifiers.Evaluation;
import weka.classifiers.functions.LibLINEAR;
import weka.classifiers.misc.InputMappedClassifier;
import weka.core.Instances;
import weka.core.converters.ConverterUtils.DataSource;

public class WekaBuildClassifier {

	private String inputFileName;
	private String outputFileName;
	
	private WekaBuildClassifier(String inputFileName, String outputFileName) {
		this.inputFileName = inputFileName;
		this.outputFileName = outputFileName;
	}
	
	public void go() {
		try {
			DataSource source = new DataSource(inputFileName);
			Instances instances = source.getDataSet();
			if (instances.classIndex() == -1) instances.setClassIndex(instances.numAttributes()-1);
			System.err.println("Num instances \t"+instances.numInstances());
			System.err.println("Num attributes \t"+instances.numAttributes());
			System.err.println("Class index \t"+instances.classIndex());
			
			LibLINEAR libLinearClassifier = new LibLINEAR();
			libLinearClassifier.setNumDecimalPlaces(6);
			InputMappedClassifier classifier = new InputMappedClassifier();
			classifier.setClassifier(libLinearClassifier);
			classifier.buildClassifier(instances);
			System.err.println(classifier.toString());
			
			weka.core.SerializationHelper.write(outputFileName, classifier);
			
			Evaluation eval = new Evaluation(instances);
			eval.evaluateModel(classifier,instances);
			System.out.println(eval.toSummaryString("\nResults on Train Data\n======\n", false));
			
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		WekaBuildClassifier tool = new WekaBuildClassifier(args[0],args[1]);
		tool.go();
	}
}
