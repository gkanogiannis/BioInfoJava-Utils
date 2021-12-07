/*
 *
 * BioInfoJava-Utils ciat.agrobio.javautils.UtilARFFClassifyFromWeka
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

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import weka.classifiers.Evaluation;
import weka.classifiers.misc.InputMappedClassifier;
import weka.core.Attribute;
import weka.core.Instances;
import weka.core.converters.ConverterUtils.DataSource;

@Parameters(commandDescription = "ARFFClassifyFromWeka")
public class UtilARFFClassifyFromWeka {

	private static UtilARFFClassifyFromWeka instance = new UtilARFFClassifyFromWeka();

	private UtilARFFClassifyFromWeka() {
	}

	public static UtilARFFClassifyFromWeka getInstance() {
		return instance;
	}

	public static String getUtilName() {
		return "ARFFClassifyFromWeka";
	}

	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(description = "Input_File", required = true)
	private String inputFileName;
	
	@Parameter(names={"--libLinearWekaModel", "-m"}, description = "LibLinear_Weka_Model_File", required = true)
	private String libLinearWekaModelFileName;
	
	@Parameter(names = { "--numberOfThreads", "-t" })
	private int numOfThreads = 1;

	@SuppressWarnings("unused")
	public void go() {
		try {
			int cpus = Runtime.getRuntime().availableProcessors();
			int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
			System.err.println("cpus=" + cpus);
			System.err.println("using=" + usingThreads);

			//Read liblinear weka model
			InputMappedClassifier classifier = (InputMappedClassifier) weka.core.SerializationHelper.read(libLinearWekaModelFileName);
			
			//Read input ARFF instances
			DataSource source = new DataSource(inputFileName);
			
			Instances unlabeled = source.getDataSet();
			unlabeled.setClassIndex(unlabeled.numAttributes()-1);
			Attribute clsUnLabeled = unlabeled.classAttribute();
			
			Instances labeled = new Instances(unlabeled);
			labeled.setClassIndex(labeled.numAttributes()-1);
			Attribute clsLabeled = labeled.classAttribute();
			
			Instances lala =classifier.getModelHeader(null);
			System.out.println(lala.classAttribute());
			
			// label instances
			for (int i = 0; i < unlabeled.numInstances(); i++) {
				double clsLabel = classifier.classifyInstance(unlabeled.instance(i));
				labeled.instance(i).setClassValue(clsLabel);
				String unlabel = unlabeled.classAttribute().value((int)unlabeled.instance(i).classValue()).trim();
				String label = labeled.classAttribute().value((int)labeled.instance(i).classValue()).trim();
				System.out.println("Sample:" + (i+1) + "\t" + unlabel + "\t" + label + "\t" + unlabel.equalsIgnoreCase(label));
			}
			
			Evaluation eval = new Evaluation(unlabeled);
			eval.evaluateModel(classifier,unlabeled);
			System.out.println(eval.toSummaryString("\nResults on Train Data\n======\n", false));
			 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
