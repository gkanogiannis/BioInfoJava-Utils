/*
 *
 * BioInfoJava-Utils ciat.agrobio.syngenta2017.WekaFeatureSelection
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
package ciat.agrobio.syngenta2017;

import java.util.Random;

import weka.attributeSelection.GreedyStepwise;
import weka.attributeSelection.WrapperSubsetEval;
import weka.classifiers.Evaluation;
import weka.classifiers.functions.SMOreg;
import weka.classifiers.meta.AttributeSelectedClassifier;
import weka.core.Instances;
import weka.core.converters.ConverterUtils.DataSource;

public class WekaFeatureSelection {

	private WekaFeatureSelection(String inputFileName) {
		this.inputFileName = inputFileName;
	}

	private String inputFileName;
	
	public void go() {
		try {
			DataSource source = new DataSource(inputFileName);
			Instances data = source.getDataSet();
			if (data.classIndex() == -1) data.setClassIndex(data.numAttributes() - 1);
			System.out.println(data.numInstances());
			System.out.println(data.numAttributes());
			//System.out.println(data.classIndex());
			
			AttributeSelectedClassifier attr_classifier = new AttributeSelectedClassifier();
			SMOreg base = new SMOreg();
			WrapperSubsetEval eval = new WrapperSubsetEval();
			GreedyStepwise search = new GreedyStepwise();
			
			attr_classifier.setDebug(true);
			base.setDebug(true);
			search.setDebuggingOutput(true);
			search.setGenerateRanking(true);
			search.setNumToSelect(1000);
			search.setNumExecutionSlots(8);
			
			eval.setClassifier(base);
			attr_classifier.setClassifier(base);
			attr_classifier.setEvaluator(eval);
			attr_classifier.setSearch(search);
			
			Evaluation evaluation = new Evaluation(data);
			evaluation.crossValidateModel(attr_classifier, data, 10, new Random(1));
			
			System.out.println(evaluation.toSummaryString());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		WekaFeatureSelection tool = new WekaFeatureSelection(args[0]);
		tool.go();
	}
}
