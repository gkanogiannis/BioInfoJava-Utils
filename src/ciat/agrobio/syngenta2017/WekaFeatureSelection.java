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
