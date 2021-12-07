/*
 *
 * BioInfoJava-Utils ciat.agrobio.bean.fingerprint.WekaAttributesMissingData
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

import weka.core.AttributeStats;
import weka.core.Instances;
import weka.core.converters.ConverterUtils.DataSource;

public class WekaAttributesMissingData {

	private String inputFileName;
	double missingData;
	
	private WekaAttributesMissingData(String inputFileName, double missingData) {
		this.inputFileName = inputFileName;
		this.missingData = missingData;
	}
	
	public void go() {
		try {
			DataSource source = new DataSource(inputFileName);
			Instances instances = source.getDataSet();
			if (instances.classIndex() == -1) instances.setClassIndex(instances.numAttributes()-1);
			System.out.println("Num instances \t"+instances.numInstances());
			System.out.println("Num attributes \t"+instances.numAttributes());
			System.out.println("Class index \t"+instances.classIndex());
			
			for(int i=0; i<instances.numAttributes(); i++) {
				AttributeStats ats = instances.attributeStats(i);
				double missing = (double)ats.missingCount / (double)instances.numInstances();
				if(missing<missingData)
					System.out.println((i+1)+"\t"+missing);
			}
			
		
			
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		WekaAttributesMissingData tool = new WekaAttributesMissingData(args[0], Double.parseDouble(args[1]));
		tool.go();
	}
}
