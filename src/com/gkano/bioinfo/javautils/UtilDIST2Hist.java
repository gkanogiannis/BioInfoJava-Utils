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

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.ArrayList;

import javax.imageio.ImageIO;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.gkano.bioinfo.var.GeneralTools;
import com.gkano.bioinfo.var.Logger;

@SuppressWarnings("FieldMayBeFinal")
@Parameters(commandDescription = "DIST2Hist")
public class UtilDIST2Hist {

    private static UtilDIST2Hist instance = new UtilDIST2Hist();

    private UtilDIST2Hist() {
    }

    public static UtilDIST2Hist getInstance() {
        return instance;
    }

    public static String getUtilName() {
        return "DIST2Hist";
    }

    @SuppressWarnings("unused")
    @Parameter(names = "--help", help = true)
    private boolean help;

    @Parameter(description = "<input file>", required = true)
    private String inputFileName;

    @Parameter(names = {"--output", "-o"}, description = "Output png image file", required = true)
    private String outputFileName = null;

    @Parameter(names = {"--settings", "-s"}, description = "Settings for bins,width,height")
    private String settings = "100,1024,768";

    public void go() {
        try {
            //Read distances matrix and sample names
            //Object[] data = readDistancesSamples(inputFileName);
            Object[] data = GeneralTools.readDistancesSamples(inputFileName);

            //Histogram
            JFreeChart hist = createHistogram((String[]) data[1], (double[][]) data[0], settings);

            //Create output png image
            BufferedImage bImage = hist.createBufferedImage(Integer.parseInt(settings.split(",")[1]),
                    Integer.parseInt(settings.split(",")[2]), null);
            try (OutputStream os = new FileOutputStream(outputFileName)) {
                ImageIO.write(bImage, "png", os);
            }
        } catch (Exception e) {
            Logger.error(this, e.getMessage());
        }
    }

    private JFreeChart createHistogram(String[] sampleNames, double[][] distances, String settings) {
        try {
            ArrayList<Double> dataList = new ArrayList<>();
            for (int i = 0; i < sampleNames.length; i++) {
                for (int j = i; j < sampleNames.length; j++) {
                    if (!Double.isNaN(distances[i][j])) {
                        dataList.add(distances[i][j]);
                    }
                }
            }

            double[] data = dataList.stream().mapToDouble(Double::doubleValue).toArray();

            HistogramDataset dataset = new HistogramDataset();
            dataset.setType(HistogramType.FREQUENCY);
            dataset.addSeries("Hist", data, Integer.parseInt(settings.split(",")[0])); // Number of bins is settings[0]

            String plotTitle = "";
            String xAxis = "Distance";
            String yAxis = "Density";
            JFreeChart chart = ChartFactory.createHistogram(plotTitle, xAxis, yAxis, dataset,
                    PlotOrientation.VERTICAL, false, false, false);
            chart.setBackgroundPaint(Color.white);

            return chart;
        } catch (NumberFormatException e) {
            Logger.error(this, e.getMessage());
            return null;
        }
    }

}
