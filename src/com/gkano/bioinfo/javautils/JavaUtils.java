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

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Date;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.MissingCommandException;
import com.beust.jcommander.ParameterException;
import com.gkano.bioinfo.var.GeneralTools;
import com.gkano.bioinfo.var.Logger;

public class JavaUtils {

    //private final GeneralTools gTools = GeneralTools.getInstance();
    public void go(String[] args) {
        JCommander jc = JCommander.newBuilder()
                //jc.addCommand(UtilVCF2ARFF.getUtilName(), new UtilVCF2ARFF());
                //jc.addCommand(UtilVCF2CSV.getUtilName(),  new UtilVCF2CSV());
                //jc.addCommand(UtilVCF2TABLE.getUtilName(),  new UtilVCF2TABLE());
                .addCommand(UtilVCF2ISTATS.getUtilName(), new UtilVCF2ISTATS()) // in R
                //jc.addCommand(UtilVCF2SVM.getUtilName(), new UtilVCF2SVM());
                .addCommand(UtilVCF2DIST.getUtilName(), new UtilVCF2DIST()) // in R
                .addCommand(UtilVCF2TREE.getUtilName(), new UtilVCF2TREE()) // in R
                //jc.addCommand(UtilVCFRemoveClones.getUtilName(), new UtilVCFRemoveClones()); //
                //jc.addCommand(UtilVCFNonrelated.getUtilName(), new UtilVCFNonrelated()); //
                //jc.addCommand(UtilVCFKeepVariants.getUtilName(), new UtilVCFKeepVariants());
                //jc.addCommand(UtilVCFKeepSamples.getUtilName(), new UtilVCFKeepSamples());
                //jc.addCommand(UtilVCFRemoveSamples.getUtilName(), new UtilVCFRemoveSamples());
                //jc.addCommand(UtilVCFClassifyFromTxtModel.getUtilName(), new UtilVCFClassifyFromTxtModel());
                //jc.addCommand(UtilVCFsIntersection.getUtilName(), new UtilVCFsIntersection()); //
                //jc.addCommand(UtilVCFFilter.getUtilName(), new UtilVCFFilter());

                .addCommand(UtilFASTA2DIST.getUtilName(), new UtilFASTA2DIST()) //

                .addCommand(UtilDIST2TREE.getUtilName(), new UtilDIST2TREE()) // in R
                .addCommand(UtilDIST2Clusters.getUtilName(), new UtilDIST2Clusters()) // in R
                //jc.addCommand(UtilDIST2Hist.getUtilName(), new UtilDIST2Hist()); // in R

                //jc.addCommand(UtilARFFClassifyFromWeka.getUtilName(), new UtilARFFClassifyFromWeka());
                //jc.addCommand(UtilARFFaddTrait.getUtilName(), new UtilARFFaddTrait());

                .build();
        try {
            jc.parse(args);

            String command = jc.getParsedCommand();
            JCommander jCommand = jc.getCommands().get(command);

            Class<?> utilClass = Class.forName(JavaUtils.class.getPackage().getName() + ".Util" + command);
            //Method getUtilInstanceMethod = utilClass.getMethod("getInstance");
            Method goMethod = utilClass.getMethod("go");
            Object util = jCommand.getObjects().get(0); // getUtilInstanceMethod.invoke(null);

            Field f = utilClass.getDeclaredField("help");
            f.setAccessible(true);
            boolean help = f.getBoolean(util);

            if (help) {
                jCommand.usage();
                //StringBuilder sb = new StringBuilder();
                //jc.usage(jc.getParsedCommand(),sb);
                //System.err.println(sb.toString());
                System.exit(0);
            }

            System.err.println(command + " : " + new Date(GeneralTools.classBuildTimeMillis(utilClass)).toString());
            goMethod.invoke(util);
        } catch (MissingCommandException | ClassNotFoundException e) {
            System.err.println("Invalid JavaUtil selection!");
            System.err.println("Use one of : " + jc.getCommands().keySet());
        } catch (ParameterException e) {
            Logger.error(this, e.getMessage());
            jc.usage();
        } catch (IllegalAccessException | IllegalArgumentException | NoSuchFieldException | NoSuchMethodException | SecurityException | InvocationTargetException e) {
            Logger.error(this, e.getMessage());
        }
    }

    public static void main(String[] args) {
        JavaUtils ju = new JavaUtils();
        ju.go(args);
    }
}
