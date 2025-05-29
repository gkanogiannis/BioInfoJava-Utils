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
package com.gkano.bioinfo.vcf;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.LockSupport;

import com.gkano.bioinfo.var.GeneralTools;
import com.gkano.bioinfo.var.Logger;

public class VariantProcessor<T> implements Runnable {

    private final VCFManager<T> vcfm;
    private final CountDownLatch startSignal;
    private final CountDownLatch doneSignal;
    private boolean verbose = false;

    private static AtomicInteger processorCount = new AtomicInteger(0);
    private final int id = processorCount.getAndIncrement();

    private int[][] localDotProd;
    private int[] localNorm;

    public static void resetCounters() {
        processorCount = new AtomicInteger(0);
    }

    public VariantProcessor(
            VCFManager<T> vcfm,
            CountDownLatch startSignal,
            CountDownLatch doneSignal
    ) {
        this(vcfm, startSignal, doneSignal, false);
    }

    public VariantProcessor(
            VCFManager<T> vcfm,
            CountDownLatch startSignal,
            CountDownLatch doneSignal,
            boolean verbose
    ) {
        this.vcfm = vcfm;
        this.startSignal = startSignal;
        this.doneSignal = doneSignal;
        this.verbose = verbose;
    }

    @Override
    public void run() {
        try {
            Logger.setVerbose(verbose);
            startSignal.await();

            int numSamples = vcfm.getNumSamples();
            int ploidy = vcfm.getPloidy();
            int maxAlleles = vcfm.getMaxAlleles();
            //int bitsPerGenotype = ploidy * maxAlleles;

            localDotProd = new int[numSamples][numSamples];
            localNorm = new int[numSamples];

            T variant;
            int[][] variantEncoded;

            while (true) {
                variant = vcfm.getNextVariantRaw();
                if (variant == null) {
                    if (!vcfm.hasMoreRaw() && vcfm.isDone()) {
                        break;
                    }
                    LockSupport.parkNanos(100_000);
                    continue;
                }
                //Process variant data
                int count = vcfm.getCurrVariantCount();
                int step = GeneralTools.getAdaptiveVariantStep(count);
                if (count % step == 0 && verbose) {
                    Logger.infoCarret(this, "VariantProcessor (" + id + "):\t" + count);
                }

                try {
                    variantEncoded = SNPEncoder.encodeSNPOneHot((String) variant, ploidy, maxAlleles, vcfm.getEncodingCache(), numSamples);    
                } catch (IllegalArgumentException e) {
                    continue;
                }
                

                for (int i = 0; i < numSamples; i++) {
                    int[] di = variantEncoded[i];
                    //float norm = (float) SNPEncoder.dotProd(di, di, bitsPerGenotype);
                    int norm = 0;
                    for (int k = 0; k < di.length; k++) {
                        norm += Integer.bitCount(di[k] & di[k]);
                    }
                    localNorm[i] += norm;

                    for (int j = i; j < numSamples; j++) {
                        int[] dj = variantEncoded[j];
                        //float dotProd = (i == j) ? norm : (float) SNPEncoder.dotProd(di, variantEncoded[j], bitsPerGenotype);
                        int dotProd = 0;
                        for (int k = 0; k < di.length; k++) {
                            dotProd += Integer.bitCount(di[k] & dj[k]);
                        }
                        localDotProd[i][j] += dotProd;
                    }
                }
            }
            doneSignal.countDown();
        } catch (InterruptedException e) {
            Logger.error(this, e.getMessage());
            doneSignal.countDown();
        }
    }

    public int getId() {
        return id;
    }

    public int[][] getLocalDotProd() {
        return localDotProd;
    }

    public int[] getLocalNorm() {
        return localNorm;
    }

}
