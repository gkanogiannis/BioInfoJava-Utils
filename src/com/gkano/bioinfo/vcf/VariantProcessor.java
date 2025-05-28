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

    private static AtomicInteger variantCount = new AtomicInteger(0);

    private static AtomicInteger taskCount = new AtomicInteger(0);

    private final int id = taskCount.getAndIncrement();
    private VariantManager<T> vm = null;
    private VCFManager<T> vcfm = null;
    private CountDownLatch startSignal = null;
    private CountDownLatch doneSignal = null;

    private int[][] localDotProd;
    private int[] localNorm;

    private boolean verbose = false;

    public static void resetCounters() {
        variantCount = new AtomicInteger(0);
        taskCount = new AtomicInteger(0);
    }

    public VariantProcessor(VariantManager<T> vm, VCFManager<T> vcfm, CountDownLatch startSignal, CountDownLatch doneSignal) {
        this(vm, vcfm, startSignal, doneSignal, false);
    }

    public VariantProcessor(VariantManager<T> vm, VCFManager<T> vcfm, CountDownLatch startSignal, CountDownLatch doneSignal, boolean verbose) {
        this.vm = vm;
        this.vcfm = vcfm;
        this.startSignal = startSignal;
        this.doneSignal = doneSignal;
        this.verbose = verbose;
    }

    public static AtomicInteger getVariantCount() {
        return variantCount;
    }

    public void setStartSignal(CountDownLatch startSignal) {
        this.startSignal = startSignal;
    }

    public void setDoneSignal(CountDownLatch doneSignal) {
        this.doneSignal = doneSignal;
    }

    public int getId() {
        return id;
    }

    @SuppressWarnings("UseSpecificCatch")
    @Override
    public void run() {
        try {
            Logger.setVerbose(verbose);
            startSignal.await();

            int numSamples = vm.getNumSamples();
            int ploidy = vm.getPloidy();
            int maxAlleles = vm.getMaxAlleles();
            //int bitsPerGenotype = ploidy * maxAlleles;

            localDotProd = new int[numSamples][numSamples];
            localNorm = new int[numSamples];

            T variant;
            int[][] variantEncoded;

            while (true) {
                variant = vm.getNextVariantRaw();
                if (variant == null) {
                    if (!vm.hasMoreRaw() && vcfm.isDone()) {
                        break;
                    }
                    LockSupport.parkNanos(100_000);
                    continue;
                }
                //Process variant data
                int count = variantCount.incrementAndGet();
                int step = GeneralTools.getAdaptiveVariantStep(count);
                if (count % step == 0 && verbose) {
                    Logger.infoCarret(this, "VariantProcessor (" + id + "):\t" + count);
                }

                variantEncoded = SNPEncoder.encodeSNPOneHot((String) variant, ploidy, maxAlleles, vm.getEncodingCache(), numSamples);

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
        } catch (Exception e) {
            Logger.error(this, e.getMessage());
            doneSignal.countDown();
        }
    }

    public int[][] getLocalDotProd() {
        return localDotProd;
    }

    public int[] getLocalNorm() {
        return localNorm;
    }

}
