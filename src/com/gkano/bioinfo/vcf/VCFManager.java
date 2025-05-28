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

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.LockSupport;
import java.util.function.Function;

import com.gkano.bioinfo.var.Logger;

@SuppressWarnings("FieldMayBeFinal")
public class VCFManager<T> implements Runnable {

    private List<String> commentData;
    private String headerData;

    private AtomicInteger currVariantId = new AtomicInteger(0);

    private VariantManager<T> vm = null;
    private final Function<String, T> variantParser;
    private List<String> inputFileNames = null;
    private boolean done = false;
    private CountDownLatch startSignal = null;
    private CountDownLatch doneSignal = null;
    private boolean verbose = false;

    public VCFManager(VariantManager<T> vm, Function<String, T> variantParser, List<String> inputFileNames, CountDownLatch startSignal, CountDownLatch doneSignal, boolean verbose) {
        this.vm = vm;
        this.variantParser = variantParser;
        this.inputFileNames = inputFileNames;
        this.startSignal = startSignal;
        this.doneSignal = doneSignal;
        this.verbose = verbose;

        this.commentData = new ArrayList<>();
    }

    public final List<String> getCommentData() {
        return commentData;
    }

    public final String getHeaderData() {
        return headerData;
    }

    public boolean isDone() {
        return this.done;
    }

    @Override
    public void run() {
        try {
            Logger.setVerbose(verbose);
            done = false;

            Logger.info(this, "START READ\tStreaming:TRUE");

            VCFDecoder decoder = new VCFDecoder();
            //byte[][] line : split at tabs, byte[] is a string between tabs
            //Stringh line : a snp line from vcf
            VCFStreamingIterator<String> iterator = new VCFStreamingIterator<>(decoder, verbose, inputFileNames);
            for (String line : iterator) {
                processVariantLine(line);
            }

            Logger.info(this, "END READ");
            vm.setNumVariants(this.currVariantId.get());
            done = true;
            doneSignal.countDown();
        } catch (Exception e) {
            Logger.error(this, e.getMessage());
            System.exit(0);
        }
    }

    private void processVariantLine(String line) {
        try {
            if (line.startsWith("##")) {
                commentData.add(line);
            } else if (line.startsWith("#")) {
                headerData = line;
                vm.setNumSamples(SNPEncoder.getNumSamplesFromHeader(headerData));
            } else if (vm.getNumSamples() > 0) {
                currVariantId.incrementAndGet();
                if (vm.getPloidy() <= 0) {
                    int[] ploidy_maxAlleles = SNPEncoder.guessPloidyAndMaxAllele(line);
                    vm.setPloidy(ploidy_maxAlleles[0]);
                    vm.setMaxAlleles(ploidy_maxAlleles[1]);
                    if (vm.getPloidy() > 0) {
                        startSignal.countDown();
                    }
                }
                T parsed = variantParser.apply(line);
                while (!vm.putVariantRaw(parsed)) {
                    LockSupport.parkNanos(100_000);
                }
            }
        } catch (Exception e) {
            Logger.error(this, e.getMessage());
        }
    }

}
