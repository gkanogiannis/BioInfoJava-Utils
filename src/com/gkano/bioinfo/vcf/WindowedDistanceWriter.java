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

import java.io.PrintStream;

import com.gkano.bioinfo.var.GeneralTools;

public interface WindowedDistanceWriter {

    void write(String chrom, int start, int end, int nVariants,
               String[] sampleNames, double[][] distances);

    static WindowedDistanceWriter concat(PrintStream ops) {
        return new ConcatTextWriter(ops);
    }

    static WindowedDistanceWriter longTsv(PrintStream ops) {
        return new LongTsvWriter(ops);
    }

    final class ConcatTextWriter implements WindowedDistanceWriter {
        private final PrintStream ops;
        ConcatTextWriter(PrintStream ops) { this.ops = ops; }

        @Override
        public void write(String chrom, int start, int end, int nVariants,
                          String[] sampleNames, double[][] distances) {
            int N = sampleNames.length;
            ops.println("# window chrom=" + chrom + " start=" + start + " end=" + end
                    + " nvariants=" + nVariants + " nsamples=" + N);
            ops.println(N + "\t" + nVariants);
            for (int i = 0; i < N; i++) {
                ops.print(sampleNames[i]);
                for (int j = 0; j < N; j++) {
                    ops.print("\t" + GeneralTools.decimalFormat.format(distances[i][j]));
                }
                ops.println();
            }
        }
    }

    final class LongTsvWriter implements WindowedDistanceWriter {
        private final PrintStream ops;
        private boolean headerWritten = false;
        LongTsvWriter(PrintStream ops) { this.ops = ops; }

        @Override
        public void write(String chrom, int start, int end, int nVariants,
                          String[] sampleNames, double[][] distances) {
            if (!headerWritten) {
                ops.println("chrom\tstart\tend\tsample_i\tsample_j\tdist");
                headerWritten = true;
            }
            int N = sampleNames.length;
            for (int i = 0; i < N; i++) {
                for (int j = i + 1; j < N; j++) {
                    ops.println(chrom + "\t" + start + "\t" + end + "\t"
                            + sampleNames[i] + "\t" + sampleNames[j] + "\t"
                            + GeneralTools.decimalFormat.format(distances[i][j]));
                }
            }
        }
    }
}
