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

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

public class SNPEncoder {

    /**
     * Encodes a SNP VCF line into one-hot genotype vectors.
     *
     * @param snpLine full tab-separated VCF line
     * @param ploidy e.g., 2 for diploid
     * @param maxNumAlleles global maximum number of alleles (for padding)
     * @return 2D array: each row is one sample's encoded genotype vector
     */
    public static int[][] encodeSNPOneHot(String snpLine, int ploidy, int maxAlleles, Map<String, int[]> encodingCache, int numSamples) {
        //String[] fields = snpLine.split("\t");
        int[][] encoded = new int[numSamples][];
        int expectedColumns = numSamples + 9;
        int fieldStart = 0, fieldEnd, col = 0;
        String[] fields = new String[expectedColumns];
        while ((fieldEnd = snpLine.indexOf('\t', fieldStart)) != -1 && col < expectedColumns) {
            fields[col++] = snpLine.substring(fieldStart, fieldEnd);
            fieldStart = fieldEnd + 1;
        }
        fields[col] = snpLine.substring(fieldStart);

        // Parse FORMAT field to locate GT index
        String formatString = fields[8];
        int gtIndex = -1;
        int fieldIdx = 0;
        int pos = 0;
        while (pos < formatString.length()) {
            int sep = formatString.indexOf(':', pos);
            if (sep == -1) sep = formatString.length();
            String token = formatString.substring(pos, sep);
            if (token.equals("GT")) {
                gtIndex = fieldIdx;
                break;
            }
            pos = sep + 1;
            fieldIdx++;
        }
        if (gtIndex == -1) {
            throw new IllegalArgumentException("FORMAT field does not contain GT");
        }

        for (int i = 0; i < numSamples; i++) {
            String sampleField = fields[i + 9];
            String[] sampleParts = sampleField.split(":");
            if (gtIndex >= sampleParts.length) {
                continue; // missing → already zero-filled
            }
            String genotype = sampleParts[gtIndex];
            int[] encoding = encodingCache.computeIfAbsent(genotype, g -> encodeGenotypeOneHot(g, ploidy, maxAlleles));
            encoded[i] = encoding;
        }

        return encoded;
    }

    public static int[] encodeGenotypeOneHot(String genotype, int ploidy, int maxAlleles) {
        //byte[] oneHot = new byte[ploidy * maxNumAlleles];
        int totalBits = ploidy * maxAlleles;
        int numInts = (totalBits + 31) / 32;
        int[] oneHot = new int[numInts];

        if (genotype == null || genotype.equals(".") || genotype.contains(".")) {
            return oneHot;
        }

        String[] alleles = genotype.split("[/|]");
        for (int i = 0; i < ploidy && i < alleles.length; i++) {
            try {
                int allele = Integer.parseInt(alleles[i]);
                if (allele >= 0 && allele < maxAlleles) {
                    int bitIndex = i * maxAlleles + allele;
                    //oneHot[i * maxAlleles + allele] = 1;
                    oneHot[bitIndex / 32] |= (1 << (bitIndex % 32));
                }
            } catch (NumberFormatException ignored) {
            }
        }

        return oneHot;
    }

    public static int[] guessPloidyAndMaxAllele(String snpLine) {
        String[] fields = snpLine.split("\t");
        if (fields.length < 10) {
            throw new IllegalArgumentException("SNP line does not contain FORMAT and sample fields.");
        }

        String[] formatFields = fields[8].split(":");
        int gtIndex = -1;
        for (int i = 0; i < formatFields.length; i++) {
            if (formatFields[i].equals("GT")) {
                gtIndex = i;
                break;
            }
        }

        if (gtIndex == -1) {
            throw new IllegalArgumentException("GT field not found in FORMAT.");
        }

        int maxPloidy = 0;
        int maxAlleleIndex = 0;

        for (int i = 9; i < fields.length; i++) {
            String[] parts = fields[i].split(":");
            if (gtIndex >= parts.length) {
                continue;
            }

            String gt = parts[gtIndex];
            if (gt == null || gt.equals(".") || gt.contains(".")) {
                continue;
            }

            String[] alleles = gt.split("[/|]");
            maxPloidy = Math.max(maxPloidy, alleles.length);

            for (String allele : alleles) {
                try {
                    int a = Integer.parseInt(allele);
                    if (a > maxAlleleIndex) {
                        maxAlleleIndex = a;
                    }
                } catch (NumberFormatException ignored) {
                }
            }
        }

        return new int[]{maxPloidy, maxAlleleIndex + 1}; // +1 to get allele count, not max index
    }

    /**
     * Extracts number of samples from the #CHROM header line of a VCF file.
     *
     * @param headerLine a VCF header line starting with "#CHROM"
     * @return number of sample fields
     */
    public static int getNumSamplesFromHeader(String headerLine) {
        if (headerLine == null || !headerLine.startsWith("#CHROM")) {
            throw new IllegalArgumentException("Line does not start with #CHROM: " + headerLine);
        }

        String[] fields = headerLine.split("\t", -1);
        return Math.max(0, fields.length - 9);
    }

    public static List<String> getSampleNamesFromHeader(String headerLine) {
        if (headerLine == null || !headerLine.startsWith("#CHROM")) {
            throw new IllegalArgumentException("Line does not start with #CHROM: " + headerLine);
        }

        String[] fields = headerLine.split("\t", -1);
        if (fields.length <= 9) {
            // No samples present
            return Collections.emptyList();
        }

        return Arrays.asList(Arrays.copyOfRange(fields, 9, fields.length));
    }

    /*
    public static final Function<String, byte[][]> StringToByteMatrixParser = line
            -> Arrays.stream(line.split("\t", -1)) // Keep empty trailing fields
                    .map(field -> field.getBytes(StandardCharsets.UTF_8))
                    .toArray(byte[][]::new);
    */

    public static final Function<String, String> StringToStringParser = line
            -> line;

   
    public static int dotProd(int[] a, int[] b, int bitsPerGenotype) {
        int sum = 0;
        for (int idx = 0; idx < bitsPerGenotype; idx++) {
            boolean isSetA = (a[idx / 32] & (1 << (idx % 32))) != 0;
            boolean isSetB = (b[idx / 32] & (1 << (idx % 32))) != 0;
            if (isSetA && isSetB) {
                sum += 1;
            }
        }
        return sum;
    }
}
