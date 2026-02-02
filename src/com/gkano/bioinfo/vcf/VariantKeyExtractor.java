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

/**
 * Extracts variant keys from VCF lines for embedding lookup.
 *
 * VCF format (tab-separated):
 * CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE1  SAMPLE2  ...
 * 0      1    2   3    4    5     6       7     8       9        10       ...
 */
public class VariantKeyExtractor {

    public enum KeyFormat {
        CHROM_POS,           // chr1:12345
        CHROM_POS_REF_ALT,   // chr1:12345:A:G
        VCF_ID               // Use ID column from VCF
    }

    private final KeyFormat format;

    public VariantKeyExtractor(KeyFormat format) {
        this.format = format;
    }

    /**
     * Extract variant key from a VCF data line.
     *
     * @param vcfLine Tab-separated VCF line (not header)
     * @return Variant key in the configured format, or null if parsing fails
     */
    public String extractKey(String vcfLine) {
        if (vcfLine == null || vcfLine.isEmpty() || vcfLine.startsWith("#")) {
            return null;
        }

        // Fast tab-based parsing without creating full String array
        int[] tabPositions = findFirstTabs(vcfLine, 5); // Need columns 0-4
        if (tabPositions == null) {
            return null;
        }

        switch (format) {
            case CHROM_POS:
                return extractChromPos(vcfLine, tabPositions);
            case CHROM_POS_REF_ALT:
                return extractChromPosRefAlt(vcfLine, tabPositions);
            case VCF_ID:
                return extractVcfId(vcfLine, tabPositions);
            default:
                return null;
        }
    }

    /**
     * Find positions of first n tabs in the string.
     *
     * @return Array of tab positions, or null if not enough tabs found
     */
    private int[] findFirstTabs(String s, int n) {
        int[] positions = new int[n];
        int count = 0;
        int pos = 0;

        while (count < n && pos < s.length()) {
            int tab = s.indexOf('\t', pos);
            if (tab < 0) break;
            positions[count++] = tab;
            pos = tab + 1;
        }

        return count >= n ? positions : null;
    }

    /**
     * Extract CHROM:POS format key.
     */
    private String extractChromPos(String line, int[] tabs) {
        String chrom = line.substring(0, tabs[0]);
        String pos = line.substring(tabs[0] + 1, tabs[1]);
        return chrom + ":" + pos;
    }

    /**
     * Extract CHROM:POS:REF:ALT format key.
     */
    private String extractChromPosRefAlt(String line, int[] tabs) {
        String chrom = line.substring(0, tabs[0]);
        String pos = line.substring(tabs[0] + 1, tabs[1]);
        String ref = line.substring(tabs[2] + 1, tabs[3]);
        String alt = line.substring(tabs[3] + 1, tabs[4]);
        return chrom + ":" + pos + ":" + ref + ":" + alt;
    }

    /**
     * Extract ID column value.
     */
    private String extractVcfId(String line, int[] tabs) {
        String id = line.substring(tabs[1] + 1, tabs[2]);
        // Return null for missing ID (".")
        if (".".equals(id)) {
            return null;
        }
        return id;
    }

    /**
     * Parse KeyFormat from string (case-insensitive).
     *
     * @param formatStr Format string: "CHROM_POS", "CHROM_POS_REF_ALT", or "VCF_ID"
     * @return Corresponding KeyFormat enum
     * @throws IllegalArgumentException if format string is invalid
     */
    public static KeyFormat parseFormat(String formatStr) {
        if (formatStr == null) {
            return KeyFormat.CHROM_POS_REF_ALT; // default
        }
        switch (formatStr.toUpperCase().replace("-", "_")) {
            case "CHROM_POS":
                return KeyFormat.CHROM_POS;
            case "CHROM_POS_REF_ALT":
                return KeyFormat.CHROM_POS_REF_ALT;
            case "VCF_ID":
            case "ID":
                return KeyFormat.VCF_ID;
            default:
                throw new IllegalArgumentException(
                    "Unknown variant key format: " + formatStr +
                    ". Valid options: CHROM_POS, CHROM_POS_REF_ALT, VCF_ID");
        }
    }
}
