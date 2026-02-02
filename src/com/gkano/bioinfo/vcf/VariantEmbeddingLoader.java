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

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import com.gkano.bioinfo.var.Logger;

/**
 * Loads variant embeddings from various file formats.
 * Supports TSV and HuggingFace JSON formats.
 */
public class VariantEmbeddingLoader {

    public enum EmbeddingFormat {
        TSV,           // Tab-separated: VARIANT_ID\tDIM_0\tDIM_1\t...
        HUGGINGFACE    // JSON format from HuggingFace models
    }

    private VariantEmbeddingLoader() {
        // Utility class
    }

    /**
     * Load embeddings from file, auto-detecting format by extension.
     *
     * @param filePath Path to embeddings file (.tsv, .txt, .json, or .gz variants)
     * @return Map from variant key (CHROM:POS:REF:ALT) to embedding vector
     * @throws IOException if file cannot be read
     */
    public static Map<String, double[]> loadEmbeddings(String filePath) throws IOException {
        EmbeddingFormat format = detectFormat(filePath);
        return loadEmbeddings(filePath, format);
    }

    /**
     * Load embeddings from file with specified format.
     *
     * @param filePath Path to embeddings file
     * @param format   File format (TSV or HUGGINGFACE)
     * @return Map from variant key to embedding vector
     * @throws IOException if file cannot be read
     */
    public static Map<String, double[]> loadEmbeddings(String filePath, EmbeddingFormat format) throws IOException {
        switch (format) {
            case TSV:
                return loadTSV(filePath);
            case HUGGINGFACE:
                return loadHuggingFace(filePath);
            default:
                throw new IllegalArgumentException("Unsupported embedding format: " + format);
        }
    }

    /**
     * Detect format from file extension.
     */
    private static EmbeddingFormat detectFormat(String filePath) {
        String lower = filePath.toLowerCase();
        // Remove compression extension first
        if (lower.endsWith(".gz")) {
            lower = lower.substring(0, lower.length() - 3);
        }
        if (lower.endsWith(".json")) {
            return EmbeddingFormat.HUGGINGFACE;
        }
        // Default to TSV for .tsv, .txt, or unknown
        return EmbeddingFormat.TSV;
    }

    /**
     * Load embeddings from TSV format.
     *
     * Format:
     * #VARIANT_ID    DIM_0    DIM_1    DIM_2    ...    (optional header)
     * chr1:12345:A:G    0.123    -0.456    0.789    ...
     * chr1:67890:C:T    0.567    0.123    -0.890    ...
     *
     * Lines starting with # are treated as comments/headers.
     */
    private static Map<String, double[]> loadTSV(String filePath) throws IOException {
        Map<String, double[]> embeddings = new HashMap<>();
        int embeddingDim = -1;
        int lineNum = 0;

        try (BufferedReader reader = createReader(filePath)) {
            String line;
            while ((line = reader.readLine()) != null) {
                lineNum++;
                line = line.trim();

                // Skip empty lines and comments/headers
                if (line.isEmpty() || line.startsWith("#")) {
                    continue;
                }

                String[] parts = line.split("\t");
                if (parts.length < 2) {
                    Logger.warn(VariantEmbeddingLoader.class,
                        "Skipping malformed line " + lineNum + ": insufficient columns");
                    continue;
                }

                String variantKey = parts[0].trim();
                int dim = parts.length - 1;

                // Validate consistent dimension
                if (embeddingDim < 0) {
                    embeddingDim = dim;
                } else if (dim != embeddingDim) {
                    Logger.warn(VariantEmbeddingLoader.class,
                        "Skipping line " + lineNum + ": dimension mismatch (expected " +
                        embeddingDim + ", got " + dim + ")");
                    continue;
                }

                // Parse embedding values
                double[] embedding = new double[dim];
                boolean valid = true;
                for (int i = 0; i < dim; i++) {
                    try {
                        embedding[i] = Double.parseDouble(parts[i + 1].trim());
                    } catch (NumberFormatException e) {
                        Logger.warn(VariantEmbeddingLoader.class,
                            "Skipping line " + lineNum + ": invalid number at column " + (i + 2));
                        valid = false;
                        break;
                    }
                }

                if (valid) {
                    embeddings.put(variantKey, embedding);
                }
            }
        }

        Logger.info(VariantEmbeddingLoader.class,
            "Loaded " + embeddings.size() + " embeddings with dimension " + embeddingDim);
        return embeddings;
    }

    /**
     * Load embeddings from HuggingFace JSON format.
     *
     * Expected format:
     * {
     *   "model_name": "...",           (optional)
     *   "embedding_dim": 768,          (optional)
     *   "variants": [
     *     {"id": "chr1:12345:A:G", "embedding": [0.123, -0.456, ...]},
     *     {"id": "chr1:67890:C:T", "embedding": [0.567, 0.123, ...]},
     *     ...
     *   ]
     * }
     *
     * Alternative flat format:
     * {
     *   "chr1:12345:A:G": [0.123, -0.456, ...],
     *   "chr1:67890:C:T": [0.567, 0.123, ...],
     *   ...
     * }
     */
    private static Map<String, double[]> loadHuggingFace(String filePath) throws IOException {
        Map<String, double[]> embeddings = new HashMap<>();

        // Read entire file content
        StringBuilder content = new StringBuilder();
        try (BufferedReader reader = createReader(filePath)) {
            String line;
            while ((line = reader.readLine()) != null) {
                content.append(line);
            }
        }

        String json = content.toString().trim();

        // Simple JSON parsing without external dependencies
        // Detect format by looking for "variants" key
        if (json.contains("\"variants\"")) {
            parseHuggingFaceNested(json, embeddings);
        } else {
            parseHuggingFaceFlat(json, embeddings);
        }

        if (!embeddings.isEmpty()) {
            int dim = embeddings.values().iterator().next().length;
            Logger.info(VariantEmbeddingLoader.class,
                "Loaded " + embeddings.size() + " embeddings with dimension " + dim);
        }

        return embeddings;
    }

    /**
     * Parse nested HuggingFace format with "variants" array.
     */
    private static void parseHuggingFaceNested(String json, Map<String, double[]> embeddings) {
        // Find variants array
        int variantsStart = json.indexOf("\"variants\"");
        if (variantsStart < 0) return;

        int arrayStart = json.indexOf('[', variantsStart);
        if (arrayStart < 0) return;

        int arrayEnd = findMatchingBracket(json, arrayStart, '[', ']');
        if (arrayEnd < 0) return;

        String variantsArray = json.substring(arrayStart + 1, arrayEnd);

        // Parse each variant object
        int pos = 0;
        while (pos < variantsArray.length()) {
            int objStart = variantsArray.indexOf('{', pos);
            if (objStart < 0) break;

            int objEnd = findMatchingBracket(variantsArray, objStart, '{', '}');
            if (objEnd < 0) break;

            String obj = variantsArray.substring(objStart, objEnd + 1);
            parseVariantObject(obj, embeddings);

            pos = objEnd + 1;
        }
    }

    /**
     * Parse a single variant object: {"id": "...", "embedding": [...]}
     */
    private static void parseVariantObject(String obj, Map<String, double[]> embeddings) {
        String id = extractJsonString(obj, "id");
        if (id == null) return;

        double[] embedding = extractJsonArray(obj, "embedding");
        if (embedding == null) return;

        embeddings.put(id, embedding);
    }

    /**
     * Parse flat HuggingFace format: {"variant_key": [embedding], ...}
     */
    private static void parseHuggingFaceFlat(String json, Map<String, double[]> embeddings) {
        // Remove outer braces
        json = json.trim();
        if (json.startsWith("{")) json = json.substring(1);
        if (json.endsWith("}")) json = json.substring(0, json.length() - 1);

        int pos = 0;
        while (pos < json.length()) {
            // Find key
            int keyStart = json.indexOf('"', pos);
            if (keyStart < 0) break;

            int keyEnd = json.indexOf('"', keyStart + 1);
            if (keyEnd < 0) break;

            String key = json.substring(keyStart + 1, keyEnd);

            // Find array value
            int arrayStart = json.indexOf('[', keyEnd);
            if (arrayStart < 0) break;

            int arrayEnd = findMatchingBracket(json, arrayStart, '[', ']');
            if (arrayEnd < 0) break;

            String arrayStr = json.substring(arrayStart + 1, arrayEnd);
            double[] embedding = parseDoubleArray(arrayStr);

            if (embedding != null && embedding.length > 0) {
                embeddings.put(key, embedding);
            }

            pos = arrayEnd + 1;
        }
    }

    /**
     * Extract a string value from JSON object.
     */
    private static String extractJsonString(String json, String key) {
        String search = "\"" + key + "\"";
        int keyPos = json.indexOf(search);
        if (keyPos < 0) return null;

        int colonPos = json.indexOf(':', keyPos + search.length());
        if (colonPos < 0) return null;

        int valueStart = json.indexOf('"', colonPos);
        if (valueStart < 0) return null;

        int valueEnd = json.indexOf('"', valueStart + 1);
        if (valueEnd < 0) return null;

        return json.substring(valueStart + 1, valueEnd);
    }

    /**
     * Extract a double array from JSON object.
     */
    private static double[] extractJsonArray(String json, String key) {
        String search = "\"" + key + "\"";
        int keyPos = json.indexOf(search);
        if (keyPos < 0) return null;

        int arrayStart = json.indexOf('[', keyPos);
        if (arrayStart < 0) return null;

        int arrayEnd = findMatchingBracket(json, arrayStart, '[', ']');
        if (arrayEnd < 0) return null;

        String arrayStr = json.substring(arrayStart + 1, arrayEnd);
        return parseDoubleArray(arrayStr);
    }

    /**
     * Parse comma-separated doubles.
     */
    private static double[] parseDoubleArray(String arrayStr) {
        String[] parts = arrayStr.split(",");
        double[] result = new double[parts.length];

        for (int i = 0; i < parts.length; i++) {
            try {
                result[i] = Double.parseDouble(parts[i].trim());
            } catch (NumberFormatException e) {
                return null;
            }
        }

        return result;
    }

    /**
     * Find matching closing bracket.
     */
    private static int findMatchingBracket(String s, int start, char open, char close) {
        int depth = 0;
        for (int i = start; i < s.length(); i++) {
            char c = s.charAt(i);
            if (c == open) depth++;
            else if (c == close) {
                depth--;
                if (depth == 0) return i;
            }
        }
        return -1;
    }

    /**
     * Create a BufferedReader with automatic gzip detection.
     */
    private static BufferedReader createReader(String filePath) throws IOException {
        InputStream is = new FileInputStream(filePath);

        // Check for gzip magic bytes
        if (filePath.toLowerCase().endsWith(".gz")) {
            is = new GZIPInputStream(is);
        }

        return new BufferedReader(new InputStreamReader(is, StandardCharsets.UTF_8));
    }

    /**
     * Get the embedding dimension from a loaded embeddings map.
     *
     * @param embeddings Map of embeddings
     * @return Embedding dimension, or -1 if map is empty
     */
    public static int getEmbeddingDimension(Map<String, double[]> embeddings) {
        if (embeddings == null || embeddings.isEmpty()) {
            return -1;
        }
        return embeddings.values().iterator().next().length;
    }
}
