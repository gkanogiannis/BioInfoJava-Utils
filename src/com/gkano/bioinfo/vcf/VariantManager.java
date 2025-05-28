package com.gkano.bioinfo.vcf;

import java.util.Collection;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;

public class VariantManager<T> {

    private int numVariants;
    private int numSamples;
    private int ploidy;
    private int maxAlleles;

    private ConcurrentLinkedQueue<T> variantRawCache = null;
    private final int maxSizeOfVariantCache;

    Map<String, int[]> genotypeEncodingCache;

    public VariantManager(int maxSizeOfVariantCache) {
        this.variantRawCache = new ConcurrentLinkedQueue<>();
        this.maxSizeOfVariantCache = maxSizeOfVariantCache;
        this.genotypeEncodingCache = new ConcurrentHashMap<>();
    }

    public Map<String, int[]> getEncodingCache() {
        return this.genotypeEncodingCache;
    }

    public void setNumSamples(int numSamples) {
        this.numSamples = numSamples;
    }

    public final int getNumSamples() {
        return numSamples;
    }

    public void setNumVariants(int numVariants) {
        this.numVariants = numVariants;
    }

    public final int getNumVariants() {
        return numVariants;
    }

    public void setPloidy(int ploidy) {
        this.ploidy = ploidy;
    }

    public int getPloidy() {
        return ploidy;
    }

    public int getMaxAlleles() {
        return maxAlleles;
    }

    public void setMaxAlleles(int maxAlleles) {
        this.maxAlleles = maxAlleles;
    }

    @SuppressWarnings("unused")
    private void clear() {
        if (variantRawCache != null) {
            variantRawCache.clear();
        }
        variantRawCache = new ConcurrentLinkedQueue<>();
    }

    public boolean hasMoreRaw() {
        return !this.variantRawCache.isEmpty();
    }

    public boolean putVariantRaw(T variant) {
        try {
            //long free = Runtime.getRuntime().freeMemory();
            //long max = Runtime.getRuntime().maxMemory();
            //long total = Runtime.getRuntime().totalMemory();
            //if(free < 48*1024*1024) {
            //	return false;
            //}
            if(this.variantRawCache.size() >= maxSizeOfVariantCache) {
                return false;
            }
            this.variantRawCache.offer(variant);
            return true;
        } catch (Exception e) {
            return false;
        }
    }

    public T getNextVariantRaw() {
        return this.variantRawCache.poll();
    }

    public float[][] reduce(Collection<VariantProcessor<String>> variantProcessors) {
        int[][] finalDotProd = new int[numSamples][numSamples];
        int[] finalNorm = new int[numSamples];

        for (VariantProcessor<String> vp : variantProcessors) {
            int[][] partialDot = vp.getLocalDotProd();
            int[] partialNorm = vp.getLocalNorm();
            for (int i = 0; i < numSamples; i++) {
                finalNorm[i] += partialNorm[i];
                for (int j = i; j < numSamples; j++) {
                    finalDotProd[i][j] += partialDot[i][j];
                }
            }
        }

        float[][] cosineDist = new float[numSamples][numSamples];
        for (int i = 0; i < numSamples; i++) {
            float normI = (float) Math.sqrt(finalNorm[i]);
            for (int j = i; j < numSamples; j++) {
                float normJ = (float) Math.sqrt(finalNorm[j]);
                float dot = finalDotProd[i][j];
                float similarity = ((normI > 0 && normJ > 0) ? (dot / (normI * normJ)) : (float)0.0);
                float dist =  (float)1.0 - similarity;
                cosineDist[i][j] = cosineDist[j][i] = dist;
            }
        }
        return cosineDist;
    }
    
}
