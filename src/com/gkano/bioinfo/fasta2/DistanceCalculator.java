/*
 *
 * BioInfoJava-Utils 
 *
 * Copyright (C) 2024 Anestis Gkanogiannis <anestis@gkanogiannis.com>
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
package com.gkano.bioinfo.fasta2;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;
import java.util.concurrent.atomic.AtomicInteger;

import com.gkano.bioinfo.var.Logger;

/**
 * Parallel D2 distance calculator using ForkJoinPool.
 * 
 * Computes pairwise D2_S distances for sequence comparison.
 * Uses work-stealing parallelism for efficient multi-core utilization.
 */
public class DistanceCalculator {
    
    private static final double EPSILON = 1E-15;
    
    private int numThreads = Runtime.getRuntime().availableProcessors();
    private boolean verbose = false;
    
    public DistanceCalculator() {
    }
    
    public DistanceCalculator(int numThreads) {
        this.numThreads = Math.max(1, numThreads);
    }
    
    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
        Logger.setVerbose(verbose);
    }
    
    /**
     * Compute pairwise D2_S distances for all sequences, with rows/cols ordered
     * by the caller-supplied id list. Callers that need matrix rows to line up
     * with a specific label order (e.g. parse order from FastaManager) MUST
     * pass their ordered id list here.
     *
     * @param sequences Map of sequence ID to SequenceD2
     * @param orderedIds Sequence ids in the order rows/cols of the returned
     *                   matrix should follow
     * @return Symmetric distance matrix indexed by orderedIds
     */
    public double[][] computeD2Distances(ConcurrentHashMap<Integer, SequenceD2> sequences,
                                         List<Integer> orderedIds) {
        if (sequences.isEmpty() || orderedIds.isEmpty()) {
            return new double[0][0];
        }

        var seqIds = new ArrayList<>(orderedIds);
        int n = seqIds.size();
        double[][] distances = new double[n][n];

        ForkJoinPool pool = new ForkJoinPool(numThreads);
        try {
            var task = new DistanceComputeTask(sequences, seqIds, distances, verbose);
            pool.execute(task);
            task.join();
        } finally {
            pool.shutdown();
        }

        return distances;
    }

    /**
     * Compute pairwise D2_S distances using the map's own key iteration order.
     * Row/col order is not stable across runs — prefer
     * {@link #computeD2Distances(ConcurrentHashMap, List)} when label alignment
     * matters.
     */
    public double[][] computeD2Distances(ConcurrentHashMap<Integer, SequenceD2> sequences) {
        return computeD2Distances(sequences, new ArrayList<>(sequences.keySet()));
    }
    
    /**
     * Compute distance between two sequences using D2_S metric.
     */
    public static double computeD2S(SequenceD2 x, SequenceD2 y) {
        double d2s = 0.0;
        double normX = 0.0;
        double normY = 0.0;
        
        // Collect all unique k-mers
        var allKmers = new HashSet<Long>();
        
        for (var it = x.getKmerIterator(); it.hasNext();) {
            it.advance();
            allKmers.add(it.key());
        }
        
        for (var it = y.getKmerIterator(); it.hasNext();) {
            it.advance();
            allKmers.add(it.key());
        }
        
        // Compute D2_S statistic
        for (long kmerCode : allKmers) {
            double cxi = x.getDoubleCountForKmerCode(kmerCode);
            double cyi = y.getDoubleCountForKmerCode(kmerCode);
            
            // Expected values
            double pxi = x.getDoubleProbForKmerCode(kmerCode);
            double pyi = y.getDoubleProbForKmerCode(kmerCode);
            
            double cxiBar = cxi - x.getTotalCounts() * pxi;
            double cyiBar = cyi - y.getTotalCounts() * pyi;
            
            double denominator = Math.sqrt(cxiBar * cxiBar + cyiBar * cyiBar);
            if (denominator < EPSILON) {
                denominator = 1.0;
            }
            
            d2s += (cxiBar * cyiBar) / denominator;
            normX += (cxiBar * cxiBar) / denominator;
            normY += (cyiBar * cyiBar) / denominator;
        }
        
        normX = Math.sqrt(normX);
        normY = Math.sqrt(normY);
        
        if (normX == 0.0 || normY == 0.0) {
            return 0.0;
        }
        
        double similarity = d2s / (normX * normY);
        return 0.5 * (1.0 - similarity);
    }
    
    /**
     * Compute distance using D2* metric (variance-based).
     */
    public static double computeD2Star(SequenceD2 x, SequenceD2 y) {
        double d2star = 0.0;
        double normX = 0.0;
        double normY = 0.0;
        
        var allKmers = new HashSet<Long>();
        
        for (var it = x.getKmerIterator(); it.hasNext();) {
            it.advance();
            allKmers.add(it.key());
        }
        
        for (var it = y.getKmerIterator(); it.hasNext();) {
            it.advance();
            allKmers.add(it.key());
        }
        
        for (long kmerCode : allKmers) {
            double cxi = x.getDoubleCountForKmerCode(kmerCode);
            double cyi = y.getDoubleCountForKmerCode(kmerCode);
            
            double pxi = x.getDoubleProbForKmerCode(kmerCode);
            double pyi = y.getDoubleProbForKmerCode(kmerCode);
            
            double cxiBar = cxi - x.getTotalCounts() * pxi;
            double cyiBar = cyi - y.getTotalCounts() * pyi;
            
            double denominator = Math.sqrt(x.getTotalCounts() * pxi) * 
                                Math.sqrt(y.getTotalCounts() * pyi);
            
            if (denominator < EPSILON) {
                denominator = 1.0;
            }
            
            d2star += (cxiBar * cyiBar) / denominator;
            normX += (cxiBar * cxiBar) / denominator;
            normY += (cyiBar * cyiBar) / denominator;
        }
        
        normX = Math.sqrt(normX);
        normY = Math.sqrt(normY);
        
        if (normX == 0.0 || normY == 0.0) {
            return 0.0;
        }
        
        double similarity = d2star / (normX * normY);
        return 0.5 * (1.0 - similarity);
    }
    
    // Inner class for parallel computation
    
    @SuppressWarnings("serial")
    private static class DistanceComputeTask extends RecursiveTask<Void> {
        
        private static final int THRESHOLD = 10;
        private static final AtomicInteger progressCounter = new AtomicInteger(0);
        
        private final ConcurrentHashMap<Integer, SequenceD2> sequences;
        private final List<Integer> seqIds;
        private final double[][] distances;
        private final int start;
        private final int end;
        private final boolean verbose;
        
        DistanceComputeTask(ConcurrentHashMap<Integer, SequenceD2> sequences,
                           List<Integer> seqIds, double[][] distances, boolean verbose) {
            this(sequences, seqIds, distances, 0, seqIds.size(), verbose);
        }
        
        DistanceComputeTask(ConcurrentHashMap<Integer, SequenceD2> sequences,
                           List<Integer> seqIds, double[][] distances,
                           int start, int end, boolean verbose) {
            this.sequences = sequences;
            this.seqIds = seqIds;
            this.distances = distances;
            this.start = start;
            this.end = end;
            this.verbose = verbose;
        }
        
        @Override
        protected Void compute() {
            if (end - start <= THRESHOLD) {
                computeDirectly();
            } else {
                int mid = (start + end) / 2;
                var left = new DistanceComputeTask(sequences, seqIds, distances, start, mid, verbose);
                var right = new DistanceComputeTask(sequences, seqIds, distances, mid, end, verbose);
                left.fork();
                right.compute();
                left.join();
            }
            return null;
        }
        
        private void computeDirectly() {
            for (int i = start; i < end; i++) {
                var seqX = sequences.get(seqIds.get(i));
                
                for (int j = i; j < seqIds.size(); j++) {
                    var seqY = sequences.get(seqIds.get(j));
                    
                    double dist = computeD2S(seqX, seqY);
                    
                    // Clean up small floating point errors
                    if (Math.abs(dist) < EPSILON) {
                        dist = 0.0;
                    }
                    
                    distances[i][j] = dist;
                    distances[j][i] = dist;
                }
                
                if (verbose) {
                    int count = progressCounter.incrementAndGet();
                    if (count % 10 == 0) {
                        //Logger.info(this, "Computed distances for {0}/{1} sequences", new Object[]{count, seqIds.size()});
                    }
                }
            }
        }
    }
}
