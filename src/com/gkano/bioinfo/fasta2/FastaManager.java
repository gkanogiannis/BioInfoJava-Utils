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
import java.util.Collections;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import com.gkano.bioinfo.var.Logger;

/**
 * Thread-safe manager for streaming FASTA/FASTQ sequences with integrated processing.
 * 
 * Usage:
 *   var manager = new FastaManager.Builder(files)
 *       .withProcessorThreads(numThreads)
 *       .build();
 *   manager.init();
 *   ConcurrentHashMap<Integer, SequenceD2> results = manager.getResults();
 *   manager.awaitCompletion();
 */
public class FastaManager implements Runnable {
    
    private final BlockingQueue<Sequence> sequenceQueue;
    private final ConcurrentHashMap<Integer, SequenceD2> processingResults;
    private final List<Integer> sequenceIds;
    private final List<String> sequenceNames;
    private final List<String> inputFiles;
    private final CountDownLatch startSignal;
    private final CountDownLatch doneSignal;
    private final FastaStreamingIterator iterator;
    
    private final ExecutorService executor;
    private int processorThreads;
    private final int kmerSize;
    private final boolean normalize;
    private final boolean verbose;
    
    private final AtomicInteger nextSequenceId = new AtomicInteger(0);
    private final AtomicInteger parsedCount = new AtomicInteger(0);
    
    private final boolean keepQualities;
    private volatile boolean isDone = false;
    private volatile Exception parseException;
    
    /**
     * Builder for more flexible FastaManager construction.
     */
    public static class Builder {
        private final List<String> files;
        private CountDownLatch startSignal;
        private CountDownLatch doneSignal;
        private int queueCapacity = 1000;
        private boolean keepQualities = false;
        private int processorThreads = 1;
        private int kmerSize = 4;
        private boolean normalize = false;
        private boolean verbose = false;
        
        public Builder(List<String> files) {
            this.files = new ArrayList<>(files);
        }
        
        public Builder withQueueCapacity(int capacity) {
            this.queueCapacity = capacity;
            return this;
        }
        
        public Builder keepQualities(boolean keep) {
            this.keepQualities = keep;
            return this;
        }
        
        public Builder withSignals(CountDownLatch start, CountDownLatch done) {
            this.startSignal = start;
            this.doneSignal = done;
            return this;
        }
        
        public Builder withProcessorThreads(int threads) {
            this.processorThreads = Math.max(1, threads);
            return this;
        }
        
        public Builder withKmerSize(int kmerSize) {
            this.kmerSize = Math.max(1, kmerSize);
            return this;
        }

        public Builder withNormalization(boolean normalize) {
            this.normalize = normalize;
            return this;
        }

        public Builder withVerbose(boolean verbose) {
            this.verbose = verbose;
            return this;
        }
        
        public FastaManager build() {
            if (startSignal == null) startSignal = new CountDownLatch(1);
            if (doneSignal == null) doneSignal = new CountDownLatch(processorThreads + 1);
            return new FastaManager(files, queueCapacity, keepQualities, startSignal, doneSignal,
                                   processorThreads, kmerSize, normalize, verbose);
        }
    }

    public FastaManager(List<String> inputFiles, boolean keepQualities,
                       CountDownLatch startSignal, CountDownLatch doneSignal) {
        this(inputFiles, 1000, keepQualities, startSignal, doneSignal, 1, 4, false, false);
    }

    public FastaManager(List<String> inputFiles, int queueCapacity, boolean keepQualities,
                       CountDownLatch startSignal, CountDownLatch doneSignal) {
        this(inputFiles, queueCapacity, keepQualities, startSignal, doneSignal, 1, 4, false, false);
    }

    public FastaManager(List<String> inputFiles, int queueCapacity, boolean keepQualities,
                       CountDownLatch startSignal, CountDownLatch doneSignal,
                       int processorThreads, int kmerSize, boolean normalize, boolean verbose) {
        this.inputFiles = new ArrayList<>(inputFiles);
        this.keepQualities = keepQualities;
        this.startSignal = startSignal;
        this.doneSignal = doneSignal;
        this.processorThreads = processorThreads;
        this.kmerSize = kmerSize;
        this.normalize = normalize;
        this.verbose = verbose;
        
        this.sequenceQueue = new LinkedBlockingQueue<>(queueCapacity);
        this.processingResults = new ConcurrentHashMap<>();
        this.sequenceIds = Collections.synchronizedList(new ArrayList<>());
        this.sequenceNames = Collections.synchronizedList(new ArrayList<>());
        this.iterator = new FastaStreamingIterator(inputFiles);
        this.executor = Executors.newFixedThreadPool(processorThreads + 1);
    }
    
    /**
     * Check if more sequences are available.
     */
    public boolean hasMore() {
        return !isDone || !sequenceQueue.isEmpty();
    }
    
    /**
     * Get the next sequence with timeout.
     */
    public Sequence getNextSequence() {
        try {
            return sequenceQueue.poll(100, TimeUnit.MILLISECONDS);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            return null;
        }
    }
    
    /**
     * Get number of sequences parsed so far.
     */
    public int getParsedCount() {
        return parsedCount.get();
    }
    
    /**
     * Get unmodifiable view of sequence IDs.
     */
    public List<Integer> getSequenceIds() {
        return Collections.unmodifiableList(sequenceIds);
    }
    
    /**
     * Get unmodifiable view of sequence names.
     */
    public List<String> getSequenceNames() {
        return Collections.unmodifiableList(sequenceNames);
    }
    
    /**
     * Check if parse encountered an error.
     */
    public Exception getParseException() {
        return parseException;
    }
    
    /**
     * Initialize processing: start all processor threads.
     * Should be called after construction and before run().
     */
    public void init() {
        Logger.setVerbose(verbose);
        if (inputFiles == null || inputFiles.isEmpty()) {
            Logger.error(this, "No FASTA input files provided.");
            throw new IllegalArgumentException("No FASTA input files provided.");
        }

        int cpus = Runtime.getRuntime().availableProcessors();
        processorThreads = (cpus < processorThreads ? cpus : processorThreads);
        Logger.info(this, "cpus=" + cpus + "\tusing=" + processorThreads);

        // Start processor threads
        for (int i = 0; i < processorThreads; i++) {
            executor.submit(new SequenceProcessor(processingResults, this, kmerSize,
                                                 normalize, startSignal, doneSignal, false));
        }
        // Start the manager itself
        executor.submit(this);
    }
    
    /**
     * Wait for all processing to complete.
     */
    public void awaitCompletion() throws InterruptedException {
        doneSignal.await();
        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
    }
    
    /**
     * Get the processed results (SequenceD2 objects with k-mers and probabilities).
     */
    public ConcurrentHashMap<Integer, SequenceD2> getResults() {
        return processingResults;
    }
    
    @Override
    public void run() {
        try {
            isDone = false;
            Logger.info(this, "START: Reading FASTA/FASTQ sequences");
            startSignal.countDown();
            
            for (var record : iterator) {
                processRecord(record);
            }
            
            Logger.info(this, "END: Read " + parsedCount.get() + " sequences");
            isDone = true;
            
        } catch (Exception e) {
            parseException = e;
            Logger.error(this, "Parse error: " + e.getMessage());
        } finally {
            doneSignal.countDown();
        }
    }
    
    private void processRecord(FastaStreamingIterator.FastaRecord record) {
        try {
            byte[] header = record.headerLine().getBytes();
            byte[] sequence = record.getFullSequence().getBytes();
            byte[] quality = record.isFastq() ? record.getFullQuality().getBytes() : null;
            
            if (header.length == 0 || sequence.length == 0) {
                return;
            }
            
            int seqId = nextSequenceId.incrementAndGet();
            var sequence_obj = quality != null && keepQualities
                ? new Sequence(seqId, header, sequence, quality)
                : new Sequence(seqId, header, sequence);
            
            // Cache name for quick access
            sequence_obj.getName();
            sequence_obj.getShortName();
            
            try {
                sequenceQueue.put(sequence_obj);
                sequenceIds.add(seqId);
                sequenceNames.add(sequence_obj.getShortName());
                parsedCount.incrementAndGet();
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                nextSequenceId.decrementAndGet();
            }
            
        } catch (Exception e) {
            Logger.error(this, "Error processing record: " + e.getMessage());
        }
    }
    
    /**
     * Clear all cached data.
     */
    public void clear() {
        sequenceQueue.clear();
        sequenceIds.clear();
        sequenceNames.clear();
        nextSequenceId.set(0);
        parsedCount.set(0);
        isDone = false;
    }
}
