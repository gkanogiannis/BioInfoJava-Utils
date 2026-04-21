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

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Level;

import com.gkano.bioinfo.var.Logger;

/**
 * Processor for k-mer extraction from sequences.
 */
public class SequenceProcessor implements Runnable {
    
    private static final AtomicInteger PROCESSOR_COUNT = new AtomicInteger(0);
    
    private final int processorId;
    private final FastaManager fastaManager;
    private final int kmerSize;
    private final boolean normalize;
    private final CountDownLatch startSignal;
    private final CountDownLatch doneSignal;
    private final boolean verbose;
    private final ConcurrentHashMap<Integer, SequenceD2> resultVectors;

    private final AtomicInteger processedCount = new AtomicInteger(0);

    public SequenceProcessor(ConcurrentHashMap<Integer, SequenceD2> resultVectors,
                            FastaManager fastaManager, int kmerSize, boolean normalize,
                            CountDownLatch startSignal, CountDownLatch doneSignal,
                            boolean verbose) {
        this.processorId = PROCESSOR_COUNT.getAndIncrement();
        this.resultVectors = resultVectors;
        this.fastaManager = fastaManager;
        this.kmerSize = kmerSize;
        this.normalize = normalize;
        this.startSignal = startSignal;
        this.doneSignal = doneSignal;
        this.verbose = verbose;

        Logger.setVerbose(verbose);
    }
    
    public static void resetCounters() {
        PROCESSOR_COUNT.set(0);
    }
    
    public int getProcessedCount() {
        return processedCount.get();
    }
    
    public int getProcessorId() {
        return processorId;
    }
    
    @Override
    public void run() {
        try {
            startSignal.await();
            
            while (fastaManager.hasMore()) {
                Sequence sequence = fastaManager.getNextSequence();
                if (sequence == null) {
                    Thread.yield();
                    continue;
                }
                
                processSequence(sequence);
            }
            
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            //Logger.error(this, "Processor {0} interrupted: {1}", new Object[]{processorId, e.getMessage()});
        } finally {
            doneSignal.countDown();
        }
    }
    
    private void processSequence(Sequence rawSeq) {
        try {
            SequenceD2 seqVector = new SequenceD2(rawSeq);

            // Extract k-mers in both forward and reverse complement directions.
            // RC is always counted to match the original fasta-package behaviour.
            for (int i = 0; i + kmerSize <= seqVector.getLength(); i++) {
                int oldAs = seqVector.as;
                int oldTs = seqVector.ts;
                int oldCs = seqVector.cs;
                int oldGs = seqVector.gs;
                // Forward strand
                long kmerCode = KmerEncoder.encodeToLong(seqVector, i, i + kmerSize, false, true);
                if (kmerCode >= 0) {
                    seqVector.insertKmerCount(kmerCode, 1);
                    seqVector.insertKmerProb(kmerCode, (short)(seqVector.as-oldAs), (short)(seqVector.ts-oldTs), (short)(seqVector.cs-oldCs), (short)(seqVector.gs-oldGs));
                }

                oldAs = seqVector.as;
                oldTs = seqVector.ts;
                oldCs = seqVector.cs;
                oldGs = seqVector.gs;
                // Reverse complement strand
                kmerCode = KmerEncoder.encodeToLong(seqVector, i, i + kmerSize, true, true);
                if (kmerCode >= 0) {
                    seqVector.insertKmerCount(kmerCode, 1);
                    seqVector.insertKmerProb(kmerCode, (short)(seqVector.as-oldAs), (short)(seqVector.ts-oldTs), (short)(seqVector.cs-oldCs), (short)(seqVector.gs-oldGs));
                }
            }
            
            //calculate per sequence probs
            // DON'T COMMENT OUT THIS LINE
			double sumprob = seqVector.calculateProbs(kmerSize);

            // Normalize probabilities if requested
            if (normalize) {
                seqVector.normalizeProbs(seqVector.getNorm());
            }
            
            // Store result
            resultVectors.put(seqVector.getSequenceId(), seqVector);
            processedCount.incrementAndGet();
            
            if (verbose) {
                System.err.printf("Processor %d: %d - %s%n", 
                    processorId, processedCount.get(), seqVector.getShortName());
            }
            
            // Clear raw sequence data to save memory
            seqVector.clearHeadSeq();

            /*
				System.err.print("\tkmer_Count="+seqVector.getTotalCounts());
				System.err.print("\tATCG="+seqVector.getTotalATCG());
				System.err.print("\tAs="+seqVector.getAs());
				System.err.print("\tTs="+seqVector.getTs());
				System.err.print("\tCs="+seqVector.getCs());
				System.err.print("\tGs="+seqVector.getGs());
				System.err.print("\tnorm="+seqVector.getNorm());
				System.err.println("\tsumprob="+sumprob);
			*/
            
        } catch (Exception e) {
            //Logger.error(this, "Error processing sequence: {0}", e.getMessage());
        }
    }
}

/**
 * Utility class for k-mer encoding/decoding.
 * Supports DNA alphabet (ATCG).
 */
final class KmerEncoder {
    
    private static final int BITS_PER_BASE = 2;
    private static final int BASES_PER_KMER = 32; // max 64 bits / 2 bits per base
    
    /**
     * Encode a k-mer to a long, optionally counting nucleotides.
     * Returns -1 if non-standard base encountered.
     */
    public static long encodeToLong(Sequence sequence, int start, int end, 
                                    boolean reverseComplement, boolean countNucl) {
        if (reverseComplement) {
            return encodeReverseComplement(sequence, start, end, countNucl);
        } else {
            return encodeForward(sequence, start, end, countNucl);
        }
    }
    
    private static long encodeForward(Sequence sequence, int start, int end, boolean countNucl) {
        long kmer = 0;
        byte a = 0, t = 0, c = 0, g = 0;
        
        for (int i = start; i < end; i++) {
            char base = sequence.charAt(i);
            kmer = kmer * 4L;
            
            switch (Character.toUpperCase(base)) {
                case 'A':
                    kmer += 0;
                    if (countNucl) a++;
                    break;
                case 'T':
                    kmer += 1;
                    if (countNucl) t++;
                    break;
                case 'C':
                    kmer += 2;
                    if (countNucl) c++;
                    break;
                case 'G':
                    kmer += 3;
                    if (countNucl) g++;
                    break;
                default:
                    return -1;
            }
        }
        
        if (countNucl && sequence instanceof SequenceD2) {
            SequenceD2 d2 = (SequenceD2) sequence;
            d2.addNucleotideCounts(a, t, c, g);
        }
        
        return kmer;
    }
    
    private static long encodeReverseComplement(Sequence sequence, int start, int end, boolean countNucl) {
        long kmer = 0;
        byte a = 0, t = 0, c = 0, g = 0;
        
        for (int i = end - 1; i >= start; i--) {
            char base = sequence.charAt(i);
            kmer = kmer * 4L;
            
            switch (Character.toUpperCase(base)) {
                case 'A':
                    kmer += 1;   // complement T
                    if (countNucl) t++;
                    break;
                case 'T':
                    kmer += 0;   // complement A
                    if (countNucl) a++;
                    break;
                case 'C':
                    kmer += 3;   // complement G
                    if (countNucl) g++;
                    break;
                case 'G':
                    kmer += 2;   // complement C
                    if (countNucl) c++;
                    break;
                default:
                    return -1;
            }
        }
        
        if (countNucl && sequence instanceof SequenceD2) {
            SequenceD2 d2 = (SequenceD2) sequence;
            d2.addNucleotideCounts(a, t, c, g);
        }
        
        return kmer;
    }
    
    /**
     * Decode a k-mer long back to forward and reverse complement strings.
     * Format: "forward/reverse"
     */
    public static String decodeFromLong(long kmer, int length) {
        StringBuilder forward = new StringBuilder();
        StringBuilder reverse = new StringBuilder();
        
        for (int i = 0; i < length; i++) {
            int base = (int) (kmer & 3);
            switch (base) {
                case 0:
                    forward.append('A');
                    reverse.append('T');
                    break;
                case 1:
                    forward.append('T');
                    reverse.append('A');
                    break;
                case 2:
                    forward.append('C');
                    reverse.append('G');
                    break;
                case 3:
                    forward.append('G');
                    reverse.append('C');
                    break;
                default:
                    throw new IllegalStateException("Invalid encoded base: " + base);
            }
            kmer >>= 2;
        }
        
        return forward.reverse().toString() + "/" + reverse.toString();
    }
    
    /**
     * Encode a k-mer string directly.
     */
    public static long encodeString(String kmer) {
        long encoded = 0;
        for (int i = 0; i < kmer.length(); i++) {
            encoded = encoded * 4L;
            switch (Character.toUpperCase(kmer.charAt(i))) {
                case 'A':
                    encoded += 0;
                    break;
                case 'T':
                    encoded += 1;
                    break;
                case 'C':
                    encoded += 2;
                    break;
                case 'G':
                    encoded += 3;
                    break;
                default:
                    return -1;
            }
        }
        return encoded;
    }
    
    private KmerEncoder() {
        // Utility class
    }
}
