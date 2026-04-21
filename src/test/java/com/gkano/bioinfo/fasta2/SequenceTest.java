package com.gkano.bioinfo.fasta2;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

class SequenceTest {

    @Test
    void nameAndShortNameParseHeader() {
        byte[] header = ">contig-123 extra description".getBytes();
        byte[] seq = "ACGT".getBytes();
        Sequence s = new Sequence(1, header, seq);
        assertEquals("contig-123 extra description", s.getName());
        assertEquals("contig-123", s.getShortName());
    }

    @Test
    void charAtAndLength() {
        Sequence s = new Sequence(1, ">x".getBytes(), "ACGT".getBytes());
        assertEquals(4, s.getLength());
        assertEquals('A', s.charAt(0));
        assertEquals('T', s.charAt(3));
    }

    @Test
    void kmerCountsInsertAndRetrieve() {
        Sequence s = new Sequence(1, ">x".getBytes(), "ACGT".getBytes());
        s.insertKmerCount(42L, 3);
        s.insertKmerCount(42L, 1);
        assertEquals(4, s.getCountForKmerCode(42L));
    }

    @Test
    void negativeKmerCodesAreIgnored() {
        Sequence s = new Sequence(1, ">x".getBytes(), "ACGT".getBytes());
        s.insertKmerCount(-1L, 5);
        assertEquals(0, s.getCountForKmerCode(-1L));
    }

    @Test
    void sequenceD2AddNucleotideCountsAccumulates() {
        SequenceD2 d2 = new SequenceD2(new Sequence(1, ">x".getBytes(), "ACGT".getBytes()));
        d2.addNucleotideCounts((byte) 2, (byte) 3, (byte) 5, (byte) 7);
        d2.addNucleotideCounts((byte) 1, (byte) 1, (byte) 1, (byte) 1);
        assertEquals(3, d2.as);
        assertEquals(4, d2.ts);
        assertEquals(6, d2.cs);
        assertEquals(8, d2.gs);
        assertEquals(21L, d2.getTotalATCG());
    }

    @Test
    void totalCountsSumsKmerCounts() {
        SequenceD2 d2 = new SequenceD2(new Sequence(1, ">x".getBytes(), "ACGT".getBytes()));
        d2.insertKmerCount(1L, 2);
        d2.insertKmerCount(2L, 3);
        d2.insertKmerCount(3L, 5);
        assertEquals(10L, d2.getTotalCounts());
    }

    @Test
    void calculateProbsIsFiniteForSimpleSequence() {
        SequenceD2 d2 = new SequenceD2(new Sequence(1, ">x".getBytes(), "ACGTACGT".getBytes()));
        // Simulate what SequenceProcessor does for k=2 on "ACGTACGT".
        int k = 2;
        byte[] seq = "ACGTACGT".getBytes();
        for (int i = 0; i + k <= seq.length; i++) {
            int oldAs = d2.as, oldTs = d2.ts, oldCs = d2.cs, oldGs = d2.gs;
            long code = KmerEncoder.encodeToLong(d2, i, i + k, false, true);
            if (code >= 0) {
                d2.insertKmerCount(code, 1);
                d2.insertKmerProb(code,
                        (short) (d2.as - oldAs),
                        (short) (d2.ts - oldTs),
                        (short) (d2.cs - oldCs),
                        (short) (d2.gs - oldGs));
            }
        }
        double sum = d2.calculateProbs(k);
        assertTrue(Double.isFinite(sum));
        assertTrue(sum > 0.0);
    }
}
