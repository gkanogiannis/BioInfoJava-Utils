package com.gkano.bioinfo.fasta2;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

class KmerEncoderTest {

    private static Sequence seqOf(String s) {
        return new Sequence(1, (">t\n").getBytes(), s.getBytes());
    }

    @Test
    void encodeStringMatchesManualCalculation() {
        // A=0, T=1, C=2, G=3. For "ATCG": 0, then *4+1=1, *4+2=6, *4+3=27.
        assertEquals(27L, KmerEncoder.encodeString("ATCG"));
        assertEquals(0L, KmerEncoder.encodeString("AAAA"));
        assertEquals(255L, KmerEncoder.encodeString("GGGG")); // 3+12+48+192
    }

    @Test
    void encodeStringIsCaseInsensitive() {
        assertEquals(KmerEncoder.encodeString("ATCG"), KmerEncoder.encodeString("atcg"));
    }

    @Test
    void encodeStringReturnsMinusOneForNonStandardBase() {
        assertEquals(-1L, KmerEncoder.encodeString("ATNG"));
        assertEquals(-1L, KmerEncoder.encodeString("----"));
    }

    @Test
    void encodeForwardMatchesEncodeString() {
        Sequence s = seqOf("ATCGATCG");
        long k = KmerEncoder.encodeToLong(s, 0, 4, false, false);
        assertEquals(KmerEncoder.encodeString("ATCG"), k);
    }

    @Test
    void encodeReverseComplementIsEncodeStringOfRevComp() {
        Sequence s = seqOf("ATCG");
        // Reverse complement of "ATCG" is "CGAT".
        long rc = KmerEncoder.encodeToLong(s, 0, 4, true, false);
        assertEquals(KmerEncoder.encodeString("CGAT"), rc);
    }

    @Test
    void encodeReturnsMinusOneOnNonStandardBase() {
        Sequence s = seqOf("ATNG");
        assertEquals(-1L, KmerEncoder.encodeToLong(s, 0, 4, false, false));
        assertEquals(-1L, KmerEncoder.encodeToLong(s, 0, 4, true, false));
    }

    @Test
    void decodeRoundTripForward() {
        long code = KmerEncoder.encodeString("ATCG");
        String decoded = KmerEncoder.decodeFromLong(code, 4);
        // Format: "forward/reverse"
        String[] parts = decoded.split("/");
        assertEquals(2, parts.length);
        assertEquals("ATCG", parts[0]);
        // Reverse strand uses a complementing map, but bit-for-bit decode only
        // complements within this function; verify both halves are 4 long.
        assertEquals(4, parts[1].length());
    }

    @Test
    void countNuclPopulatesSequenceD2Counts() {
        SequenceD2 d2 = new SequenceD2(seqOf("AATCG"));
        KmerEncoder.encodeToLong(d2, 0, 4, false, true); // AATC
        assertEquals(2, d2.as);
        assertEquals(1, d2.ts);
        assertEquals(1, d2.cs);
        assertEquals(0, d2.gs);
    }

    @Test
    void reverseComplementEncodingDiffersFromForwardForNonPalindrome() {
        // "AAGT" is not its own reverse complement (RC is "ACTT")
        Sequence s = seqOf("AAGT");
        long fwd = KmerEncoder.encodeToLong(s, 0, 4, false, false);
        long rev = KmerEncoder.encodeToLong(s, 0, 4, true, false);
        assertNotEquals(fwd, rev);
        assertTrue(fwd >= 0 && rev >= 0);
        assertEquals(KmerEncoder.encodeString("AAGT"), fwd);
        assertEquals(KmerEncoder.encodeString("ACTT"), rev);
    }

    @Test
    void reverseComplementOfPalindromeEqualsForward() {
        // "ACGT" IS its own reverse complement — sanity-check the invariant.
        Sequence s = seqOf("ACGT");
        long fwd = KmerEncoder.encodeToLong(s, 0, 4, false, false);
        long rev = KmerEncoder.encodeToLong(s, 0, 4, true, false);
        assertEquals(fwd, rev);
    }
}
