package com.gkano.bioinfo.vcf;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

class WindowPolicyTest {

    @Test
    void bpModeBoundariesAtMultiplesOfSize() {
        WindowPolicy p = new WindowPolicy(WindowPolicy.Mode.BP, 100, 100, 1);
        assertFalse(p.advance("chr1", 50));
        // pos 150 is in [100,200): boundary fires before incorporating it
        assertTrue(p.advance("chr1", 150));
        WindowPolicy.Window w = p.consumeClosedWindow();
        assertNotNull(w);
        assertEquals("chr1", w.chrom);
        assertEquals(0, w.start);
        assertEquals(100, w.end);
        assertEquals(1, w.count);

        // pos 250 is in [200,300): boundary again
        assertTrue(p.advance("chr1", 250));
        w = p.consumeClosedWindow();
        assertEquals(100, w.start);
        assertEquals(200, w.end);
        assertEquals(1, w.count);
    }

    @Test
    void variantsModeFiresEveryNthVariant() {
        WindowPolicy p = new WindowPolicy(WindowPolicy.Mode.VARIANTS, 3, 3, 1);
        assertFalse(p.advance("chr1", 10));
        assertFalse(p.advance("chr1", 20));
        assertFalse(p.advance("chr1", 30));
        // 4th variant -> boundary fires before incorporating it
        assertTrue(p.advance("chr1", 40));
        WindowPolicy.Window w = p.consumeClosedWindow();
        assertNotNull(w);
        assertEquals("chr1", w.chrom);
        assertEquals(3, w.count);
        assertEquals(10, w.start);
        assertEquals(31, w.end); // last pos + 1
    }

    @Test
    void contigChangeForcesBoundaryInBpMode() {
        WindowPolicy p = new WindowPolicy(WindowPolicy.Mode.BP, 10000, 10000, 1);
        assertFalse(p.advance("chr1", 100));
        assertFalse(p.advance("chr1", 200));
        assertTrue(p.advance("chr2", 100));
        WindowPolicy.Window w = p.consumeClosedWindow();
        assertEquals("chr1", w.chrom);
        assertEquals(2, w.count);
    }

    @Test
    void contigChangeForcesBoundaryInVariantsMode() {
        WindowPolicy p = new WindowPolicy(WindowPolicy.Mode.VARIANTS, 100, 100, 1);
        assertFalse(p.advance("chr1", 100));
        assertFalse(p.advance("chr1", 200));
        assertTrue(p.advance("chr2", 50));
        WindowPolicy.Window w = p.consumeClosedWindow();
        assertEquals("chr1", w.chrom);
        assertEquals(2, w.count);
    }

    @Test
    void slidingStepRejected() {
        assertThrows(UnsupportedOperationException.class,
                () -> new WindowPolicy(WindowPolicy.Mode.BP, 100, 50, 1));
    }

    @Test
    void invalidSizeRejected() {
        assertThrows(IllegalArgumentException.class,
                () -> new WindowPolicy(WindowPolicy.Mode.BP, 0, 0, 1));
    }

    @Test
    void minVariantsDoesNotAffectBoundaryDecision() {
        // boundary still fires regardless of min-variants
        WindowPolicy p = new WindowPolicy(WindowPolicy.Mode.VARIANTS, 2, 2, 5);
        assertFalse(p.advance("chr1", 1));
        assertFalse(p.advance("chr1", 2));
        assertTrue(p.advance("chr1", 3));
        // min-variants only affects whether the caller decides to emit
        assertEquals(5, p.getMinVariants());
    }

    @Test
    void finalizeCurrentEmitsLastOpenWindow() {
        WindowPolicy p = new WindowPolicy(WindowPolicy.Mode.VARIANTS, 100, 100, 1);
        p.advance("chr1", 10);
        p.advance("chr1", 20);
        WindowPolicy.Window w = p.finalizeCurrent();
        assertNotNull(w);
        assertEquals(2, w.count);
        assertNull(p.finalizeCurrent());
    }
}
