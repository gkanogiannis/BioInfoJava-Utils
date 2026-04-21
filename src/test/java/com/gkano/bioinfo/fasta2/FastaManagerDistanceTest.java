package com.gkano.bioinfo.fasta2;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

/**
 * End-to-end tests exercising the fasta2 pipeline:
 * FastaManager -> SequenceProcessor -> DistanceCalculator.
 */
class FastaManagerDistanceTest {

    private static Path writeFasta(Path dir, String name, String content) throws IOException {
        Path p = dir.resolve(name);
        Files.writeString(p, content);
        return p;
    }

    private static ConcurrentHashMap<Integer, SequenceD2> processFile(String path, int k) throws InterruptedException {
        var manager = new FastaManager.Builder(List.of(path))
                .withProcessorThreads(1)
                .withKmerSize(k)
                .withNormalization(false)
                .build();
        manager.init();
        manager.awaitCompletion();
        return manager.getResults();
    }

    @Test
    void managerProducesOneResultPerSequence(@TempDir Path tmp) throws Exception {
        Path f = writeFasta(tmp, "a.fa",
                ">s1\nACGTACGTACGTACGT\n>s2\nTTTTTTTTTTTTTTTT\n>s3\nGGGGCCCCAAAATTTT\n");
        var results = processFile(f.toString(), 4);
        assertEquals(3, results.size());
        for (SequenceD2 d2 : results.values()) {
            assertTrue(d2.getTotalCounts() > 0,
                    "processed sequence should have non-zero k-mer counts");
        }
    }

    @Test
    void selfDistanceIsZero(@TempDir Path tmp) throws Exception {
        Path f = writeFasta(tmp, "a.fa",
                ">s1\nACGTACGTACGTACGTACGTACGT\n>s2\nAAAATTTTGGGGCCCCAAAATTTT\n");
        var results = processFile(f.toString(), 4);
        double[][] m = new DistanceCalculator(1).computeD2Distances(results);
        assertEquals(2, m.length);
        assertEquals(0.0, m[0][0], 1e-9, "self-distance must be 0");
        assertEquals(0.0, m[1][1], 1e-9, "self-distance must be 0");
    }

    @Test
    void distanceMatrixIsSymmetric(@TempDir Path tmp) throws Exception {
        Path f = writeFasta(tmp, "a.fa",
                ">s1\nACGTACGTACGTACGTACGT\n"
              + ">s2\nAAAATTTTGGGGCCCCAAAA\n"
              + ">s3\nCCGGCCGGTTAACCGGAATT\n");
        var results = processFile(f.toString(), 4);
        double[][] m = new DistanceCalculator(1).computeD2Distances(results);
        int n = m.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                assertEquals(m[i][j], m[j][i], 1e-12,
                        "distance matrix should be symmetric at (" + i + "," + j + ")");
            }
        }
    }

    @Test
    void identicalSequencesHaveZeroDistance(@TempDir Path tmp) throws Exception {
        Path f = writeFasta(tmp, "a.fa",
                ">s1\nACGTACGTACGTACGTACGT\n>s2\nACGTACGTACGTACGTACGT\n");
        var results = processFile(f.toString(), 4);
        double[][] m = new DistanceCalculator(1).computeD2Distances(results);
        assertEquals(0.0, m[0][1], 1e-9);
        assertEquals(0.0, m[1][0], 1e-9);
    }

    @Test
    void distancesAreBoundedInZeroToOne(@TempDir Path tmp) throws Exception {
        Path f = writeFasta(tmp, "a.fa",
                ">s1\nACGTACGTACGTACGTACGTACGT\n"
              + ">s2\nAAAATTTTGGGGCCCCAAAATTTT\n"
              + ">s3\nCCGGCCGGTTAACCGGAATTCCGG\n");
        var results = processFile(f.toString(), 4);
        double[][] m = new DistanceCalculator(1).computeD2Distances(results);
        for (double[] row : m) {
            for (double v : row) {
                assertTrue(v >= -1e-12 && v <= 1.0 + 1e-12,
                        "D2S-derived distance out of expected [0,1] range: " + v);
            }
        }
    }
}
