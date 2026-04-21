package com.gkano.bioinfo.fasta2;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class FastaStreamingIteratorTest {

    private static Path writeFile(Path dir, String name, String content) throws IOException {
        Path p = dir.resolve(name);
        Files.writeString(p, content);
        return p;
    }

    private static Path writeGzipFile(Path dir, String name, String content) throws IOException {
        Path p = dir.resolve(name);
        try (GZIPOutputStream out = new GZIPOutputStream(new FileOutputStream(p.toFile()))) {
            out.write(content.getBytes());
        }
        return p;
    }

    private static List<FastaStreamingIterator.FastaRecord> collect(FastaStreamingIterator it) {
        List<FastaStreamingIterator.FastaRecord> out = new ArrayList<>();
        for (var r : it) {
            out.add(r);
        }
        return out;
    }

    @Test
    void parsesSimpleFasta(@TempDir Path tmp) throws IOException {
        Path f = writeFile(tmp, "a.fa",
                ">seq1\nACGT\n>seq2\nTTTT\n");
        try (var it = new FastaStreamingIterator(List.of(f.toString()))) {
            var recs = collect(it);
            assertEquals(2, recs.size());
            assertEquals(">seq1", recs.get(0).headerLine());
            assertEquals("ACGT", recs.get(0).getFullSequence());
            assertTrue(recs.get(0).isFasta());
            assertFalse(recs.get(0).isFastq());
            assertEquals(">seq2", recs.get(1).headerLine());
            assertEquals("TTTT", recs.get(1).getFullSequence());
        }
    }

    @Test
    void parsesMultiLineSequence(@TempDir Path tmp) throws IOException {
        Path f = writeFile(tmp, "multi.fa",
                ">seq1\nACGT\nACGT\nACGT\n>seq2\nTT\nTT\n");
        try (var it = new FastaStreamingIterator(List.of(f.toString()))) {
            var recs = collect(it);
            assertEquals(2, recs.size());
            assertEquals("ACGTACGTACGT", recs.get(0).getFullSequence());
            assertEquals("TTTT", recs.get(1).getFullSequence());
        }
    }

    @Test
    void parsesFastq(@TempDir Path tmp) throws IOException {
        Path f = writeFile(tmp, "q.fq",
                "@seq1\nACGT\n+\n!!!!\n@seq2\nTTTT\n+\nIIII\n");
        try (var it = new FastaStreamingIterator(List.of(f.toString()))) {
            var recs = collect(it);
            assertEquals(2, recs.size());
            assertTrue(recs.get(0).isFastq());
            assertEquals("ACGT", recs.get(0).getFullSequence());
            assertEquals("!!!!", recs.get(0).getFullQuality());
            assertEquals("TTTT", recs.get(1).getFullSequence());
            assertEquals("IIII", recs.get(1).getFullQuality());
        }
    }

    @Test
    void readsGzippedInput(@TempDir Path tmp) throws IOException {
        Path f = writeGzipFile(tmp, "a.fa.gz",
                ">seq1\nACGT\n>seq2\nGGGG\n");
        try (var it = new FastaStreamingIterator(List.of(f.toString()))) {
            var recs = collect(it);
            assertEquals(2, recs.size());
            assertEquals("ACGT", recs.get(0).getFullSequence());
            assertEquals("GGGG", recs.get(1).getFullSequence());
        }
    }

    @Test
    void concatenatesMultipleFiles(@TempDir Path tmp) throws IOException {
        Path a = writeFile(tmp, "a.fa", ">a1\nACGT\n");
        Path b = writeFile(tmp, "b.fa", ">b1\nTTTT\n>b2\nGGGG\n");
        try (var it = new FastaStreamingIterator(List.of(a.toString(), b.toString()))) {
            var recs = collect(it);
            assertEquals(3, recs.size());
            assertEquals(">a1", recs.get(0).headerLine());
            assertEquals(">b1", recs.get(1).headerLine());
            assertEquals(">b2", recs.get(2).headerLine());
        }
    }

    @Test
    void handlesEmptyLinesBetweenRecords(@TempDir Path tmp) throws IOException {
        Path f = writeFile(tmp, "blank.fa",
                "\n>seq1\nACGT\n\n\n>seq2\nTTTT\n\n");
        try (var it = new FastaStreamingIterator(List.of(f.toString()))) {
            var recs = collect(it);
            assertEquals(2, recs.size());
            assertEquals("ACGT", recs.get(0).getFullSequence());
            assertEquals("TTTT", recs.get(1).getFullSequence());
        }
    }

    @Test
    void emptyFileYieldsNoRecords(@TempDir Path tmp) throws IOException {
        Path f = writeFile(tmp, "empty.fa", "");
        try (var it = new FastaStreamingIterator(List.of(f.toString()))) {
            assertEquals(0, collect(it).size());
        }
    }
}
