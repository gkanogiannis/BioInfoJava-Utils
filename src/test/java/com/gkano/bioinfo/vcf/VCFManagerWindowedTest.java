package com.gkano.bioinfo.vcf;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class VCFManagerWindowedTest {

    private static final String HEADER =
            "##fileformat=VCFv4.2\n"
          + "##contig=<ID=chr1>\n"
          + "##contig=<ID=chr2>\n"
          + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n";

    private static final String VCF_SINGLE_CHROM = HEADER
          + "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n"
          + "chr1\t200\t.\tG\tC\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0\n"
          + "chr1\t300\t.\tC\tA\t.\tPASS\t.\tGT\t1/1\t0/0\t0/1\n"
          + "chr1\t400\t.\tT\tG\t.\tPASS\t.\tGT\t0/0\t1/1\t0/1\n"
          + "chr1\t500\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\t0/1\t0/0\n";

    private static final String VCF_TWO_CHROMS = HEADER
          + "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n"
          + "chr1\t200\t.\tG\tC\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0\n"
          + "chr2\t100\t.\tT\tG\t.\tPASS\t.\tGT\t0/0\t1/1\t0/1\n"
          + "chr2\t200\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\t0/1\t0/0\n";

    private static final class CapturedWindow {
        final String chrom;
        final int start;
        final int end;
        final int count;
        final String[] names;
        final double[][] dist;

        CapturedWindow(String chrom, int start, int end, int count,
                       String[] names, double[][] dist) {
            this.chrom = chrom;
            this.start = start;
            this.end = end;
            this.count = count;
            this.names = names;
            // Defensive copy because the manager zeroes its own arrays after the call.
            this.dist = new double[dist.length][];
            for (int i = 0; i < dist.length; i++) {
                this.dist[i] = dist[i].clone();
            }
        }
    }

    private static List<CapturedWindow> runWindowed(String vcfPath, WindowPolicy policy)
            throws Exception {
        VCFManager vcfm = new VCFManager(
                List.of(vcfPath), 2, SNPEncoder.StringToStringParser, false);
        List<CapturedWindow> captured = new ArrayList<>();
        WindowedDistanceWriter writer = (chrom, start, end, n, names, d) ->
                captured.add(new CapturedWindow(chrom, start, end, n, names, d));
        vcfm.setWindowing(policy, writer);
        vcfm.init();
        new Thread(vcfm).start();
        vcfm.awaitFinalization();
        return captured;
    }

    private static double[][] runStandard(String vcfPath) throws Exception {
        VCFManager vcfm = new VCFManager(
                List.of(vcfPath), 2, SNPEncoder.StringToStringParser, false);
        vcfm.init();
        new Thread(vcfm).start();
        vcfm.awaitFinalization();
        return vcfm.reduceDotProdToDistances();
    }

    private static Path writeVcf(Path dir, String name, String content) throws Exception {
        Path p = dir.resolve(name);
        Files.writeString(p, content);
        return p;
    }

    @Test
    void singleWindowMatchesNonWindowedByteForByte(@TempDir Path tmp) throws Exception {
        Path vcf = writeVcf(tmp, "single.vcf", VCF_SINGLE_CHROM);

        double[][] expected = runStandard(vcf.toString());
        // Window large enough that every variant (5) lands in the same window.
        WindowPolicy policy = new WindowPolicy(WindowPolicy.Mode.VARIANTS, 1000, 1000, 1);
        List<CapturedWindow> got = runWindowed(vcf.toString(), policy);
        assertEquals(1, got.size(), "expected exactly one window");
        CapturedWindow w = got.get(0);
        assertEquals(5, w.count);
        for (int i = 0; i < expected.length; i++) {
            for (int j = 0; j < expected.length; j++) {
                assertEquals(expected[i][j], w.dist[i][j], 1e-12,
                        "mismatch at (" + i + "," + j + ")");
            }
        }
    }

    @Test
    void contigChangeProducesSeparateWindows(@TempDir Path tmp) throws Exception {
        Path vcf = writeVcf(tmp, "two.vcf", VCF_TWO_CHROMS);
        // Large window; only contig change can break it.
        WindowPolicy policy = new WindowPolicy(WindowPolicy.Mode.VARIANTS, 1000, 1000, 1);
        List<CapturedWindow> got = runWindowed(vcf.toString(), policy);
        assertEquals(2, got.size(), "contig change must close the window");
        assertEquals("chr1", got.get(0).chrom);
        assertEquals(2, got.get(0).count);
        assertEquals("chr2", got.get(1).chrom);
        assertEquals(2, got.get(1).count);
    }

    @Test
    void bpWindowsTileGenome(@TempDir Path tmp) throws Exception {
        Path vcf = writeVcf(tmp, "single.vcf", VCF_SINGLE_CHROM);
        // 200bp tiles -> [0,200) {100}, [200,400) {200,300}, [400,600) {400,500}
        WindowPolicy policy = new WindowPolicy(WindowPolicy.Mode.BP, 200, 200, 1);
        List<CapturedWindow> got = runWindowed(vcf.toString(), policy);
        assertEquals(3, got.size());
        assertEquals(0, got.get(0).start);
        assertEquals(200, got.get(0).end);
        assertEquals(1, got.get(0).count);
        assertEquals(200, got.get(1).start);
        assertEquals(400, got.get(1).end);
        assertEquals(2, got.get(1).count);
        assertEquals(400, got.get(2).start);
        assertEquals(600, got.get(2).end);
        assertEquals(2, got.get(2).count);
    }

    @Test
    void longTsvWriterEmitsExpectedRows(@TempDir Path tmp) throws Exception {
        Path vcf = writeVcf(tmp, "single.vcf", VCF_SINGLE_CHROM);

        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try (PrintStream ps = new PrintStream(baos)) {
            VCFManager vcfm = new VCFManager(
                    List.of(vcf.toString()), 2, SNPEncoder.StringToStringParser, false);
            // 200bp tiles -> 3 windows; with 3 samples -> 3 i<j pairs each -> 9 data rows + 1 header
            WindowPolicy policy = new WindowPolicy(WindowPolicy.Mode.BP, 200, 200, 1);
            vcfm.setWindowing(policy, WindowedDistanceWriter.longTsv(ps));
            vcfm.init();
            new Thread(vcfm).start();
            vcfm.awaitFinalization();
        }
        String[] lines = baos.toString().split("\n");
        assertEquals(10, lines.length, "expected 1 header + 9 rows; got: " + baos);
        assertEquals("chrom\tstart\tend\tsample_i\tsample_j\tdist", lines[0]);
        for (int i = 1; i < lines.length; i++) {
            assertEquals(6, lines[i].split("\t", -1).length, "row " + i + " column count");
        }
    }

    @Test
    void perWindowTreesEmittedAndWellFormed(@TempDir Path tmp) throws Exception {
        Path vcf = writeVcf(tmp, "single.vcf", VCF_SINGLE_CHROM);

        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try (PrintStream ps = new PrintStream(baos)) {
            VCFManager vcfm = new VCFManager(
                    List.of(vcf.toString()), 2, SNPEncoder.StringToStringParser, false);
            WindowedDistanceWriter writer = (chrom, start, end, n, names, d) -> {
                ps.println("# window chrom=" + chrom + " start=" + start + " end=" + end
                        + " nvariants=" + n + " nsamples=" + names.length);
                try {
                    ps.println((String) new com.gkano.bioinfo.tree.HierarchicalCluster(false)
                            .hclusteringTree(names, d, null)[0]);
                } catch (Exception e) {
                    throw new RuntimeException(e);
                }
            };
            // 200bp tiles -> 3 windows
            vcfm.setWindowing(new WindowPolicy(WindowPolicy.Mode.BP, 200, 200, 1), writer);
            vcfm.init();
            new Thread(vcfm).start();
            vcfm.awaitFinalization();
        }

        String out = baos.toString();
        long headerCount = out.lines().filter(l -> l.startsWith("# window")).count();
        long newickCount = out.lines().filter(l -> l.endsWith(";")).count();
        assertEquals(3, headerCount);
        assertEquals(3, newickCount);
        assertFalse(out.isEmpty());
        assertTrue(out.contains("s1"));
    }
}
