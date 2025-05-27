/*
 *
 * BioInfoJava-Utils 
 *
 * Copyright (C) 2021 Anestis Gkanogiannis <anestis@gkanogiannis.com>
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
package com.gkano.bioinfo.fasta;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.zip.GZIPInputStream;

public class FastaStreamingIterator implements Iterable<List<byte[]>> {
    private final List<String> inputPaths;

    public FastaStreamingIterator(List<String> inputPaths) {
        this.inputPaths = inputPaths;
    }

    public FastaStreamingIterator(String... inputPaths) {
        this(Arrays.asList(inputPaths));
    }

    @Override
    public Iterator<List<byte[]>> iterator() {
        return new Iterator<List<byte[]>>() {

            private BufferedReader reader = null;
            private int currentPathIndex = -1;
            private boolean usingStdin = false;
            private String line = null;

            private void advanceFile() {
                try {
                    if (reader != null && !usingStdin) {
                        reader.close();
                    }
                    currentPathIndex++;
                    if (currentPathIndex >= inputPaths.size()) {
                        reader = null;
                        return;
                    }

                    String path = inputPaths.get(currentPathIndex);
                    InputStream in;
                    if (path.equals("-")) {
                        in = System.in;
                        usingStdin = true;
                    } else {
                        File f = new File(path);
                        in = new FileInputStream(f);
                        if (isGzipped(f)) {
                            in = new GZIPInputStream(in);
                        }
                    }

                    reader = new BufferedReader(new InputStreamReader(in));
                } catch (IOException e) {
                    throw new RuntimeException("Error opening input: " + inputPaths.get(currentPathIndex), e);
                }
            }

            private boolean isGzipped(File file) throws IOException {
                try (InputStream in = new FileInputStream(file)) {
                    byte[] sig = new byte[2];
                    int n = in.read(sig);
                    return n == 2 && sig[0] == (byte) 0x1F && sig[1] == (byte) 0x8B;
                }
            }

            //Return a list of lines of full reads from the fasta/fastq file (possibly multiline)
            private List<byte[]> nextRecord() {
                try {
                    List<byte[]> chunk = new ArrayList<>();
                    while (true) {
                        if (reader == null) {
                            advanceFile();
                            if (reader == null) return null;
                        }

                        if (line == null) {
                            line = reader.readLine();
                        }

                        while (line != null && line.trim().isEmpty()) {
                            line = reader.readLine();
                        }

                        if (line == null) {
                            advanceFile(); // advance to next file
                            continue;
                        }

                        // FASTA or FASTQ header
                        if (line.startsWith(">") || line.startsWith("@")) {
                            chunk.add(line.getBytes());
                            String prefix = line.startsWith("@") ? "@" : ">";
                            line = reader.readLine();

                            while (line != null && !line.startsWith(">") && !line.startsWith("@") && !line.startsWith("+")) {
                                chunk.add(line.getBytes());
                                line = reader.readLine();
                            }

                            // If it's FASTQ, add "+" and the quality lines
                            if (prefix.equals("@") && line != null && line.startsWith("+")) {
                                chunk.add(line.getBytes());
                                line = reader.readLine();
                                while (line != null && !line.startsWith("@") && !line.startsWith(">")) {
                                    chunk.add(line.getBytes());
                                    line = reader.readLine();
                                }
                            }

                            return chunk;
                        } else {
                            line = reader.readLine();
                        }
                    }
                } catch (IOException e) {
                    throw new RuntimeException("Error reading FASTA/FASTQ input", e);
                }
            }

            private List<byte[]> nextChunk = null;

            @Override
            public boolean hasNext() {
                if (nextChunk != null) return true;
                nextChunk = nextRecord();
                return nextChunk != null;
            }

            @Override
            public List<byte[]> next() {
                if (!hasNext()) throw new NoSuchElementException();
                List<byte[]> result = nextChunk;
                nextChunk = null;
                return result;
            }
        };
    }
}
