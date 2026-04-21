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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.zip.GZIPInputStream;

/**
 * Streaming iterator for FASTA/FASTQ files with support for:
 * - Multiple input files
 * - GZIP compression detection
 * - Multiline sequences
 * - stdin input
 */
public class FastaStreamingIterator implements Iterable<com.gkano.bioinfo.fasta2.FastaStreamingIterator.FastaRecord>, AutoCloseable {
    
    private final List<String> inputPaths;
    private BufferedReader currentReader;
    private int currentPathIndex = -1;
    private boolean usingStdin = false;
    
    public FastaStreamingIterator(List<String> inputPaths) {
        this.inputPaths = new ArrayList<>(inputPaths);
    }
    
    public FastaStreamingIterator(String... inputPaths) {
        this(List.of(inputPaths));
    }
    
    @Override
    public Iterator<FastaRecord> iterator() {
        return new FastaIterator(inputPaths);
    }
    
    @Override
    public void close() throws IOException {
        if (currentReader != null && !usingStdin) {
            currentReader.close();
        }
    }
    
    /**
     * Represents a single sequence record (header + sequence + optional quality).
     */
    static final class FastaRecord {
        private final String headerLine;
        private final List<String> sequenceLines;
        private final List<String> qualityLines;
        
        public FastaRecord(String headerLine, List<String> sequenceLines, List<String> qualityLines) {
            this.headerLine = headerLine;
            this.sequenceLines = sequenceLines;
            this.qualityLines = qualityLines;
        }
        
        public String headerLine() {
            return headerLine;
        }
        
        public List<String> sequenceLines() {
            return sequenceLines;
        }
        
        public List<String> qualityLines() {
            return qualityLines;
        }

        /**
         * Concatenate all sequence lines into a single string.
         */
        public String getFullSequence() {
            return String.join("", sequenceLines);
        }
        
        /**
         * Concatenate all quality lines into a single string.
         */
        public String getFullQuality() {
            return qualityLines.isEmpty() ? "" : String.join("", qualityLines);
        }

        public boolean isFastq() {
            return !qualityLines.isEmpty();
        }

        public boolean isFasta() {
            return qualityLines.isEmpty();
        }
    }
    
    // Private inner iterator implementation
    private class FastaIterator implements Iterator<FastaRecord> {
        private final List<String> paths;
        private BufferedReader reader;
        private int pathIndex = -1;
        private String currentLine;
        private FastaRecord nextRecord;
        private boolean closed = false;
        
        FastaIterator(List<String> paths) {
            this.paths = paths;
        }
        
        @Override
        public boolean hasNext() {
            if (closed) return false;
            if (nextRecord != null) return true;
            nextRecord = fetchNextRecord();
            return nextRecord != null;
        }
        
        @Override
        public FastaRecord next() {
            if (!hasNext()) {
                throw new NoSuchElementException("No more records");
            }
            FastaRecord result = nextRecord;
            nextRecord = null;
            return result;
        }
        
        private FastaRecord fetchNextRecord() {
            try {
                while (true) {
                    if (reader == null) {
                        if (!advanceToNextFile()) {
                            closed = true;
                            return null;
                        }
                    }
                    
                    // Skip empty lines
                    while (currentLine != null && currentLine.trim().isEmpty()) {
                        currentLine = reader.readLine();
                    }
                    
                    if (currentLine == null) {
                        closeCurrentReader();
                        reader = null;
                        continue;
                    }
                    
                    // Check for FASTA/FASTQ header
                    if (currentLine.startsWith(">") || currentLine.startsWith("@")) {
                        return parseRecord();
                    }
                    
                    currentLine = reader.readLine();
                }
            } catch (IOException e) {
                throw new RuntimeException("Error reading sequence file", e);
            }
        }
        
        private boolean advanceToNextFile() throws IOException {
            closeCurrentReader();
            pathIndex++;
            
            if (pathIndex >= paths.size()) {
                return false;
            }
            
            String path = paths.get(pathIndex);
            InputStream in = "-".equals(path) 
                ? System.in 
                : openFileStream(path);
            // Track whether we're reading from stdin so we don't close System.in accidentally
            usingStdin = "-".equals(path);
            
            reader = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8));
            // Prime the buffer by reading the first non-null line; this ensures fetchNextRecord
            // sees the first header and doesn't immediately close the reader.
            currentLine = reader.readLine();
            return true;
        }
        
        private InputStream openFileStream(String path) throws IOException {
            File file = new File(path);
            InputStream in = new FileInputStream(file);
            
            if (isGzipped(file)) {
                in = new GZIPInputStream(in);
            }
            
            return in;
        }
        
        private boolean isGzipped(File file) throws IOException {
            try (InputStream in = new FileInputStream(file)) {
                byte[] sig = new byte[2];
                int n = in.read(sig);
                return n == 2 && sig[0] == (byte) 0x1F && sig[1] == (byte) 0x8B;
            }
        }
        
        private FastaRecord parseRecord() throws IOException {
            String headerLine = currentLine;
            boolean isFastq = headerLine.startsWith("@");
            
            List<String> seqLines = new ArrayList<>();
            List<String> qualLines = new ArrayList<>();
            
            currentLine = reader.readLine();
            
            // Read sequence lines
            while (currentLine != null && !currentLine.startsWith(">") && 
                   !currentLine.startsWith("@") && !currentLine.startsWith("+")) {
                if (!currentLine.trim().isEmpty()) {
                    seqLines.add(currentLine);
                }
                currentLine = reader.readLine();
            }
            
            // For FASTQ, read quality separator and quality lines
            if (isFastq && currentLine != null && currentLine.startsWith("+")) {
                currentLine = reader.readLine();
                
                while (currentLine != null && !currentLine.startsWith("@") && 
                       !currentLine.startsWith(">")) {
                    if (!currentLine.trim().isEmpty()) {
                        qualLines.add(currentLine);
                    }
                    currentLine = reader.readLine();
                }
            }
            
            return new FastaRecord(headerLine, seqLines, qualLines);
        }
        
        private void closeCurrentReader() throws IOException {
            if (reader != null && !usingStdin) {
                reader.close();
            }
        }
    }
}
