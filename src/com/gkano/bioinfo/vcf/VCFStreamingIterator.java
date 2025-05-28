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
package com.gkano.bioinfo.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PushbackInputStream;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.zip.GZIPInputStream;

public class VCFStreamingIterator<T> implements Iterator<T>, Iterable<T> {

    private final List<String> inputPaths;
    private final VCFDecoderInterface<T> decoder;
    private BufferedReader currentReader;
    private int currentPathIndex;
    private String nextLine;
    private T nextDecoded;
    private boolean usingStdin = false;
    private boolean verbose = false;

    public VCFStreamingIterator(VCFDecoderInterface<T> decoder, boolean verbose, List<String> inputPaths) {
        this.inputPaths = inputPaths;
        this.decoder = decoder;
        this.verbose = verbose;
        this.currentPathIndex = -1;
        advanceFile();  // Open first input
        advance();      // Read first valid line
    }

    public VCFStreamingIterator(VCFDecoderInterface<T> decoder, boolean verbose, String... inputPaths) {
        this(decoder, verbose, Arrays.asList(inputPaths));
    }

    private void advanceFile() {
        try {
            if (currentReader != null && !usingStdin) {
                currentReader.close();
            }

            currentPathIndex++;
            if (currentPathIndex < inputPaths.size()) {
                String path = inputPaths.get(currentPathIndex);
                if(verbose) System.err.println("Reading from: " + path);

                InputStream in;
                if (path.equals("-")) {
                    in = System.in;
                    usingStdin = true;
                } else {
                    File file = new File(path);
                    in = new FileInputStream(file);
                }

                // Automatically detect and wrap gzip if needed
                in = detectAndWrap(in);
                currentReader = new BufferedReader(new InputStreamReader(in));

            } else {
                currentReader = null;
            }

        } catch (IOException e) {
            throw new RuntimeException("Error opening input: " + inputPaths.get(currentPathIndex), e);
        }
    }

    private InputStream detectAndWrap(InputStream in) throws IOException {
        PushbackInputStream pb = new PushbackInputStream(in, 2);
        byte[] signature = new byte[2];
        int read = pb.read(signature);
        pb.unread(signature, 0, read);  // push back bytes

        if (read == 2 && signature[0] == (byte) 0x1F && signature[1] == (byte) 0x8B) {
            return new GZIPInputStream(pb);
        } else {
            return pb;
        }
    }

    private void advance() {
        try {
            while (currentReader != null) {
                nextLine = currentReader.readLine();
                if (nextLine != null) {
                    nextDecoded = decoder.decode(nextLine);
                    if (nextDecoded != null) return;
                } else {
                    advanceFile();  // move to next file or end
                }
            }
            nextDecoded = null;
        } catch (IOException e) {
            throw new RuntimeException("Error reading input", e);
        }
    }

    @Override
    public boolean hasNext() {
        return nextDecoded != null;
    }

    @Override
    public T next() {
        if (!hasNext()) throw new NoSuchElementException();
        T result = nextDecoded;
        advance();
        return result;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("remove() not supported.");
    }

    @Override
    public Iterator<T> iterator() {
        return this;
    }

    public void close() {
        try {
            if (currentReader != null && !usingStdin) {
                currentReader.close();
            }
        } catch (IOException e) {
            // suppress
        }
    }
}
