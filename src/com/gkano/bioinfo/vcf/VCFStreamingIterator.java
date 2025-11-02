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

import org.itadaki.bzip2.BZip2InputStream;
import org.tukaani.xz.XZInputStream;

import com.gkano.bioinfo.var.Logger;

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
        Logger.setVerbose(this.verbose);
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
                Logger.info(this, "Reading from: " + path);

                InputStream in;
                if (path.equals("-")) {
                    in = System.in;
                    usingStdin = true;
                } else {
                    File file = new File(path);
                    in = new FileInputStream(file);
                }

                // Automatically detect and wrap compressed stream if needed
                in = detectAndWrap(in);
                currentReader = new BufferedReader(new InputStreamReader(in));

            } else {
                currentReader = null;
            }

        } catch (IOException e) {
            Logger.error(this, "Error opening input: " + inputPaths.get(currentPathIndex));
            Logger.error(this, e.getMessage());
            System.exit(1);
        }
    }

    private InputStream detectAndWrap(InputStream in) throws IOException {
        final int MAX_MAGIC = 6;
        PushbackInputStream pb = new PushbackInputStream(in, MAX_MAGIC);

        byte[] signature = new byte[MAX_MAGIC];
        int read = pb.read(signature, 0, signature.length);

        // If we read anything, push it back so downstream consumers see a full stream
        if (read > 0) {
            pb.unread(signature, 0, read);
        } else if (read == -1) {
            // Empty stream
            Logger.info(this, "Empty stream; no compression detected.");
            return pb;
        }

        // gzip magic: 1F 8B
        if (read >= 2
                && signature[0] == (byte) 0x1F
                && signature[1] == (byte) 0x8B) {
            Logger.info(this, "Detected gzip compression.");
            return new GZIPInputStream(pb);
        }

        // bzip2 magic: 42 5A 68 ("BZh")
        if (read >= 3
                && signature[0] == (byte) 0x42
                && signature[1] == (byte) 0x5A
                && signature[2] == (byte) 0x68) {
            Logger.info(this, "Detected bzip2 compression.");
            // 'true' to decompress concatenated streams if present
            return new BZip2InputStream(pb, false);
        }

        // xz magic: FD 37 7A 58 5A 00
        if (read >= 6
                && signature[0] == (byte) 0xFD
                && signature[1] == (byte) 0x37
                && signature[2] == (byte) 0x7A
                && signature[3] == (byte) 0x58
                && signature[4] == (byte) 0x5A
                && signature[5] == (byte) 0x00) {
            Logger.info(this, "Detected xz compression.");
            return new XZInputStream(pb);
        }

        Logger.info(this, "No compression detected.");
        return pb;
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
            Logger.error(this, "Error reading input");
            System.exit(1);
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
