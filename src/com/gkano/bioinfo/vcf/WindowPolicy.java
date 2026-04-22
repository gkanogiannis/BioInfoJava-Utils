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

public final class WindowPolicy {

    public enum Mode { NONE, BP, VARIANTS }

    public static final class Window {
        public final String chrom;
        public final int start;
        public final int end;
        public final int count;
        public Window(String chrom, int start, int end, int count) {
            this.chrom = chrom;
            this.start = start;
            this.end = end;
            this.count = count;
        }
    }

    private final Mode mode;
    private final int size;
    private final int step;
    private final int minVariants;

    private String curChrom;
    private int curStart;
    private int curEnd;
    private int curCount;
    private Window lastClosed;

    public WindowPolicy(Mode mode, int size, int step, int minVariants) {
        if (mode == null) throw new IllegalArgumentException("mode must not be null");
        if (mode != Mode.NONE) {
            if (size <= 0) throw new IllegalArgumentException("window size must be > 0");
            if (step <= 0) throw new IllegalArgumentException("step must be > 0");
            if (step != size) {
                throw new UnsupportedOperationException(
                        "sliding windows (--step != --window-*) not yet implemented");
            }
            if (minVariants < 1) throw new IllegalArgumentException("--min-variants must be >= 1");
        }
        this.mode = mode;
        this.size = size;
        this.step = step;
        this.minVariants = minVariants;
    }

    public Mode getMode() { return mode; }
    public int getSize() { return size; }
    public int getStep() { return step; }
    public int getMinVariants() { return minVariants; }

    /**
     * Decide if a window boundary fires *before* incorporating this variant.
     * When it does, the just-closed window is held and retrievable via
     * {@link #consumeClosedWindow()}, and a fresh window is opened around (chrom, pos).
     *
     * @return true iff a boundary fired before this variant
     */
    public boolean advance(String chrom, int pos) {
        if (mode == Mode.NONE) {
            curCount++;
            return false;
        }

        if (curChrom == null) {
            openWindow(chrom, pos);
            curCount = 1;
            return false;
        }

        boolean boundary;
        if (!curChrom.equals(chrom)) {
            boundary = true;
        } else if (mode == Mode.BP) {
            boundary = (pos >= curEnd) || (pos < curStart);
        } else { // VARIANTS
            boundary = (curCount >= size);
        }

        if (boundary) {
            lastClosed = new Window(curChrom, curStart, curEnd, curCount);
            openWindow(chrom, pos);
            curCount = 1;
            return true;
        }

        curCount++;
        if (mode == Mode.VARIANTS) {
            curEnd = pos + 1;
        }
        return false;
    }

    private void openWindow(String chrom, int pos) {
        curChrom = chrom;
        if (mode == Mode.BP) {
            curStart = (pos / size) * size;
            curEnd = curStart + size;
        } else if (mode == Mode.VARIANTS) {
            curStart = pos;
            curEnd = pos + 1;
        } else {
            curStart = 0;
            curEnd = 0;
        }
    }

    /**
     * Retrieve (and clear) the window that was just closed by the most recent
     * {@link #advance(String, int)} call that returned true.
     */
    public Window consumeClosedWindow() {
        Window w = lastClosed;
        lastClosed = null;
        return w;
    }

    /**
     * Close and return the currently-open window. Used by the producer at EOF
     * to emit the final window. Returns null if no window is open or it is empty.
     */
    public Window finalizeCurrent() {
        if (curChrom == null || curCount == 0) return null;
        Window w = new Window(curChrom, curStart, curEnd, curCount);
        curChrom = null;
        curCount = 0;
        curStart = 0;
        curEnd = 0;
        return w;
    }
}
