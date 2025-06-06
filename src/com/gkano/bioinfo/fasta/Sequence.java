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

import java.nio.charset.StandardCharsets;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import com.gkano.bioinfo.var.GeneralTools;

import gnu.trove.iterator.TLongIntIterator;
import gnu.trove.map.hash.TLongIntHashMap;

public class Sequence {

    private final int sequenceId;

    private byte[] header;
    private byte[] seq;
    private byte[] qual;

    private String name;
    private String shortName;

    private int abundanceCluster = -1;

    private TLongIntHashMap kmerCounts;

    private int[] ranks = null;

    public Sequence(int sequenceId, byte[] header, byte[] seq) {
        this.sequenceId = sequenceId;
        this.header = header;
        this.seq = seq;
        kmerCounts = new TLongIntHashMap();
    }

    public Sequence(int sequenceId, byte[] header, byte[] seq, byte[] qual) {
        this.sequenceId = sequenceId;
        this.header = header;
        this.seq = seq;
        this.qual = qual;
        kmerCounts = new TLongIntHashMap();
    }

    public void initRanks(Map<Long, Integer> spaceRanks) {
        if (ranks == null) {
            int spaceSize = spaceRanks.size();
            ranks = new int[spaceSize];

            //put present kmers in a hashmap
            TreeMap<Long, Integer> unsorted = new TreeMap<>();
            for (TLongIntIterator it = iteratorCounts(); it.hasNext();) {
                it.advance();
                long kmerCode = it.key();
                int kmerCount = it.value();
                unsorted.put(kmerCode, kmerCount);
            }
            //System.out.println("unsorted1 "+unsorted.size());
            //System.out.println(unsorted);

            //insert zeros for non present kmers
            for (long kmerCode : spaceRanks.keySet()) {
                if (!unsorted.containsKey(kmerCode)) {
                    unsorted.put(kmerCode, 0);
                }
            }
            //System.out.println("unsorted2 "+unsorted.size());
            //System.out.println(unsorted);

            //sort the map
            Map<Long, Integer> sorted = GeneralTools.sortByValue(unsorted);
            //System.out.println("sorted "+sorted.size());
            //System.out.println(sorted);

            //write ranks
            int current = 0;
            for (Entry<Long, Integer> entry : sorted.entrySet()) {
                ranks[spaceRanks.get(entry.getKey())] = current++;
            }
            //System.out.println("ranks "+ranks.length);
            //System.out.println(Arrays.toString(ranks));
            //System.out.println("space "+spaceRanks.size());
            //System.out.println(spaceRanks.toString());

            unsorted.clear();
            sorted.clear();
        }
    }

    public int[] getRanks() {
        return ranks;
    }

    /*
	public void remove(ArrayList<Long> removed){
		for(long remove : removed){
			kmers.remove(remove);
		}
	}
     */
    public void insertKmerCount(long kmerCode, int count) {
        if (kmerCode < 0L) {
            return;
        }
        kmerCounts.adjustOrPutValue(kmerCode, count, count);
    }

    public int getCountForKmerCode(long kmerCode) {
        return kmerCounts.get(kmerCode);
    }

    public void clearHead() {
        header = null;
    }

    public void clearSeq() {
        seq = null;
        qual = null;
    }

    public void clearHeadSeq() {
        header = null;
        seq = null;
        qual = null;
    }

    public void clearVec() {
        if (kmerCounts != null) {
            kmerCounts.clear();
            kmerCounts = null;
        }
    }

    public void clearFull() {
        header = null;
        seq = null;
        qual = null;
        clearVec();
    }

    public Integer getSequenceId() {
        return sequenceId;
    }

    public byte[] getHeader() {
        return header;
    }

    public byte[] getSeq() {
        return seq;
    }

    public byte[] getQual() {
        return qual;
    }

    public String getName() {
        if (name == null) {
            try {
                name = new String(getHeader()).substring(1).trim().
                        replaceAll("[^A-Za-z0-9_\\.\\+\\-\\=\\|]", " ").
                        replaceAll("\\s+", " ");
            } catch (Exception e) {
                name = null;
            }
        }
        return name;
    }

    public String getShortName() {
        if (shortName == null) {
            try {
                shortName = new String(getHeader()).substring(1).trim().split("\\s+")[0].
                        replaceAll("[^A-Za-z0-9_\\.\\+\\-\\=\\|]", " ").
                        replaceAll("\\s+", " ");
            } catch (Exception e) {
                shortName = null;
            }
        }
        return shortName;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        if (header != null) {
            sb.append(new String(header, StandardCharsets.UTF_8));
        }
        if (seq != null) {
            sb.append("\n").append(new String(seq, StandardCharsets.UTF_8));
        }
        if (qual != null) {
            sb.append("\n+\n").append(new String(qual, StandardCharsets.UTF_8));
        }
        return sb.toString();
    }

    public int getAbundanceCluster() {
        return abundanceCluster;
    }

    public void setAbundanceCluster(int abundanceCluster) {
        this.abundanceCluster = abundanceCluster;
    }

    public char charAt(int index) {
        return (char) (seq[index] & 0xff);
    }

    public int getLength() {
        return seq == null ? 0 : seq.length;
    }

    public void addWith(Sequence other) {
        for (TLongIntIterator it = other.iteratorCounts(); it.hasNext();) {
            it.advance();
            insertKmerCount(it.key(), it.value());
        }
    }

    public TLongIntIterator iteratorCounts() {
        return new MyIterator();
    }

    private class MyIterator implements TLongIntIterator {

        private final TLongIntIterator iterator;
        private int current;
        private boolean hasNext;

        MyIterator() {
            iterator = kmerCounts.iterator();
            current = kmerCounts.size();
        }

        @Override
        public boolean hasNext() {
            hasNext = (current > 0);
            return hasNext;
        }

        @Override
        public void advance() {
            if (!hasNext) {
                return;
            }
            iterator.advance();
            current--;
        }

        @Override
        public long key() {
            return iterator.key();
        }

        @Override
        public int value() {
            return iterator.value();
        }

        @Override
        public void remove() {
        }

        @Override
        public int setValue(int arg0) {
            return 0;
        }

    }

}
