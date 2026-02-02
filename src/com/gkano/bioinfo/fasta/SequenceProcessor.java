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

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;

public class SequenceProcessor implements Runnable {

	private static AtomicInteger sequenceCount = new AtomicInteger(0);
	
	private static AtomicInteger taskCount = new AtomicInteger(0);
	
	private final int id = taskCount.getAndIncrement();
	private FastaManager frm = null;
	private int k;
	private boolean normalize;
	private CountDownLatch startSignal = null;
	private CountDownLatch doneSignal = null;
	private boolean verbose = false;

	private ConcurrentHashMap<Integer, SequenceD2> seqVectors;
	
	public static void resetCounters(){
		sequenceCount = new AtomicInteger(0);
		taskCount = new AtomicInteger(0);
	}
	
	public SequenceProcessor(ConcurrentHashMap<Integer, SequenceD2> seqVectors, FastaManager frm, int k, boolean normalize, CountDownLatch startSignal, CountDownLatch doneSignal, boolean verbose) {
		this.seqVectors = seqVectors;
		this.frm = frm;
		this.k = k;
		this.normalize = normalize;
		this.startSignal = startSignal;
		this.doneSignal = doneSignal;
		this.verbose = verbose;
	}
	
	public static AtomicInteger getSequenceCount() {
		return sequenceCount;
	}

	public void setStartSignal(CountDownLatch startSignal) {
		this.startSignal = startSignal;
	}

	public void setDoneSignal(CountDownLatch doneSignal) {
		this.doneSignal = doneSignal;
	}

	public int getId() {
		return id;
	}

	@SuppressWarnings("unused")
	public void run() {
		try{
			startSignal.await();
			boolean done = false;
			while(!done){
				Sequence sequence = frm.getNextSequence();
				SequenceD2 sequenceD2 = null;
				if(sequence!=null){
					sequenceD2 = new SequenceD2(sequence);
					sequenceD2.getName(); sequenceD2.getShortName();
				}
				if(sequenceD2==null){
					if(!frm.hasMore()){
						done = true;
						break;
					}
					continue;
				}
				sequenceCount.incrementAndGet();
				
				long kmerCode;
				for(int i=0; i+k<=sequenceD2.getLength(); i++){	
					int oldAs = sequenceD2.as;
					int oldTs = sequenceD2.ts;
					int oldCs = sequenceD2.cs;
					int oldGs = sequenceD2.gs;
					kmerCode = SequenceProcessor.encodeToLong(sequenceD2, i, i+k, false, true);
					if(kmerCode>=0L){
						sequenceD2.insertKmerCount(kmerCode, 1);
						sequenceD2.insertKmerProb(kmerCode, (short)(sequenceD2.as-oldAs), (short)(sequenceD2.ts-oldTs), (short)(sequenceD2.cs-oldCs), (short)(sequenceD2.gs-oldGs));
					}
					//reverse
					oldAs = sequenceD2.as;
					oldTs = sequenceD2.ts;
					oldCs = sequenceD2.cs;
					oldGs = sequenceD2.gs;
					kmerCode = SequenceProcessor.encodeToLong(sequenceD2, i, i+k, true, true);
					if(kmerCode>=0L){
						sequenceD2.insertKmerCount(kmerCode, 1);
						sequenceD2.insertKmerProb(kmerCode, (short)(sequenceD2.as-oldAs), (short)(sequenceD2.ts-oldTs), (short)(sequenceD2.cs-oldCs), (short)(sequenceD2.gs-oldGs));
					}
				}
				//calculate per sequence probs
				double sumprob = sequenceD2.calculateProbs(k);
				if(normalize) {
					sequenceD2.normalizeProbs(sequenceD2.getNorm());
				}
				SequenceD2 seqVector = sequenceD2;
				//ReadD2Centroid seqVector = new ReadD2Centroid(read);
				seqVectors.put(seqVector.getSequenceId(), seqVector);
								
				if(verbose) System.err.println(sequenceCount.get()+"\t"+seqVector.getShortName());
				
				//Clean sequence
				sequenceD2.clearHeadSeq();

				/*
				System.err.print("\tkmer_Count="+seqVector.getTotalCounts());
				System.err.print("\tATCG="+seqVector.getTotalATCG());
				System.err.print("\tAs="+seqVector.getAs());
				System.err.print("\tTs="+seqVector.getTs());
				System.err.print("\tCs="+seqVector.getCs());
				System.err.print("\tGs="+seqVector.getGs());
				System.err.print("\tnorm="+seqVector.getNorm());
				System.err.println("\tsumprob="+sumprob);
				*/
				//System.err.println(GeneralTools.time()+" END of Vector Creation for Sequence="+seqName);
			}
			doneSignal.countDown();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static long encodeToLong(Sequence sequence, int start, int end, boolean reverse){
		if(!reverse){	
			long num = 0L;
			for(int i=start; i<end; i++){
				num = num * 4L;
				switch (sequence.charAt(i)){
					case 'A':
					case 'a':
						num = num + 0L;
					break;
	
					case 'T':
					case 't':
						num = num + 1L;
					break;
	
					case 'C':
					case 'c':
						num = num + 2L;
					break;
	
					case 'G':
					case 'g':
						num = num + 3L;
					break;
	
					default:
						num = -1L;
					break;
				}
				if (num == -1L){
					break;
				}
			}
			return num;
		}
		else{
			long num = 0L;
			for(int i=end-1; i>=start; i--){
				num = num * 4L;
				switch (sequence.charAt(i)){
					case 'A':
					case 'a':
						num = num + 1L;
					break;
	
					case 'T':
					case 't':
						num = num + 0L;
					break;
	
					case 'C':
					case 'c':
						num = num + 3L;
					break;
	
					case 'G':
					case 'g':
						num = num + 2L;
					break;
	
					default:
						num = -1L;
					break;
				}
				if (num == -1L){
					break;
				}
			}
			return num;
		}
	}
	
	public static long encodeToLong(SequenceD2 sequence, int start, int end, boolean reverse, boolean countNucl){
		if(!reverse){	
			long num = 0L;
			short a=0;
			short t=0;
			short c=0;
			short g=0;
			for(int i=start; i<end; i++){
				num = num * 4L;
				switch (sequence.charAt(i)){
					case 'A':
					case 'a':
						num = num + 0L;
						if(countNucl)
							a++;
					break;
	
					case 'T':
					case 't':
						num = num + 1L;
						if(countNucl)
							t++;
					break;
	
					case 'C':
					case 'c':
						num = num + 2L;
						if(countNucl)
							c++;
					break;
	
					case 'G':
					case 'g':
						num = num + 3L;
						if(countNucl)
							g++;
					break;
	
					default:
						num = -1L;
					break;
				}
				if (num == -1L){
					break;
				}
			}
			if (countNucl && num != -1L){
				sequence.as += a;
				sequence.ts += t;
				sequence.cs += c;
				sequence.gs += g;
			}
			return num;
		}
		else{
			long num = 0L;
			short a=0;
			short t=0;
			short c=0;
			short g=0;
			for(int i=end-1; i>=start; i--){
				num = num * 4L;
				switch (sequence.charAt(i)){
					case 'A':
					case 'a':
						num = num + 1L;
						if(countNucl)
							t++;
					break;
	
					case 'T':
					case 't':
						num = num + 0L;
						if(countNucl)
							a++;
					break;
	
					case 'C':
					case 'c':
						num = num + 3L;
						if(countNucl)
							g++;
					break;
	
					case 'G':
					case 'g':
						num = num + 2L;
						if(countNucl)
							c++;
					break;
	
					default:
						num = -1L;
					break;
				}
				if (num == -1L){
					break;
				}
			}
			if (countNucl && num != -1L){
				sequence.as += a;
				sequence.ts += t;
				sequence.cs += c;
				sequence.gs += g;
			}
			return num;
		}
	}
	
	public static String decodeFromLong(long encoded, int length){
		StringBuilder sbf = new StringBuilder();
		StringBuilder sbr = new StringBuilder();
		for(int i=0; i<length; i++){
			int mod = (int)(encoded%4L);
			switch(mod){
				case 0:
					sbf.append('A');
					sbr.append('T');
					encoded -= 0L;
					break;
				case 1:	
					sbf.append('T');
					sbr.append('A');
					encoded -= 1L;
					break;
				case 2:	
					sbf.append('C');
					sbr.append('G');
					encoded -= 2L;
					break;
				case 3:	
					sbf.append('G');
					sbr.append('C');
					encoded -= 3L;
					break;
				default:	
					break;
			}
			encoded /= 4L;
		}
		String f =  sbf.reverse().toString();
		String r = sbr.toString();
		return f+"/"+r;
	}
	
	public static long encodeToLong(String kmerS, int k){
		long num = 0L;
		for(int i=0; i<k; i++){
			num = num * 4L;
			switch (kmerS.charAt(i)){
				case 'A':
				case 'a':
					num = num + 0L;
				break;

				case 'T':
				case 't':
					num = num + 1L;
				break;

				case 'C':
				case 'c':
					num = num + 2L;
				break;

				case 'G':
				case 'g':
					num = num + 3L;
				break;

				default:
					num = -1L;
				break;
			}
			if (num == -1L){
				break;
			}
		}
		return num;
	}
	
	/*
	public static byte[] encode(Read read, int start, int end){
		int size = ((end-start-1)>>2) + 1;
		byte[] encoded = new byte[size];
		int mod = (end-start)%4;
		int d1,d2,d3,d4;
		byte b;
		for(int i=0; i<size; i++){
			if(i < size-1){
				d1 = DNA.toIndex(read.charAt((i<<2) + 0 + start));
				d2 = DNA.toIndex(read.charAt((i<<2) + 1 + start));
				d3 = DNA.toIndex(read.charAt((i<<2) + 2 + start));
				d4 = DNA.toIndex(read.charAt((i<<2) + 3 + start));
				b =  (byte) (d1 + (d2<<2) + (d3<<4) + (d4<<6));
				encoded[i] = b;
				continue;
			}
			switch (mod) {
			case 0:
				d1 = DNA.toIndex(read.charAt((i<<2) + 0 + start));
				d2 = DNA.toIndex(read.charAt((i<<2) + 1 + start));
				d3 = DNA.toIndex(read.charAt((i<<2) + 2 + start));
				d4 = DNA.toIndex(read.charAt((i<<2) + 3 + start));
				b =  (byte) (d1 + (d2<<2) + (d3<<4) + (d4<<6));
				encoded[i] = b;
				break;
			case 1:
				d1 = DNA.toIndex(read.charAt((i<<2) + 0 + start));
				b = (byte)d1;
				encoded[i] = b;
				break;
			case 2:
				d1 = DNA.toIndex(read.charAt((i<<2) + 0 + start));
				d2 = DNA.toIndex(read.charAt((i<<2) + 1 + start));
				b =  (byte) (d1 + (d2<<2));
				encoded[i] = b;
				break;
			case 3:
				d1 = DNA.toIndex(read.charAt((i<<2) + 0 + start));
				d2 = DNA.toIndex(read.charAt((i<<2) + 1 + start));
				d3 = DNA.toIndex(read.charAt((i<<2) + 2 + start));
				b =  (byte) (d1 + (d2<<2) + (d3<<4));
				encoded[i] = b;
				break;
			default:
				break;
			}
		}
		return encoded;
	}
	*/
	
	/*
	public static char[] decode(byte[] encoded, int k) {
		char[] decoded = new char[k];
		int size = encoded.length;
		int mod = k % 4;
		int d1, d2, d3, d4;
		byte b;
		for (int i = 0; i < size; i++) {
			b = encoded[i];
			if (i < size - 1) {
				d1 = (int) (b & (3));
				d2 = (int) ((b & (3 << 2)) >> 2);
				d3 = (int) ((b & (3 << 4)) >> 4);
				d4 = (int) ((b & (3 << 6)) >> 6);
				decoded[(i << 2) + 0] = DNA.toChar(d1);
				decoded[(i << 2) + 1] = DNA.toChar(d2);
				decoded[(i << 2) + 2] = DNA.toChar(d3);
				decoded[(i << 2) + 3] = DNA.toChar(d4);
				continue;
			}
			switch (mod) {
			case 0:
				d1 = (int) (b & (3));
				d2 = (int) ((b & (3 << 2)) >> 2);
				d3 = (int) ((b & (3 << 4)) >> 4);
				d4 = (int) ((b & (3 << 6)) >> 6);
				decoded[(i << 2) + 0] = DNA.toChar(d1);
				decoded[(i << 2) + 1] = DNA.toChar(d2);
				decoded[(i << 2) + 2] = DNA.toChar(d3);
				decoded[(i << 2) + 3] = DNA.toChar(d4);
				break;
			case 1:
				d1 = (int) (b & (3));
				decoded[(i << 2) + 0] = DNA.toChar(d1);
				break;
			case 2:
				d1 = (int) (b & (3));
				d2 = (int) ((b & (3 << 2)) >> 2);
				decoded[(i << 2) + 0] = DNA.toChar(d1);
				decoded[(i << 2) + 1] = DNA.toChar(d2);
				break;
			case 3:
				d1 = (int) (b & (3));
				d2 = (int) ((b & (3 << 2)) >> 2);
				d3 = (int) ((b & (3 << 4)) >> 4);
				decoded[(i << 2) + 0] = DNA.toChar(d1);
				decoded[(i << 2) + 1] = DNA.toChar(d2);
				decoded[(i << 2) + 2] = DNA.toChar(d3);
				break;
			default:
				break;
			}
		}
		return decoded;
	}
	*/	
}
