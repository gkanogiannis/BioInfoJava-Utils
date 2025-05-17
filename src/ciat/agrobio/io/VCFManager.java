/*
 *
 * BioInfoJava-Utils ciat.agrobio.io.VCFManager
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
package ciat.agrobio.io;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import ciat.agrobio.core.GeneralTools;
import ciat.agrobio.core.Logger;
import ciat.agrobio.core.Variant;
import ciat.agrobio.core.VariantManager;

public class VCFManager implements Runnable{
	private List<byte[]> commentData;
	private byte[][] headerData;
	
	private AtomicInteger currVariantId = new AtomicInteger(0);
	
	private VariantManager vm = null;
	private List<String> inputFileNames = null;
	private boolean done = false;
	private CountDownLatch startSignal = null;
	private CountDownLatch doneSignal = null;

	private boolean useMappedBuffer = false;
	private boolean verbose = false;
	
	public VCFManager(VariantManager vm, List<String> inputFileNames, CountDownLatch startSignal, CountDownLatch doneSignal, boolean useMappedBuffer, boolean verbose) {
		this.vm = vm;
		this.inputFileNames = inputFileNames;
		this.startSignal = startSignal;
		this.doneSignal = doneSignal;
		this.useMappedBuffer = useMappedBuffer;
		this.verbose = verbose;
	
		this.commentData = new ArrayList<byte[]>();
	}
	
	public final List<byte[]> getCommentData() {
		return commentData;
	}
	
	public final byte[][] getHeaderData() {
		return headerData;
	}
	
	public boolean isDone() {
		return this.done;
	}
	
	public void run() {
		try{
			Logger.setVerbose(verbose);
			done = false;
			startSignal.countDown();
			
			Logger.info(this, "START READ\tStreaming:"+(!useMappedBuffer));
		    
		    VCFDecoder decoder = new VCFDecoder();
		    //byte[][] line : split at tabs, byte[] is a string between tabs
			if(useMappedBuffer){
				VCFIterator<byte[][]> iterator = new VCFIterator<>(decoder, inputFileNames);
				for (List<byte[][]> chunk : iterator) {
					//System.err.println(GeneralTools.time()+" VCFManager: Chunk with "+chunk.size()+" lines.");
					//for(byte[][] line : chunk) {
					for(int i=0; i<chunk.size(); i++) {
						byte[][] line = chunk.get(i);
						processVariantLine(line);
					}
					chunk.clear();
					chunk = null;
				}
			}
			else {
				VCFStreamingIterator<byte[][]> iterator = new VCFStreamingIterator<>(decoder, verbose, inputFileNames);
				for (byte[][] line : iterator) {
					processVariantLine(line);	
				}				
			}
		    
			Logger.info(this, "END READ");
			vm.setNumVariants(this.currVariantId.get());
			done = true;
			doneSignal.countDown();
		}
		catch(Exception e){
			e.printStackTrace();
			System.exit(0);
		}
	}

	private void processVariantLine(byte[][] line) {
		try {
			if(line[0][0]=='#' && line[0][1]=='#') {
				commentData.add(line[0]);
			}
			else if(line[0][0]=='#') {
				headerData = line;
				vm.setNumSamples(headerData.length - 9);
			}
			else if(line.length>=10) {
				Variant variant = new Variant(currVariantId.incrementAndGet(), line);
				
				if(vm.getPloidy()<=0) {
					String format = variant.getFormat();
					int indexGT = Arrays.asList(format.split(":")).indexOf("GT");
					String GT = new String(variant.getDataRaw()[9]).split(":")[indexGT];
					if(GT.contains("/")) {
						vm.setPloidy((byte)GT.split("/").length);
					}
					else if(GT.contains("|")) {
						vm.setPloidy((byte)GT.split("\\|").length);
					}
					else {
						vm.setPloidy((byte)1);
						//System.err.println("Cannot understand GT field.");
						//System.exit(vm.getPloidy());
					}
				}
					
				while(!vm.putVariantRaw(variant)) {
					//System.err.println(" VCFManager: Error in inserting "+currVariantId.get()+" variant.");
					//currVariantId.decrementAndGet();
					//i--;
					//variant = null;
					TimeUnit.MILLISECONDS.sleep(500);
				}

				line = null;
			}
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}

}
