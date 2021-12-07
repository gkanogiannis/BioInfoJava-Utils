/*
 *
 * BioInfoJava-Utils ciat.agrobio.core.VariantProcessor
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
package ciat.agrobio.core;

import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;

import ciat.agrobio.io.VCFManager;

public class VariantProcessor implements Runnable {

	private static AtomicInteger variantCount = new AtomicInteger(0);
	
	private static AtomicInteger taskCount = new AtomicInteger(0);
	
	private final int id = taskCount.getAndIncrement();
	private VCFManager vcfm = null;
	private CountDownLatch startSignal = null;
	private CountDownLatch doneSignal = null;
	
	private ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>> individualsToVariantsData;
	//private ConcurrentSkipListSet<String> variantNames;
	
	public static void resetCounters(){
		variantCount = new AtomicInteger(0);
		taskCount = new AtomicInteger(0);
	}
	
	public VariantProcessor(ConcurrentHashMap<String, ConcurrentHashMap<Integer, Object>> individualsToVariantsData, VCFManager vcfm, CountDownLatch startSignal, CountDownLatch doneSignal) {
		this.individualsToVariantsData = individualsToVariantsData;
		this.vcfm = vcfm;
		this.startSignal = startSignal;
		this.doneSignal = doneSignal;
	}
	
	public static AtomicInteger getVariantCount() {
		return variantCount;
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

	public void run() {
		try{
			startSignal.await();
			boolean done = false;
			while(!done){
				Variant variant = vcfm.getNextVariant();
				if(variant==null){
					if(!vcfm.hasMore()){
						done = true;
						break;
					}
					continue;
				}
				//Process variant data
				int count = variantCount.incrementAndGet(); 
				if(count % 10000 == 0){
					System.err.println(GeneralTools.time()+" VariantProcessor: "+id+"\t"+count);
				}
				List<String> variantData = Arrays.asList(variant.getData());
				int variantId = variant.getVariantId();
				//variantNames.add(variantName);
				String format = variantData.get(8);
				int indexGT = Arrays.asList(format.split(":")).indexOf("GT");
				for(int i=9; i< vcfm.getHeaderData().size(); i++) {
					String GT = variantData.get(i).split(":")[indexGT];
					Object GTCode;
					GTCode = GenotypeEncoder.encodeGTComplex(GT);
					individualsToVariantsData.get(vcfm.getHeaderData().get(i)).put(variantId, GTCode);
				}
			}
			doneSignal.countDown();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	
}
