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
package com.gkano.bioinfo.tree;

import com.gkano.bioinfo.var.Logger;

import gnu.trove.iterator.TLongDoubleIterator;
import gnu.trove.map.hash.TLongDoubleHashMap;

public class ClusterVectorTrove {
	
	private double norm = 0.0;
	
	private TLongDoubleHashMap kmersLo1;
	private TLongDoubleHashMap kmersLo2;
	private TLongDoubleHashMap kmersHi1;
	private TLongDoubleHashMap kmersHi2;
	private long splitKey;
	
	public ClusterVectorTrove(long initialCapacityLong, long splitKey) {
		int initialCapacityInt = (int)((double)(initialCapacityLong/4L)*1.1);
		if(initialCapacityInt>1000000000) initialCapacityInt = 1000000000;
		Logger.info(this, "\tTrove size/4="+initialCapacityInt);
		
		//
		//initialCapacityInt = 1024;
		//
		
		kmersLo1 = new TLongDoubleHashMap(initialCapacityInt);
		kmersLo1.ensureCapacity(initialCapacityInt);
		kmersLo2 = new TLongDoubleHashMap(initialCapacityInt);
		kmersLo2.ensureCapacity(initialCapacityInt);
		kmersHi1 = new TLongDoubleHashMap(initialCapacityInt);
		kmersHi1.ensureCapacity(initialCapacityInt);
		kmersHi2 = new TLongDoubleHashMap(initialCapacityInt);
		kmersHi2.ensureCapacity(initialCapacityInt);
		this.splitKey = splitKey;
	}
	
	public void clear(){
		if(kmersLo1!=null){
			kmersLo1.clear();
			kmersLo1 = null;
		}
		if(kmersLo2!=null){
			kmersLo2.clear();
			kmersLo2 = null;
		}
		if(kmersHi1!=null){
			kmersHi1.clear();
			kmersHi1 = null;
		}
		if(kmersHi2!=null){
			kmersHi2.clear();
			kmersHi2 = null;
		}
	}
	
	public boolean isEmpty() {
		if((kmersLo1==null || kmersLo1.isEmpty()) && (kmersLo2==null || kmersLo2.isEmpty()) && (kmersHi1==null || kmersHi1.isEmpty()) && (kmersHi2==null || kmersHi2.isEmpty())){
			return true;
		}
		return false;
	}
	
	public void insertKMer(long kmerCode, double count){
		try{	
			if(kmerCode<0L){
				return;
			}
			if(kmerCode>splitKey){
				if(kmerCode>((splitKey*15L)/10L)){
					kmersHi2.adjustOrPutValue(kmerCode, count, count);
				}
				else{
					kmersHi1.adjustOrPutValue(kmerCode, count, count);
				}
			}
			else{
				if(kmerCode>((splitKey*5L)/10L)){
					kmersLo2.adjustOrPutValue(kmerCode, count, count);
				}
				else{
					kmersLo1.adjustOrPutValue(kmerCode, count, count);
				}
			}
		}
			
		catch(Exception e){
			System.out.println("\tTrove full"+"\tLo1="+kmersLo1.size()+"\tLo2="+kmersLo2.size()+"\tHi1="+kmersHi1.size()+"\tHi2="+kmersHi2.size());
			//e.printStackTrace();
			//System.exit(1);
		}
		
	}
	
	public void replaceKMer(long kmerCode, double count){
		if(kmerCode<0L){
			return;
		}
		if(kmerCode>splitKey){
			if(kmerCode>((splitKey*15)/10)){
				kmersHi2.put(kmerCode, count);
			}
			else{
				kmersHi1.put(kmerCode, count);
			}
		}
		else{
			if(kmerCode>((splitKey*5)/10)){
				kmersLo2.put(kmerCode, count);
			}
			else{
				kmersLo1.put(kmerCode, count);
			}
		}
	}
	
	public double getCountForKmerCode(long kmerCode){
		if(kmerCode>splitKey){
			if(kmerCode>((splitKey*15)/10)){
				return kmersHi2.get(kmerCode);
			}
			else{
				return kmersHi1.get(kmerCode);
			}
		}
		else{
			if(kmerCode>((splitKey*5)/10)){
				return kmersLo2.get(kmerCode);
			}
			else{
				return kmersLo1.get(kmerCode);
			}
		}
	}
	
	public long size(){
		return (long)kmersLo1.size()+(long)kmersHi1.size()+(long)kmersLo2.size()+(long)kmersHi2.size();
	}
	
	public void normalize(double norm){
		if(isEmpty() || norm==0.0){
			return;
		}
		for ( TLongDoubleIterator it = iterator(); it.hasNext(); ) {
			it.advance();
			long kmerCode = it.key();
			double kmerCount = it.value();
			replaceKMer(kmerCode, kmerCount / norm);
		}
	}
	
	public double getNorm(){
		double n = norm;
		if(n==0.0){
			n = computeNorm();
			norm = n;
		}
		return n;
	}
	
	private double computeNorm(){
		double norm = 0.0;
		try {
			if(isEmpty()){
				return 0.0;
			}
			for ( TLongDoubleIterator it = iterator(); it.hasNext(); ) {
				it.advance();
				norm += it.value()*it.value();
			}
			return Math.sqrt(norm);
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
		return Math.sqrt(norm);
	}
	
	/*
	public double similarityWithClusterVector(ClusterVectorTrove thatVector, Dictionary dictionary) {
		if(thatVector == null || thatVector.isEmpty() || isEmpty()){
			return Double.NEGATIVE_INFINITY;
		}
		
		double similarity = 0.0;
		double thisNorm = getNorm();
		double thatNorm = thatVector.getNorm();
		
		for ( TLongDoubleIterator thatIterator = thatVector.iterator(); thatIterator.hasNext(); ) {
			thatIterator.advance();
			long thatKmerCode = thatIterator.key();
			double thatKmerCount = thatIterator.value();
			double thisKmerCount = getCountForKmerCode(thatKmerCode);
			if(thisKmerCount != 0.0){
				//Integer globalCount = dictionary.getGlobalCountFor(entry.getKey());
				double thisWeight = thisKmerCount; //weightFunction(thisLocalCount, globalCount);
				double thatWeight = thatKmerCount; //weightFunction(thatLocalCount, globalCount);
				similarity += thisWeight * thatWeight;
			}
		}
		
		if(thisNorm==0) thisNorm=1.0;
		if(thatNorm==0) thatNorm=1.0;
		return similarity / (thisNorm * thatNorm);
	}
	*/
	
	public TLongDoubleIterator iterator() {
		return new MyIterator();
	}
	
	private class MyIterator implements TLongDoubleIterator {

		private final TLongDoubleIterator iteratorLo1;
        private final TLongDoubleIterator iteratorLo2;
        private final TLongDoubleIterator iteratorHi1;
        private final TLongDoubleIterator iteratorHi2;
        private int currentLo1;
        private int currentLo2;
        private int currentHi1;
        private int currentHi2;
        private boolean isOnLo1;
        private boolean isOnLo2;
        private boolean isOnHi1;
        private boolean isOnHi2;
        private boolean hasNext;
        
        MyIterator() {
            iteratorLo1 = kmersLo1.iterator();
            iteratorLo2 = kmersLo2.iterator();
            iteratorHi1 = kmersHi1.iterator();
            iteratorHi2 = kmersHi2.iterator();
            currentLo1 = kmersLo1.size();
            currentLo2 = kmersLo2.size();
            currentHi1 = kmersHi1.size();
            currentHi2 = kmersHi2.size();
            isOnLo1 = true;
            isOnLo2 = false;
            isOnHi1 = false;
            isOnHi2 = false;
        }

		public boolean hasNext() {
			if(isOnLo1){
        		hasNext = (currentLo1 > 0);
        		if(!hasNext){
        			isOnLo1 = false;
        			isOnLo2 = true;
        			return hasNext();
        		}
        		return hasNext;
        	}
        	else if(isOnLo2){
        		hasNext = (currentLo2 > 0);
        		if(!hasNext){
        			isOnLo2 = false;
        			isOnHi1 = true;
        			return hasNext();
        		}
        		return hasNext;
        	}
        	else if(isOnHi1){
        		hasNext = (currentHi1 > 0);
        		if(!hasNext){
        			isOnHi1 = false;
        			isOnHi2 = true;
        			return hasNext();
        		}
        		return hasNext;
        	}
        	else if(isOnHi2){
        		hasNext = (currentHi2 > 0);
        		if(!hasNext){
        			isOnHi2 = false;
        			return false;
        		}
        	}
        	return hasNext;
		}

		public void advance() {
			if(!hasNext){
    			return;
    		}
        	
        	if(isOnLo1){
        		iteratorLo1.advance();
                currentLo1--;
        	}
        	else if(isOnLo2){
        		iteratorLo2.advance();
                currentLo2--;
        	}
        	else if(isOnHi1){
        		iteratorHi1.advance();
                currentHi1--;
        	}
        	else if(isOnHi2){
        		iteratorHi2.advance();
                currentHi2--;
        	}
		}
		
		public long key() {
			if(isOnLo1){
        		return iteratorLo1.key();
        	}
        	else if(isOnLo2){
        		return iteratorLo2.key();
        	}
        	else if(isOnHi1){
        		return iteratorHi1.key();
        	}
        	else if(isOnHi2){
        		return iteratorHi2.key();
        	}
			return 0;
		}
		
		public double value() {
			if(isOnLo1){
        		return iteratorLo1.value();
        	}
        	else if(isOnLo2){
        		return iteratorLo2.value();
        	}
        	else if(isOnHi1){
        		return iteratorHi1.value();
        	}
        	else if(isOnHi2){
        		return iteratorHi2.value();
        	}
			return 0.0;
		}

		
		//Unused
		public void remove() {
		}
		public double setValue(double arg0) {
			return 0.0;
		}	
		
	}
	
	/*
	private class MyIterator2 {

        private final TLongDoubleIterator iteratorLo1;
        private final TLongDoubleIterator iteratorLo2;
        private final TLongDoubleIterator iteratorHi1;
        private final TLongDoubleIterator iteratorHi2;
        private int currentLo1;
        private int currentLo2;
        private int currentHi1;
        private int currentHi2;
        private boolean hasNext;
        private boolean isOnLo1;
        private boolean isOnLo2;
        private boolean isOnHi1;
        private boolean isOnHi2;

        MyIterator2() {
            iteratorLo1 = kmersLo1.iterator();
            iteratorLo2 = kmersLo2.iterator();
            iteratorHi1 = kmersHi1.iterator();
            iteratorHi2 = kmersHi2.iterator();
            currentLo1 = kmersLo1.size();
            currentLo2 = kmersLo2.size();
            currentHi1 = kmersHi1.size();
            currentHi2 = kmersHi2.size();
            isOnLo1 = true;
            isOnLo2 = false;
            isOnHi1 = false;
            isOnHi2 = false;
        }

        public boolean hasNext() {
        	if(isOnLo1){
        		hasNext = (currentLo1 > 0);
        		if(!hasNext){
        			isOnLo1 = false;
        			isOnLo2 = true;
        			return hasNext();
        		}
        		return hasNext;
        	}
        	else if(isOnLo2){
        		hasNext = (currentLo2 > 0);
        		if(!hasNext){
        			isOnLo2 = false;
        			isOnHi1 = true;
        			return hasNext();
        		}
        		return hasNext;
        	}
        	else if(isOnHi1){
        		hasNext = (currentHi1 > 0);
        		if(!hasNext){
        			isOnHi1 = false;
        			isOnHi2 = true;
        			return hasNext();
        		}
        		return hasNext;
        	}
        	else if(isOnHi2){
        		hasNext = (currentHi2 > 0);
        		if(!hasNext){
        			isOnHi2 = false;
        			return false;
        		}
        	}
        	return hasNext;
        }

        @Override
        public Pair<Long, Double> next() {
        	if(!hasNext){
    			return null;
    		}
        	
        	if(isOnLo1){
        		iteratorLo1.advance();
        		Pair<Long, Double> pair = new Pair<Long, Double>(iteratorLo1.key(), iteratorLo1.value());
                currentLo1--;
                return pair;
        	}
        	else if(isOnLo2){
        		iteratorLo2.advance();
        		Pair<Long, Double> pair = new Pair<Long, Double>(iteratorLo2.key(), iteratorLo2.value());
                currentLo2--;
                return pair;
        	}
        	else if(isOnHi1){
        		iteratorHi1.advance();
        		Pair<Long, Double> pair = new Pair<Long, Double>(iteratorHi1.key(), iteratorHi1.value());
                currentHi1--;
                return pair;
        	}
        	else if(isOnHi2){
        		iteratorHi2.advance();
        		Pair<Long, Double> pair = new Pair<Long, Double>(iteratorHi2.key(), iteratorHi2.value());
                currentHi2--;
                return pair;
        	}
        	
        	return null;
        }

        @Override
        public void remove() {
        }
    }
    */
}
