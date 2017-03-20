package SoilJ.tools;

/**
 *SoilJ.tools is a collection of classes for SoilJ, 
 *a collection of ImageJ plugins for the semi-automatized processing of 3-D X-ray images of soil columns
 *Copyright 2014 2015 2016 2017 John Koestel
 *
 *This program is free software: you can redistribute it and/or modify
 *it under the terms of the GNU General Public License as published by
 *the Free Software Foundation, either version 3 of the License, or
 *(at your option) any later version.
 *
 *This program is distributed in the hope that it will be useful,
 *but WITHOUT ANY WARRANTY; without even the implied warranty of
 *MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *GNU General Public License for more details.
 *
 *You should have received a copy of the GNU General Public License
 *along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.util.Date;

import ij.plugin.PlugIn;

/** 
 * RollerCaster is a SoilJ class with subroutines for casting arrays.
 * 
 * @author John Koestel
 *
 */

public class RollerCaster implements PlugIn {
	
	final double MilliSecondsOfOneDay = 24*60*60*1000;

	
	public void run(String arg) {
				//ok, this is not needed..
	}

	public float[] intVector2Float(int[] myVector){
		
		float[] outVector =  new float[myVector.length];
		
		for (int i = 0 ; i < myVector.length ; i++) {
			outVector[i] = (float)myVector[i];
		}
		
		return outVector;
	}

	public double[] intVector2Double(int[] myVector){
		
		double[] outVector =  new double[myVector.length];
		
		for (int i = 0 ; i < myVector.length ; i++) {
			outVector[i] = (double)myVector[i];
		}
		
		return outVector;
	}

	public double[] castInt2Double(int[] myArray) {
		
		double[] out = new double[myArray.length];
		
		for (int i = 0 ; i < out.length ; i++) out[i] = (double)myArray[i];
		
		return out;
	}

	public float[] castInt2Float(int[] myArray) {
		
		float[] out = new float[myArray.length];
		
		for (int i = 0 ; i < out.length ; i++) out[i] = (float)myArray[i];
		
		return out;
	}

	public float[] castDouble2Float(double[] myArray) {
		
		float[] out = new float[myArray.length];
		
		for (int i = 0 ; i < out.length ; i++) out[i] = (float)myArray[i];
		
		return out;
	}
	
	public int[] castDouble2Int(double[] myArray) {
		
		int[] out = new int[myArray.length];
		
		for (int i = 0 ; i < out.length ; i++) out[i] = (int)Math.round(myArray[i]);
		
		return out;
	}
	
	public double date2Days(Date date){
	    //  convert a date to an integer and back again
	    double currentTime=date.getTime();
	    currentTime = currentTime/MilliSecondsOfOneDay;
	    return currentTime; 
	}

	public double[] reverseDoubleArray(double[] myArray) {
		
		double[] oldArray = myArray;
		for (int i = 0 ; i < myArray.length ; i++) myArray[i] = oldArray[myArray.length - i - 1];
		
		return myArray;
		
	}
	
}