package soilj.tools;

import java.util.Date;

import ij.plugin.PlugIn;


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