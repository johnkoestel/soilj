package soilj.tools;

import java.awt.Color;

import org.apache.commons.math3.stat.StatUtils;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class TailoredMaths implements PlugIn {

	
	public void run(String arg) {
				//ok, this is not needed..
	}

	public double[] oneDMedianFilter(double[] greyAtThisAngle, int footprintOfMedianFilter) {
		
		int j, k;
		
		//apply a median filter to greyAtThisAngle
		double[] mGreyAtThisAngle = new double[greyAtThisAngle.length];
		double[] withinTheFilter = new double[footprintOfMedianFilter];
		for (j = 0 ; j < greyAtThisAngle.length ; j++) mGreyAtThisAngle[j] = 0; //nullning av median filtered radial profile 
		for (j = (footprintOfMedianFilter - 1) / 2 ; j < greyAtThisAngle.length - (footprintOfMedianFilter - 1) / 2; j++) {
			for (k = 0 ; k < footprintOfMedianFilter ; k++) withinTheFilter[k] = greyAtThisAngle[j - (footprintOfMedianFilter - 1) / 2 + k]; 
			mGreyAtThisAngle[j] = StatUtils.percentile(withinTheFilter, 50);
		}
		
		return mGreyAtThisAngle;
	}

	public double[] oneDFirstDerivative(double[] mGreyAtThisAngle, int footprintOfMedianFilter) {
		
		int j;		
		double[] dMG = new double[mGreyAtThisAngle.length];
		for (j = 0 ; j < mGreyAtThisAngle.length ; j++) dMG[j] = 0; //nullning av first derivative
		for (j = (footprintOfMedianFilter - 1) / 2 ; j < mGreyAtThisAngle.length - (footprintOfMedianFilter - 1) / 2 - 1; j++) {
			dMG[j] = mGreyAtThisAngle[j + 1] - mGreyAtThisAngle[j];
		}
		
		return dMG;
	}

	public double[] LinearLOESSFilter(double[] data, int windowHalfSize) {
				
		FitStuff fs =  new FitStuff();
		
		double[] smoothie = new double[data.length];
		
		//create prolonged data series for smoothing
		double[] xData = new double[data.length + 2 * windowHalfSize];
		
		for (int i = 0 ; i < windowHalfSize ; i++) {
			xData[i] = data[0];
			xData[windowHalfSize + data.length + i] = data[data.length - 1];
		}
		
		for (int i = 0 ; i < data.length ; i++) xData[windowHalfSize + i] = data[i];
		
		//go and smooth!
		for (int i = 0 ; i < data.length ; i++) {
					
			double[] nowWindow = new double[2 * windowHalfSize];
			double[] xWindow = new double[2 * windowHalfSize];
			
			for (int j = 0 ; j < 2 * windowHalfSize ; j++) nowWindow[j] = xData[i + j];
			
			for (int j = 0 ; j < nowWindow.length ; j++) xWindow[j] = j;
			
			smoothie[i] = fs.evalLinearFunction(fs.fitLinearFunction(xWindow, nowWindow), windowHalfSize);				
					
		}
		
		return smoothie;
	}

	public double[][] getXYOfEllipseFromAngle(double[] angleAtThisAngle, double xCenter, double yCenter, double major, double minor, double theta) {
		
		int j;
	
		double[][] xy = new double[angleAtThisAngle.length][2];
	
		for (j = 0 ; j < angleAtThisAngle.length ; j++) {			
			double alpha = angleAtThisAngle[j] - theta + Math.PI / 2;
			double a = major;
			double b = minor;			
			xy[j][0] = xCenter - a * Math.cos(alpha) * Math.cos(theta) + b * Math.sin(alpha) * Math.sin(theta);
			xy[j][1] = yCenter + a * Math.cos(alpha) * Math.sin(theta) + b * Math.sin(alpha) * Math.cos(theta);			
		}
	
		return xy;
	}

	public float[] logitCorrectionFunctionCalculator(int i, double[][][] bhc, int bhcEntry, float[] sR) {
		
		float[] ftot = new float[sR.length];
		int standardRadius = sR.length;
		
		double my = bhc[i][0][bhcEntry];
		double sd = bhc[i][1][bhcEntry];
		double thresh = bhc[i][2][bhcEntry];
		double skipslast = bhc[i][3][bhcEntry];
		double r = bhc[i][4][bhcEntry];
		double dy = bhc[i][5][bhcEntry];		
		double[] f = new double[(int)Math.round(standardRadius - skipslast - thresh)];
		
		for (int j = 0 ; j < (int)Math.round(standardRadius - skipslast - thresh); j++) {
			f[(int)Math.round(standardRadius - skipslast - thresh) - j - 1] = dy * (1 - 1 / (1 + Math.exp(-(j - my) / sd))) + r;	 //switch vector around!!!! important!!
		}
		
		int cc=0;
		for (int j = 0 ; j < skipslast ; j++) {
			ftot[cc] = (float)StatUtils.min(f);
			cc++;
		}
		for (int j = 0 ; j < f.length; j++) {
			ftot[cc] = (float)f[j];
			cc++;
		}
		for (int j = 0 ; j < thresh ; j++) {
			ftot[cc] = (float)StatUtils.max(f);
			cc++;
		}
			
		return ftot;
	}

	public double min(double[] myArray) {		
		double myMin = StatUtils.min(myArray);		
		return myMin;
	}
	
	public double minNoZero(double[] myArray) {
		double myMax = StatUtils.max(myArray);
		for (int i = 0 ; i < myArray.length ; i++) if (myArray[i] == 0) myArray[i] = myMax;
		double myMin = StatUtils.min(myArray);		
		return myMin;
	}

	public double max(double[] myArray) {		
		double myMax = StatUtils.max(myArray);		
		return myMax;
	}

	public static int factorial(int n) {
	    int fact = 1; // this  will be the result
	    for (int i = 1; i <= n; i++) {
	        fact *= i;
	    }
	    return fact;
	}
	
}
	