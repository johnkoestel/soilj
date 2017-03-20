package SoilJ.tools;

import SoilJ.tools.ObjectDetector.ColCoords3D;

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

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.PolygonRoi;
import ij.plugin.PlugIn;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;

/** 
 * HistogramStuff is a SoilJ class contains all subroutines that deal with histogram information.
 * 
 * @author John Koestel
 *
 */

public class HistogramStuff implements PlugIn {

	
	public void run(String arg) {
				//ok, this is not needed..
	}

	public class IlluminationInfo {
		
		public int min;
		public int quantile01;
		public int quantile05;
		public int quantile10;
		public int lowerQuartile;
		public int median;
		public int upperQuartile;
		public int quantile90;
		public int quantile95;
		public int quantile99;
		public int max;
		
		public int specialQuantile;
		
		public int mean;
		
	}
	
	public int findQuantileFromHistogram(int[] myHist, double quantile) {
		
		double[] cumHist = calcCumulativeHistogram(myHist);
		int myQuantile = findPercentileFromCumHist(cumHist, quantile);
		
		return myQuantile;
	}

	public int findPercentileFromCumHist(double[] cumHist, double threshold) {		
		
		int myThresh = 0;		
		double maxCumHist = cumHist[cumHist.length - 1]; 
		
		for (int i = 1 ; i < cumHist.length ; i++) {
			if (cumHist[i] / maxCumHist >= threshold) {
				myThresh = i;
				break;
			}
		}
		
		return myThresh;
		
	}
	
	public int findModeFromHistogram(int[] myHist) {
		
		int max = 0;
		int maxVal = 0;
		
		for (int i = 1 ; i < myHist.length ; i++) {
			if (myHist[i] > maxVal) {
				maxVal = myHist[i];
				max = i;
			}
		}		
		
		return max;
		
	}

	public int findPercentileFromTiff(ImageStack nowStack, double threshold) {
		
		ImagePlus nowTiff = new ImagePlus();
		nowTiff.setStack(nowStack);
		
		int[] myHist = sampleHistogram(nowTiff);
		double[] cumHist = calcCumulativeHistogram(myHist);		
		int myP = findPercentileFromCumHist(cumHist,threshold);
		
		return myP;
	}

	public int findNonZeroMinimumFromHistogram(int[] myHist, int largerThan) {
		
		int myMinimum = 0;
		int cc = 1;
		
		while (myMinimum == 0) {
			if (myHist[cc] > largerThan) {	
				myMinimum = cc;				
			}
			cc++;
		}
		
		return myMinimum;
	}

	public int findNonZeroModeFromHistogram(int[] myHist) {
		
		int myModeNumber = 0;
		int myMode = 0;
		for (int i = 1 ; i < myHist.length ; i++) if (myHist[i] > myModeNumber) {
			myModeNumber = myHist[i];
			myMode = i;
		}
		
		return myMode;
	}

	public int findMinFromHistogram(int[] myHist) {
		
		int myMinimum = 0;
		for (int i = 0 ; i < myHist.length ; i++) if (myHist[i] > 0) {
			myMinimum = i;
			break;
		}
		
		return myMinimum;
	}

	public float findMeanFromHistogram(int[] myHist) {
		
		float myMean = 0;
		float myWeightedSum = 0;
		float myNumberOfEntries = 0;
		for (int i = 0 ; i < myHist.length ; i++) if (myHist[i] > 0) {
			myNumberOfEntries += myHist[i];
			myWeightedSum += i * myHist[i];
		}
		
		myMean = myWeightedSum / myNumberOfEntries;
		
		return myMean;
	}

	public double findStdFromHistogram(int[] myHist, float myMean) {
				
		float myWeightedSDSum = 0;
		float myNumberOfEntries = 0;
		for (int i = 0 ; i < myHist.length ; i++) if (myHist[i] > 0) {
			myNumberOfEntries += myHist[i];
			myWeightedSDSum += (i - myMean) * (i - myMean) * myHist[i];
		}
		
		double myStd = Math.sqrt(myWeightedSDSum / (myNumberOfEntries - 1));
		
		return myStd;
	}
	
	public int findMaxFromHistogram(int[] myHist) {
				
		int myMaximum = 0;
		for (int i = 0 ; i < myHist.length ; i++) if (myHist[i] > 0) myMaximum = i;
		
		return myMaximum;
	}

	public int findMedianFromHistogram(int[] myHist) {
		
		double[] cumHist = calcCumulativeHistogram(myHist);
		int myMedian = findPercentileFromCumHist(cumHist, 0.5);
		
		return myMedian;
	}

	public double[] calcCumulativeHistogram(int[] myHist) {
		
		double[] cumHist =  new double[myHist.length];
	
		cumHist[0]=0;
		for (int i = 1 ; i < myHist.length ; i++) {
			cumHist[i] = cumHist[i-1] + myHist[i];
		}	
		
		return cumHist;
	}

	public int[] sampleHistogram(ImagePlus nowTiff) {
	
		int numberOfGreyValues = 0;
		switch (nowTiff.getBitDepth()){
			case 8 : numberOfGreyValues = 256;
			case 16 : numberOfGreyValues = 65536;
			case 32 : {
				numberOfGreyValues = 65536;
				ImageConverter.setDoScaling(true);
				ImageConverter sIC = new ImageConverter(nowTiff);		
				sIC.convertToGray16();
			}
		}	
		
		//nowTiff.show();
		
		int[] myHist = new int[numberOfGreyValues];
		ImageProcessor outIP = null;
		
		//sample histogram
		for (int i = 1 ; i < nowTiff.getStackSize() + 1 ; i++) {  
			nowTiff.setPosition(i);			
			if (nowTiff.getBitDepth()==32) outIP = nowTiff.getProcessor();
			else outIP = nowTiff.getProcessor();
			int[] newHist=outIP.getHistogram();	
			for (int j = 1 ; j < myHist.length - 1; j++) {
				myHist[j] = myHist[j] + newHist[j];
			}			
		}
		
		return myHist;
	}
	
	public double[] convert16to8BitHistogram(double[] my16Hist) {
		
		double[] my8Hist = new double[256];
		
		int cc = 0;
		for (int i = 0 ; i < my16Hist.length ; i += 256) {			
			
			double mysum = 0;
			
			for (int j = 0 ; j < 256 ; j++) {
				mysum += my16Hist[i + j];
			}	
			
			my8Hist[cc] = mysum / 256; 
			cc++;
		}
		
		return my8Hist;
		
	}

	public int[] sampleThisHistogram16(ImagePlus nowTiff, String myName, InputOutput.MyFolderCollection mFC, MenuWaiter.HistogramMenuReturn hMR, int[] myGandS) {
	
		InputOutput jIO = new InputOutput();	
		RoiHandler rH = new RoiHandler();
							
		int[] myHist = new int[256 * 256];
		int[] newHist;	
		int i, j;
				
		//read gauge file if desired
		PolygonRoi[] pRoi = new PolygonRoi[nowTiff.getNSlices()];		
		if (hMR.useInnerCircle) {
			
			//read inner circle
			ObjectDetector jOD = new ObjectDetector();
			ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();
			int versio = jIO.checkInnerCircleFileVersion(mFC.myInnerCircleFiles[myGandS[0]]);			
			if (versio == 0) jCO = jIO.readInnerCircleVer0(mFC.myInnerCircleFiles[myGandS[0]]);	
			else jCO = jIO.readInnerCircleVer1(mFC.myInnerCircleFiles[myGandS[0]]);
			
			pRoi = rH.makeMeAPolygonRoiStack("inner", "exact", jCO, 0);
		}
		
		//get the stacked histogram		
		for (i = 1 ; i < nowTiff.getNSlices() + 1 ; i++) {
			
			IJ.showStatus("Getting 16-bit histogram of slice #" + i + "/" + (nowTiff.getNSlices()));
			
			nowTiff.setPosition(i);		
			
			ImageProcessor myIP = nowTiff.getProcessor();
			ImageProcessor modIP = myIP.duplicate();
			
			//cut out everything outside column
			if (hMR.useInnerCircle) {
				modIP.setRoi(pRoi[i - 1]);
				modIP.setColor(0);
				modIP.fillOutside(pRoi[i - 1]);				
			}
			
			newHist=modIP.getHistogram();
			newHist[0] = 0; //set zero GV to zero
			for (j = 0 ; j < myHist.length ; j++) {
				myHist[j] = myHist[j] + newHist[j];
			}			
			
		}		
		
		return myHist;
		
	}

	public int findTheKnee(int[] myHist) {
		
		int myKnee = 0;
		
		//find maximum and last entry
		double maxVal = 0;
		int myMax = 0;		
		int myLast = 0;
		for (int i = 0 ; i < myHist.length ; i++) {
			if (myHist[i] > maxVal) {
				maxVal = myHist[i];
				myMax = i;
			}
			if (myHist[i] > 0) myLast = i;
		}
		
		//find the coefficients of the spanned string
		double ssA = -maxVal / (myLast - myMax);
		double ssAlpha = - Math.atan(ssA);
		double ssB = Math.tan(ssAlpha) * myLast;
	
		//prepare plot
		double[] xAxis = new double[myHist.length];
		double[] doubleHist = new double[myHist.length];
		double[] spannedString = new double[myHist.length];
		for (int x = 0 ; x < myLast ; x++) {
			xAxis[x] = x;
			doubleHist[x] = myHist[x];
			spannedString[x] = ssA * x + ssB; 
		}
		
		//find line perpendicular to spanned string
		double pBeta = Math.PI / 2 - ssAlpha;
		double pA = Math.tan(pBeta);
		
		//find the longest distance between myHist and spanned string
		double theLongest = 0;		
		for (int x = myMax + 1; x < myLast ; x++){
			double s0 = myHist[x];
			double xNew = (ssB + pA * x - s0) / (pA - ssA);
			double x1 = xNew - x;
			double y = ssA * xNew + ssB;
			double y1 = y - s0; 
			double l =  Math.sqrt(x1*x1 + y1*y1);
			if (l > theLongest) {
				theLongest = l;
				myKnee = x;
			}
			
			//plot it for verification
			double[] lx = new double[(int)Math.round(x1)];
			double[] ly = new double[(int)Math.round(x1)];
			for (int xx = x ; xx < (int)Math.round(xNew) ; xx++) {
				lx[xx - x] = xx;
				ly[xx - x] = s0 + pA * (xx - x); 
			}
			Plot myPlot = new Plot("verification", "greyscale","frequency", xAxis, doubleHist, 4);			
			myPlot.draw();
			myPlot.addPoints(xAxis, spannedString,2);
			myPlot.draw();
			myPlot.addPoints(lx, ly,1);
			//myPlot.setLimits(x - 5, x1 + 5, 0, maxVal);
			myPlot.draw();			
			PlotWindow nowWindow = myPlot.show();
			nowWindow.close();
	
		}
		
		return myKnee;
	}
	
}