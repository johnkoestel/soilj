package soilj.tools;

import java.awt.Color;
import java.io.File;

import org.apache.commons.math3.stat.StatUtils;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class RoiHandler implements PlugIn {

	
	public void run(String arg) {
				//ok, this is not needed..
	}
	
	public PolygonRoi makeRoiFromFittedEllipse(FitStuff.FittedEllipse fE) {
		
		TailoredMaths math = new TailoredMaths();
		
		int j, cc;
		int maxAlpha = 360;
		int dAlpha = 5;	
		double[] myAngles = new double[maxAlpha/dAlpha];
		float[] x = new float[maxAlpha/dAlpha];
		float[] y = new float[maxAlpha/dAlpha];
		double[][] xy = new double[maxAlpha/dAlpha][maxAlpha/dAlpha];	
						
		cc = 0;
		for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {myAngles[cc] = angle; cc++;}
			
		//create xy coordinates of ellipse..
		xy = math.getXYOfEllipseFromAngle(myAngles, fE.xCenter, fE.yCenter, fE.majorRadius, fE.minorRadius, fE.theta);
		
		//cast it into two arrays..
		for (j = 0 ; j < myAngles.length ; j++) x[j] = (float)xy[j][0];
		for (j = 0 ; j < myAngles.length ; j++) y[j] = (float)xy[j][1]; 
					
		PolygonRoi pRoi = new PolygonRoi(x, y, Roi.POLYGON);		
		
		return pRoi;		
		
	}
	
	public PolygonRoi[] makeMeAPolygonRoiStack(String innerOrOuter, String wideTightOrExact, ObjectDetector.ColumnCoordinates preciseCC, int cutAwayFromWalls) {

		TailoredMaths math = new TailoredMaths();
		
		int i, j, cc;
		int maxAlpha = 360;
		int dAlpha = 5;	
		double[] myAngles = new double[maxAlpha/dAlpha];
		float[] x = new float[maxAlpha/dAlpha];
		float[] y = new float[maxAlpha/dAlpha];
		double[][] xy = new double[maxAlpha/dAlpha][maxAlpha/dAlpha];	
		double majRad = 0;
		double minRad = 0;
		double clipper = -100;
		
		PolygonRoi[] pRoi = new PolygonRoi[preciseCC.heightOfColumn];
		
		cc = 0;
		for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {myAngles[cc] = angle; cc++;}
		
		if (wideTightOrExact.contentEquals("wide")) clipper = 1;
		if (wideTightOrExact.contentEquals("tight")) clipper = -1;
		if (wideTightOrExact.contentEquals("exact")) clipper = 0;
		if (wideTightOrExact.contentEquals("manual")) clipper = cutAwayFromWalls;
		
		for (i = 0 ; i < preciseCC.heightOfColumn ; i++) { 
			
			if (innerOrOuter.contentEquals("inner")) {
				majRad = preciseCC.innerMajorRadius[i];
				minRad = preciseCC.innerMinorRadius[i];
			} 
			if (innerOrOuter.contentEquals("outer")) {
				majRad = preciseCC.outerMajorRadius[i];
				minRad = preciseCC.outerMinorRadius[i];
				clipper = -clipper;
			}
			if (innerOrOuter.contentEquals("outerFromInner")) {
				majRad = preciseCC.innerMajorRadius[i];
				minRad = preciseCC.innerMinorRadius[i];
				if (preciseCC.wallThickness[i] - 10 > 3) clipper = preciseCC.wallThickness[i] - clipper ;
				else clipper = 3;
			}
			
			xy = math.getXYOfEllipseFromAngle(myAngles, preciseCC.xmid[i], preciseCC.ymid[i], majRad + clipper, minRad + clipper, preciseCC.theta[i]);
		
			for (j = 0 ; j < myAngles.length ; j++) x[j] = (float)xy[j][0];
			for (j = 0 ; j < myAngles.length ; j++) y[j] = (float)xy[j][1]; 
					
			pRoi[i - preciseCC.topOfColumn] = new PolygonRoi(x, y, Roi.POLYGON);
		}
		
		return pRoi;
	}
	
	public PolygonRoi makeMeAnIndependentRoi(int[] imageDimensions, MenuWaiter.PoreSpaceAnalyserOptions mPSA) {

		TailoredMaths math = new TailoredMaths();
		
		int maxAlpha = 360;
		int dAlpha = 1;	
		double[] myAngles = new double[maxAlpha/dAlpha];
		float[] x = new float[maxAlpha/dAlpha];
		float[] y = new float[maxAlpha/dAlpha];	
		double[][] xy = new double[maxAlpha/dAlpha][2];	
		int xMid = imageDimensions[0] / 2;
		int yMid = imageDimensions[1] / 2;		
		
		if (mPSA.choiceOfRoi.equals("Cuboid")) {
			float halfEdgeLength = mPSA.edgeLengthOfRoi / 2;
			float[] subX = new float[]{halfEdgeLength, halfEdgeLength, -halfEdgeLength, -halfEdgeLength};
			float[] subY = new float[]{halfEdgeLength, -halfEdgeLength, -halfEdgeLength, halfEdgeLength};
			for (int i = 0 ; i < 4 ; i++) {
				x[i] = xMid - subX[i];
				y[i] = yMid - subY[i];
			}			
		} else if (mPSA.choiceOfRoi.equals("Cylinder")) {
			float radius = mPSA.edgeLengthOfRoi / 2;
			
			int cc = 0;
			for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {
				myAngles[cc] = angle; 
				cc++;
			}		
			xy = math.getXYOfEllipseFromAngle(myAngles, xMid, yMid, radius, radius, 0);
			
			for (int i = 0 ; i < myAngles.length ; i++) {
				x[i] = (float)xy[i][0];
				y[i] = (float)xy[i][1];
			}	
		}
		
		PolygonRoi pRoi = new PolygonRoi(x, y, Roi.POLYGON);
		
		return pRoi;
	}
	
	public PolygonRoi[] makeMeASplineRoiStack(String innerOrOuter, String wideTightOrExact, ObjectDetector.ColumnCoordinates preciseCC, double cutAway) {

		TailoredMaths math = new TailoredMaths();
		
		int i, j, cc;
		int maxAlpha = 360;
		int dAlpha = 5;	
		double[] myAngles = new double[maxAlpha/dAlpha];
		float[] x = new float[maxAlpha/dAlpha];
		float[] y = new float[maxAlpha/dAlpha];
		double[][] xy = new double[maxAlpha/dAlpha][maxAlpha/dAlpha];	
		double majRad, minRad;
		double clipper = -100;
		
		PolygonRoi[] pRoi = new PolygonRoi[preciseCC.heightOfColumn];
		
		cc = 0;
		for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {myAngles[cc] = angle; cc++;}
		
		if (wideTightOrExact.contentEquals("wide")) clipper = 5;
		if (wideTightOrExact.contentEquals("tight")) clipper = -3;
		if (wideTightOrExact.contentEquals("exact")) clipper = 0;
		if (wideTightOrExact.contentEquals("manual")) clipper = -cutAway;
		
		for (i = 0 ; i < preciseCC.heightOfColumn ; i++) { 
			
			if (innerOrOuter.contentEquals("inner")) {
				majRad = preciseCC.innerMajorRadius[i];
				minRad = preciseCC.innerMinorRadius[i];
			} 
			else {
				majRad = preciseCC.outerMajorRadius[i];
				minRad = preciseCC.outerMinorRadius[i];
			}
			
			xy = math.getXYOfEllipseFromAngle(myAngles, preciseCC.xmid[i], preciseCC.ymid[i], majRad + clipper, minRad + clipper, preciseCC.theta[i]);
		
			for (j = 0 ; j < myAngles.length ; j++) x[j] = (float)xy[j][0];
			for (j = 0 ; j < myAngles.length ; j++) y[j] = (float)xy[j][1]; 
					
			pRoi[i] = new PolygonRoi(x, y, Roi.POLYGON);
		}
		
		return pRoi;
	}

	public void showMeMyRoi(ImagePlus nowTiff, ImageProcessor myIP, PolygonRoi pRoi, int colorCode) {
		
		//ContrastEnhancer myBC = new ContrastEnhancer();
		//myBC.equalize(myIP);
		
		if (pRoi != null) {
			Overlay myO = new Overlay(pRoi);		
			myO.setStrokeColor(Color.YELLOW); // if 1
			if (colorCode == 0) myO.setStrokeColor(Color.WHITE);		
			if (colorCode == 2) myO.setStrokeColor(Color.BLUE);
			if (colorCode == 3) myO.setStrokeColor(Color.RED);		
			nowTiff.setOverlay(myO);
		} else nowTiff.setHideOverlay(true);			
		
		
		nowTiff.updateAndDraw();
		nowTiff.show();
		
		IJ.wait(2500);
	
	}

	public PolygonRoi[] makeMeAPolygonRoiStack4RawImage(String innerOrOuter, String wideTightOrExact, ObjectDetector.ColumnCoordinates preciseCC, double cutAway) {

		//init units
		TailoredMaths math = new TailoredMaths();
		
		//init variables
		int i, j, cc;
		int maxAlpha = 360;
		int dAlpha = 5;	
		double[] myAngles = new double[maxAlpha/dAlpha];
		float[] x = new float[maxAlpha/dAlpha];
		float[] y = new float[maxAlpha/dAlpha];
		double[][] xy = new double[maxAlpha/dAlpha][maxAlpha/dAlpha];	
		double majRad, minRad;
		double clipper = -100;
		
		PolygonRoi[] pRoi = new PolygonRoi[preciseCC.heightOfColumn];
		
		cc = 0;
		for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {myAngles[cc] = angle; cc++;}
		
		if (wideTightOrExact.contentEquals("wide")) clipper = 1;
		if (wideTightOrExact.contentEquals("tight")) clipper = -1;
		if (wideTightOrExact.contentEquals("exact")) clipper = 0;
		if (wideTightOrExact.contentEquals("manual")) clipper = -cutAway;
		
		for (i = preciseCC.topOfColumn ; i < preciseCC.topOfColumn + preciseCC.heightOfColumn ; i++) { 
			
			if (innerOrOuter.contentEquals("inner")) {
				majRad = preciseCC.innerMajorRadius[i];
				minRad = preciseCC.innerMinorRadius[i];
			} 
			else {
				majRad = preciseCC.outerMajorRadius[i];
				minRad = preciseCC.outerMinorRadius[i];
			}
			
			xy = math.getXYOfEllipseFromAngle(myAngles, preciseCC.xmid[i], preciseCC.ymid[i], majRad + clipper, minRad + clipper, preciseCC.theta[i]);
		
			for (j = 0 ; j < myAngles.length ; j++) x[j] = (float)xy[j][0];
			for (j = 0 ; j < myAngles.length ; j++) y[j] = (float)xy[j][1]; 
					
			pRoi[i - preciseCC.topOfColumn] = new PolygonRoi(x, y, Roi.POLYGON);
		}
		
		return pRoi;
	}

	public ImagePlus[] prepareDesiredRoi(int imageNumber, InputOutput.MyFolderCollection mFC, ImagePlus nowTiff, MenuWaiter.PoreSpaceAnalyserOptions mPSA) {
		
		//init units
		InputOutput jIO = new InputOutput();
		ObjectDetector jOD = new ObjectDetector();
		HistogramStuff hist = new HistogramStuff();
		RoiHandler roi = new RoiHandler();
		
		//init some imortant variables
		int startSlice = 0;
		int stopSlice = nowTiff.getNSlices();		
		ImagePlus[] outTiffs = new ImagePlus[2];		
		String imgName = mFC.myTiffs[imageNumber];
		
		//init surface corrected binaries..
		ImageStack corrected4SurfaceStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());		
		ImageStack corr4ThicknessAnalyses = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
				
		//and here we go..
		if (mPSA.choiceOfRoi.equals("RealSample")) {
			
			ImagePlus soilSurface = new ImagePlus();
			PolygonRoi[] pRoi = null;
						
			//select the correct gauge and surface files			
			int[] myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(imgName, mFC.myGaugePaths, mFC.mySurfaceFileNames);
			
			//read gauge file
			String nowGaugePath = mFC.myGaugePaths[myGandS[0]];
			ObjectDetector.ColumnCoordinates jCO = jIO.readGaugeFile(nowGaugePath);	
			pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, mPSA.cutAwayFromWall);
			
			//load surfaces
			soilSurface = jIO.openTiff3D(mFC.mySurfaceFolder + "\\" + mFC.mySurfaceFileNames[myGandS[1]]);
			
			//in case the roi is a sub-sample
			if (mPSA.cutAwayFromTop > 0) {
				
				//determine the depth of the soil surface				
				int topSurfaceLocation = jOD.findMedianSoilSurfacePosition(soilSurface);	
				
				startSlice = (int)Math.round(topSurfaceLocation + mPSA.cutAwayFromTop);
				stopSlice = (int)Math.round(topSurfaceLocation + mPSA.cutAwayFromTop + mPSA.heightOfRoi);
				
				for (int i = startSlice ; i < stopSlice ; i++) {
					
					nowTiff.setPosition(i + 1);
					ImageProcessor nowIP = nowTiff.getProcessor();
					IJ.showStatus("Compiling cut-out slice #" + (i + 1 - startSlice) + "/" + mPSA.heightOfRoi);
				
					ImageProcessor modIP = nowIP.duplicate();				
				
					modIP.setRoi(pRoi[i - 1]);
					modIP.setColor(0);
					modIP.fillOutside(pRoi[i - 1]);	
				
					corrected4SurfaceStack.addSlice(modIP);					
				}
				
				ImagePlus whyTheFuck = new ImagePlus();
				whyTheFuck.setStack(corrected4SurfaceStack);
				outTiffs[0] = whyTheFuck;
				outTiffs[1] = null; 			
			} 
			
			else {						// if the whole sample should be analyzed
			
				//assign soil top and bottom surfaces		
				ImageStack surStack = soilSurface.getStack();
				ImageProcessor topIP = surStack.getProcessor(1);		
				ImageProcessor botIP = surStack.getProcessor(2);
				
				//get max surface depressions				
				int[] topSurfHist = topIP.getHistogram();
				int maxTopDepression = hist.findMaxFromHistogram(topSurfHist);			
				int[] botSurfHist = botIP.getHistogram();
				int maxBotDepression = hist.findMaxFromHistogram(botSurfHist);
				
				for (int i = startSlice ; i < stopSlice ; i++) {
					
					nowTiff.setPosition(i + 1);
					ImageProcessor nowIP = nowTiff.getProcessor();
					
					if (i < maxTopDepression) {		
						
						IJ.showStatus("Correcting for soil top surface in slice #" + (i + 1) + "/" + maxTopDepression);				
						
						ImageProcessor modIP = nowIP.duplicate();
						ImageProcessor forThIP = nowIP.duplicate();
						
						for (int x = 0 ; x < nowIP.getWidth() ; x++) {
							for (int y = 0 ; y < nowIP.getHeight() ; y++) {					
								int nowPix = nowIP.getPixel(x, y);
								int surPix = topIP.getPixel(x, y);
								if (i <= surPix) modIP.putPixelValue(x, y, 0);
								else modIP.putPixelValue(x, y, nowPix);
								if (i <= surPix - 100) forThIP.putPixelValue(x, y, 0);		//to prevent underestimation of pore diameters reaching the top..
								else forThIP.putPixelValue(x, y, nowPix);	
							}			
						}			
						
						modIP.setRoi(pRoi[i]);
						modIP.setColor(0);
						modIP.fillOutside(pRoi[i]);			
						
						forThIP.setRoi(pRoi[i]);
						forThIP.setColor(0);
						forThIP.fillOutside(pRoi[i]);		
												
						corrected4SurfaceStack.addSlice(modIP);
						corr4ThicknessAnalyses.addSlice(forThIP);
						
					}
					
					else if (i > nowTiff.getNSlices() - maxBotDepression - 1) {
						
						IJ.showStatus("Correcting for soil bottom surface in slice #" + (i + 1 - maxBotDepression -  maxTopDepression) + "/" + maxBotDepression);				
						
						ImageProcessor modIP = nowIP.duplicate();
						ImageProcessor forThIP = nowIP.duplicate();
						
						for (int x = 0 ; x < nowIP.getWidth() ; x++) {
							for (int y = 0 ; y < nowIP.getHeight() ; y++) {					
								int nowPix = nowIP.getPixel(x, y);
								int surPix = nowTiff.getNSlices() - botIP.getPixel(x, y);
								if (i >= surPix) modIP.putPixelValue(x, y, 0);
								else modIP.putPixelValue(x, y, nowPix);
								if (i >= surPix + 100) forThIP.putPixelValue(x, y, 0);		//to prevent underestimation of pore diameters reaching the top..
								else forThIP.putPixelValue(x, y, nowPix);	
							}			
						}		
						
						modIP.setRoi(pRoi[i]);
						modIP.setColor(0);
						modIP.fillOutside(pRoi[i]);	
						
						forThIP.setRoi(pRoi[i]);
						forThIP.setColor(0);
						forThIP.fillOutside(pRoi[i]);						
						
						corrected4SurfaceStack.addSlice(modIP);
						corr4ThicknessAnalyses.addSlice(forThIP);
						
					}
					
					else {
						
						IJ.showStatus("Adding in-between slice #" + (i + 1 - maxBotDepression -  maxTopDepression) + "/" + (nowTiff.getNSlices() - maxBotDepression -  maxTopDepression));
						
						ImageProcessor modIP = nowIP.duplicate();			
						ImageProcessor forThIP = nowIP.duplicate();
						
						modIP.setRoi(pRoi[i]);
						modIP.setColor(0);
						modIP.fillOutside(pRoi[i]);	
						
						forThIP.setRoi(pRoi[i]);
						forThIP.setColor(0);
						forThIP.fillOutside(pRoi[i]);						
						
						corrected4SurfaceStack.addSlice(modIP);		
						corr4ThicknessAnalyses.addSlice(forThIP);
						
					}
				}			
				
				ImagePlus corrTiff = new ImagePlus();
				ImagePlus corrThTiff = new ImagePlus();
				corrTiff.setStack(corrected4SurfaceStack);
				corrThTiff.setStack(corr4ThicknessAnalyses);		
				
				outTiffs[0] = corrTiff;
				outTiffs[1] = corrThTiff;
			}			
		
		}
		
		if (!mPSA.choiceOfRoi.equals("RealSample")) {
			outTiffs[0] = nowTiff;
			outTiffs[1] = null; 	
		}
		
		//save the cut-out samples
		if (mPSA.cutAwayFromTop > 0 | mPSA.cutAwayFromWall > 0 | mPSA.choiceOfRoi.equals("RealSample")) {
			jIO.tiffSaver(mFC.myPreOutFolder, imgName, outTiffs[0]);
		}
		if (mPSA.cutAwayFromTop == 0 & mPSA.choiceOfRoi.equals("RealSample")) {
			String thicknessTiffLocation = mFC.myPreOutFolder + "//Tiff4ThicknessAnalyses";
			new File(thicknessTiffLocation).mkdir();
			jIO.tiffSaver(thicknessTiffLocation, imgName, outTiffs[1]);		
		} 
		
		return outTiffs;
	}
}