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

import java.io.File;

import org.apache.commons.math3.stat.StatUtils;

import SoilJ.tools.ObjectDetector.ColCoords3D;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

/** 
 * RoiHandler is a SoilJ class containing all subroutines for the creation of 
 * regions of interests (ROI) overlays or overlay stacks.
 * 
 * @author John Koestel
 *
 */

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
	
	public PolygonRoi[] makeMeAPolygonRoiStack(String innerOrOuter, String wideTightOrExact, ObjectDetector.ColCoords3D preciseCC, int cutAwayFromWalls) {

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
		double clipper0 = -100;		
		double clipper = clipper0;
		
		PolygonRoi[] nRoi = new PolygonRoi[preciseCC.heightOfColumn];
		
		cc = 0;
		for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {myAngles[cc] = angle; cc++;}
		
		if (wideTightOrExact.contentEquals("wide")) clipper0 = 1;
		if (wideTightOrExact.contentEquals("tight")) clipper0 = -1;
		if (wideTightOrExact.contentEquals("exact")) clipper0 = 0;
		if (wideTightOrExact.contentEquals("manual")) clipper0 = cutAwayFromWalls;
		
		for (i = 0 ; i < preciseCC.heightOfColumn ; i++) { 
			
			if (innerOrOuter.contentEquals("inner")) {
				majRad = preciseCC.innerMajorRadius[i];
				minRad = preciseCC.innerMinorRadius[i];
				clipper = clipper0;
			} 
			if (innerOrOuter.contentEquals("outer")) {
				majRad = preciseCC.outerMajorRadius[i];
				minRad = preciseCC.outerMinorRadius[i];
				clipper = -clipper0;
			}
			if (innerOrOuter.contentEquals("outerFromInner")) {
				majRad = preciseCC.innerMajorRadius[i];
				minRad = preciseCC.innerMinorRadius[i];
				if (preciseCC.wallThickness[i] - 10 > 3) clipper = preciseCC.wallThickness[i] - clipper0 ;
				else clipper = 3;
			}
			
			xy = math.getXYOfEllipseFromAngle(myAngles, preciseCC.xmid[i], preciseCC.ymid[i], majRad + clipper, minRad + clipper, preciseCC.theta[i]);
		
			for (j = 0 ; j < myAngles.length ; j++) x[j] = (float)xy[j][0];
			for (j = 0 ; j < myAngles.length ; j++) y[j] = (float)xy[j][1]; 
					
			nRoi[i - preciseCC.topOfColumn] = new PolygonRoi(x, y, Roi.POLYGON);
		}
		
		return nRoi;
	}
	
	public PolygonRoi makeMeAnIndependentRoi(int[] imageDimensions, MenuWaiter.PoreSpaceAnalyzerOptions mPSA) {

		TailoredMaths math = new TailoredMaths();
		
		int maxAlpha = 360;
		int dAlpha = 1;	
		double[] myAngles = new double[maxAlpha/dAlpha];		
		double[][] xy = new double[maxAlpha/dAlpha][2];
		PolygonRoi pRoi = null;
		
		if (mPSA.choiceOfRoi.equals("Cuboid")) {	
			
			float x1 = mPSA.cubeX1;
			float y1 = mPSA.cubeY1;
			float x2 = mPSA.cubeX2;
			float y2 = mPSA.cubeY2;
			
			if (mPSA.cubeX1 == 0) {
				x1 = 1;
			}
			
			if (mPSA.cubeY1 == 0) {
				y1 = 1;
			}
			
			if (mPSA.cubeX2 == 0) {
				x2 = imageDimensions[0];
			}
			
			if (mPSA.cubeY2 == 0) {
				y2 = imageDimensions[1];
			}
			
			float[] x = new float[]{x1, x2, x2, x1};
			float[] y = new float[]{y1, y1, y2, y2};	
			
			pRoi = new PolygonRoi(x, y, Roi.POLYGON);
			
		} else if (mPSA.choiceOfRoi.equals("Cylinder")) {
			
			int cc = 0;
			for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {
				myAngles[cc] = angle; 
				cc++;
			}		
			xy = math.getXYOfEllipseFromAngle(myAngles, mPSA.cylX, mPSA.cylY, mPSA.cylRadius, mPSA.cylRadius, 0);
			
			float[] x = new float[maxAlpha/dAlpha];
			float[] y = new float[maxAlpha/dAlpha];
			
			for (int i = 0 ; i < myAngles.length ; i++) {
				x[i] = (float)xy[i][0];
				y[i] = (float)xy[i][1];
			}	
			
			pRoi = new PolygonRoi(x, y, Roi.POLYGON);
		}
		
		return pRoi;
	}
	
	public PolygonRoi[] makeMeASplineRoiStack(String innerOrOuter, String wideTightOrExact, ObjectDetector.ColCoords3D preciseCC, double cutAway) {

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

	public PolygonRoi[] makeMeAPolygonRoiStack4RawImage(String innerOrOuter, String wideTightOrExact, ObjectDetector.ColCoords3D preciseCC, double cutAway) {

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

	public ImagePlus[] prepareDesiredRoi(InputOutput.MyFolderCollection mFC, ImagePlus nowTiff, MenuWaiter.PoreSpaceAnalyzerOptions mPSA) {
		
		//init units
		InputOutput jIO = new InputOutput();
		ObjectDetector jOD = new ObjectDetector();
		ImageManipulator jIM = new ImageManipulator();
		HistogramStuff hist = new HistogramStuff();
		RoiHandler roi = new RoiHandler();
		RollerCaster rC = new RollerCaster();
		
		//init some important variables
		int startSlice = 0;
		int stopSlice = nowTiff.getNSlices();		
		ImagePlus[] outTiffs = new ImagePlus[2];		
		String imgName = mFC.colName;
		double area = 0;
		
		//init surface corrected binaries..
		ImageStack corrected4SurfaceStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());		
		ImageStack corr4ThicknessAnalyses = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
				
		//and here we go..
		if (mPSA.choiceOfRoi.equals("RealSample")) {
			
			ImagePlus soilSurface = new ImagePlus();
						
			//select the correct gauge and surface files
			int[] myGandS = new int[2];
			if (mPSA.includeSurfaceTopography) myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(imgName, mFC.myInnerCircleFiles, mFC.mySurfaceFileNames);
			else myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(imgName, mFC.myInnerCircleFiles, null);
			
			//read InnerCircle file
			String nowGaugePath = mFC.myInnerCircleFiles[myGandS[0]];			
			ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();
			int versio = jIO.checkInnerCircleFileVersion(nowGaugePath);			
			if (versio == 0) jCO = jIO.readInnerCircleVer0(nowGaugePath);	
			else jCO = jIO.readInnerCircleVer1(nowGaugePath);
			
			int averageRadius = (int)Math.round((StatUtils.mean(jCO.innerMajorRadius) + StatUtils.mean(jCO.innerMinorRadius)) / 2);
			PolygonRoi[] pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, -mPSA.cutAwayFromWall);
			PolygonRoi[] iRoi = null;
			if (mPSA.cutAwayFromWall == 0) roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 10);
			if (mPSA.cutAwayFromCenter > 0) iRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, -(averageRadius - mPSA.cutAwayFromCenter));
			
			//calculate average area
			double[] avgRadius = new double[jCO.innerMajorRadius.length];
			for (int radiusFinder = 0 ; radiusFinder < jCO.innerMajorRadius.length ; radiusFinder++) {
				avgRadius[radiusFinder] = (jCO.innerMajorRadius[radiusFinder] + jCO.innerMinorRadius[radiusFinder]) / 2;
			}
			double myRadius = StatUtils.mean(avgRadius);
			area = myRadius * myRadius * Math.PI;
			
			//load surfaces
			if (mPSA.includeSurfaceTopography) soilSurface = jIO.openTiff3D(mFC.mySurfaceFolder + "\\" + mFC.mySurfaceFileNames[myGandS[1]]);
			
			//if neither soil surface should be included
			if ((mPSA.cutAwayFromTop > 0 & mPSA.cutAwayFromBottom > 0 ) | !mPSA.includeSurfaceTopography) {
				
				ImageStack[] tempStack = new ImageStack[2];
				ImageStack dimStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
				tempStack[0] = dimStack;tempStack[1] = dimStack;
				
				outTiffs = jIM.addColumnRegion2Stack(nowTiff, 0, nowTiff.getNSlices(), pRoi, iRoi);			
				outTiffs[1] = null;
				
			}
			
			//if only bottom voxels nbeed to be removed below surface
			if ((mPSA.cutAwayFromTop > 0 & mPSA.cutAwayFromBottom == 0 ) & mPSA.includeSurfaceTopography) {
				
				//assign soil top and bottom surfaces		
				ImageStack surStack = soilSurface.getStack();						
				ImageProcessor botIP = surStack.getProcessor(2);
							
				int[] botSurfHist = botIP.getHistogram();
				int maxBotDepression = hist.findMaxFromHistogram(botSurfHist);
				
				ImageStack[] tempStack = new ImageStack[2];
				ImageStack dimStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
				tempStack[0] = dimStack;tempStack[1] = dimStack;
				ImagePlus[] topHalf = jIM.addColumnRegion2Stack(nowTiff, 0, maxBotDepression, pRoi, iRoi);
				ImagePlus[] botHalf = jIM.removePoresBelowSurface(nowTiff,maxBotDepression, botIP, pRoi, iRoi);
				
				outTiffs = jIM.stackImagePlusArrays(topHalf, botHalf);
				
			}
			
			//if only top voxels nbeed to be removed below surface
			if ((mPSA.cutAwayFromTop == 0 & mPSA.cutAwayFromBottom > 0 ) & mPSA.includeSurfaceTopography) {
				
				//assign soil top and bottom surfaces		
				ImageStack surStack = soilSurface.getStack();
				ImageProcessor topIP = surStack.getProcessor(1);		
				
				//get max surface depressions				
				int[] topSurfHist = topIP.getHistogram();
				int maxTopDepression = hist.findMaxFromHistogram(topSurfHist);			
				
				ImagePlus[] topHalf = jIM.removePoresAboveSurface(nowTiff, maxTopDepression, topIP, pRoi, iRoi);
				ImagePlus[] botHalf = jIM.addColumnRegion2Stack(nowTiff,maxTopDepression, nowTiff.getNSlices(), pRoi, iRoi);
				
				outTiffs = jIM.stackImagePlusArrays(topHalf, botHalf);
				
			}
			
			//if both soil surfaces should be included
			if (mPSA.cutAwayFromTop == 0 & mPSA.includeSurfaceTopography & mPSA.cutAwayFromBottom == 0) {
				
				//assign soil top and bottom surfaces		
				ImageStack surStack = soilSurface.getStack();
				ImageProcessor topIP = surStack.getProcessor(1);		
				ImageProcessor botIP = surStack.getProcessor(2);
				
				//get max surface depressions				
				int[] topSurfHist = topIP.getHistogram();
				int maxTopDepression = hist.findMaxFromHistogram(topSurfHist);			
				int[] botSurfHist = botIP.getHistogram();
				int maxBotDepression = hist.findMaxFromHistogram(botSurfHist);
				
				ImagePlus[] topOfC = jIM.removePoresAboveSurface(nowTiff, maxTopDepression, topIP, pRoi, iRoi);
				ImagePlus[] middleOfC= jIM.addColumnRegion2Stack(nowTiff,maxTopDepression, maxBotDepression, pRoi, iRoi);
				ImagePlus[] botOfC = jIM.removePoresBelowSurface(nowTiff,maxBotDepression, botIP, pRoi, iRoi);
			
				outTiffs = jIM.stackImagePlusArrays(topOfC, middleOfC);
				outTiffs = jIM.stackImagePlusArrays(outTiffs, botOfC);				
			
			}			
		
			//cut image in XY-plane
			int minX = nowTiff.getWidth();
			int minY = nowTiff.getHeight();
			int maxX = 0;
			int maxY = 0;
			for (int i = 0 ; i < pRoi.length ; i++) {
				int[] nowX = pRoi[i].getXCoordinates();
				int[] nowY = pRoi[i].getXCoordinates();
				
				double nmiX = StatUtils.min(rC.castInt2Double(nowX)) + pRoi[i].getXBase();
				double nmiY = StatUtils.min(rC.castInt2Double(nowY)) + pRoi[i].getYBase();
				double nmaX = StatUtils.max(rC.castInt2Double(nowX)) + pRoi[i].getXBase();
				double nmaY = StatUtils.max(rC.castInt2Double(nowY)) + pRoi[i].getYBase();
				
				if (nmiX < minX) minX = (int)Math.round(nmiX);
				if (nmiY < minY) minY = (int)Math.round(nmiY);
				if (nmaX > maxX) maxX = (int)Math.round(nmaX);
				if (nmaY > maxY) maxY = (int)Math.round(nmaY);
				
			}
		
			//do the cut
			int[] imageDimensions = {maxX-minX, maxY-minY}; 
			mPSA.choiceOfRoi = "Cuboid";
			mPSA.cubeX1 = minX;
			mPSA.cubeY1 = minY;
			mPSA.cubeX2 = maxX;
			mPSA.cubeY2 = maxY;
			PolygonRoi cutRoi = makeMeAnIndependentRoi(imageDimensions, mPSA);			
			ImageStack cutStack0 = new ImageStack(maxX-minX, maxY-minY);
			ImageStack cutStack1 = new ImageStack(maxX-minX, maxY-minY);
		
			for (int i = 0 ; i < outTiffs[0].getNSlices() ; i++) {
				
				outTiffs[0].setPosition(i+1);
				ImageProcessor nowIP = outTiffs[0].getProcessor();
				
				nowIP.setRoi(cutRoi);
				ImageProcessor cutIP = nowIP.crop();
				
				cutStack0.addSlice(cutIP);
				
				if (outTiffs[1] != null) {
					
					outTiffs[1].setPosition(i+1);
					ImageProcessor now1IP = outTiffs[1].getProcessor();
					
					now1IP.setRoi(cutRoi);
					ImageProcessor cutIP1 = now1IP.crop();
					
					cutStack1.addSlice(cutIP1);					
				}
			}
			
			outTiffs[0].setStack(cutStack0);
			if (outTiffs[1] != null) outTiffs[1].setStack(cutStack1);
			
			mPSA.choiceOfRoi = "RealSample";
			
			//if there is a surface File, also cut this one..
			if (mPSA.includeSurfaceTopography) {
				
				ImageStack cutStackS = new ImageStack(maxX-minX, maxY-minY);
				
				for (int i = 0 ; i < 2 ; i++) {
					
					soilSurface.setPosition(i + 1);					
					ImageProcessor nowIP = soilSurface.getProcessor();
					
					nowIP.setRoi(cutRoi);
					ImageProcessor cutIP = nowIP.crop();
					
					cutStackS.addSlice(cutIP);
				}
				
				ImagePlus cutSurfacefile = new ImagePlus();
				cutSurfacefile.setStack(cutStackS);
				
				//save it
				jIO.tiffSaver(mFC.myCutSurfaceFolder, mFC.fileName, cutSurfacefile);				
			}
			
		}
		
		if (!mPSA.choiceOfRoi.equals("RealSample")) {
			
			if (mPSA.choiceOfRoi.equals("Cuboid")) {
				startSlice = (int)Math.round(mPSA.cubeZ1);
				stopSlice = (int)Math.round(mPSA.cubeZ2);
				area = (mPSA.cubeX2 - mPSA.cubeX1) *  (mPSA.cubeY2 - mPSA.cubeY1);
			}
			
			if (mPSA.choiceOfRoi.equals("Cylinder")) {
				startSlice = (int)Math.round(mPSA.cylZ1);
				stopSlice = (int)Math.round(mPSA.cylZ2);
				area = mPSA.cylRadius * mPSA.cylRadius * Math.PI; 
			}
			
			if (mPSA.choiceOfRoi.equals("Everything!") | mPSA.choiceOfRoi.equals("TopOfEveryThing") | mPSA.choiceOfRoi.equals("BottomOfEveryThing")) {
				
				if (mPSA.areaOfInterest == 0) area = nowTiff.getWidth() * nowTiff.getHeight();
				else area = mPSA.areaOfInterest;
				
				startSlice = 0;
				stopSlice = nowTiff.getNSlices();
				
				if (mPSA.choiceOfRoi.equals("TopOfEveryThing")) stopSlice = (int)Math.round(stopSlice / 2);
				if (mPSA.choiceOfRoi.equals("BottomOfEveryThing")) startSlice = (int)Math.round(stopSlice / 2);
				
			}
			
			int[] imageDimensions = {nowTiff.getWidth(), nowTiff.getHeight()};
			if (mPSA.cutCanvas) {
				imageDimensions[0] = mPSA.cubeX2 - mPSA.cubeX1; //at the moment only possible using the cuboid ROI
				imageDimensions[1] = mPSA.cubeY2 - mPSA.cubeY1; //at the moment only possible using the cuboid ROI
			}
			
			ImageStack roiStack = new ImageStack(imageDimensions[0], imageDimensions[1]);
			
			PolygonRoi pRoi = makeMeAnIndependentRoi(imageDimensions, mPSA);
			
			for (int i = startSlice ; i < stopSlice ; i++) {
				
				nowTiff.setPosition(i + 1);
				ImageProcessor nowIP = nowTiff.getProcessor();
				IJ.showStatus("Compiling cut-out slice #" + (i + 1 - startSlice));
			
				ImageProcessor modIP = nowIP.duplicate();				
			
				if (!mPSA.choiceOfRoi.equals("TopOfEveryThing") & !mPSA.choiceOfRoi.equals("BottomOfEveryThing") & !mPSA.choiceOfRoi.equals("EveryThing!")) {
					modIP.setRoi(pRoi);
					modIP.setColor(0);
					modIP.fillOutside(pRoi);
				}
				
				if (mPSA.cutCanvas) {
					ImageProcessor cutIP = modIP.crop();
					roiStack.addSlice(cutIP);
				}
				else roiStack.addSlice(modIP);					
			}
			
			ImagePlus whyTheFuck = new ImagePlus();
			whyTheFuck.setStack(roiStack);
			outTiffs[0] = whyTheFuck;
			outTiffs[1] = null;
		}
		
		//save the cut-out samples
		if (mPSA.cutAwayFromTop > 0 | mPSA.cutAwayFromWall > 0 | mPSA.choiceOfRoi.equals("RealSample")) {
			jIO.tiffSaver(mFC.myPreOutFolder, imgName + ".tif", outTiffs[0]);						
		}
		if (mPSA.cutAwayFromTop == 0 & mPSA.choiceOfRoi.equals("RealSample")) {
			String thicknessTiffLocation = mFC.myPreOutFolder + "//Tiff4ThicknessAnalyses";
			new File(thicknessTiffLocation).mkdir();
			jIO.tiffSaver(thicknessTiffLocation, imgName + ".tif", outTiffs[1]);
		} 
		
		if (!mPSA.choiceOfRoi.equals("RealSample")) {
			jIO.tiffSaver(mFC.myPreOutFolder, imgName + ".tif", outTiffs[0]);
		}
		
		//save also the cross-sectional area to allow for macroporosity calculations.. 
		mPSA.areaOfInterest = area;
		String path4Area = mFC.myPreOutFolder + "\\" + "area.asc";
		jIO.writeStringIntoAsciiFile(path4Area, "" +  String.format("%8.0f\t",area));
		return outTiffs;
	}
}