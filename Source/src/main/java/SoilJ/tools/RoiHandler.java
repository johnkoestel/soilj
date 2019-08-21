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
	
	public class ColumnRoi {
		
		public double area;
		public ObjectDetector.ColCoords3D jCO;
		public PolygonRoi[] pRoi;
		public PolygonRoi[] iRoi;
		public ImagePlus nowTiff;
		public ImagePlus surfaceNotCut;
		
	}
	
	public ObjectDetector.ColCoords3D scaleColumnCoordinates(double[] iDims, ObjectDetector.ColCoords3D inCo, MenuWaiter.PoreSpaceAnalyzerOptions mPSA, double scalingFactor, double[] oddVoxelContribution) {
		
		ObjectDetector jOD = new ObjectDetector();
		ObjectDetector.ColCoords3D outCo = jOD.new ColCoords3D();
		
		int cFT = mPSA.mRSO.cutAwayFromTop;
		int iHeight = mPSA.mRSO.heightOfRoi;
		
		//check if height of ROI is defined as whole or half column length
		if (iHeight == 0) {
			if (mPSA.mRSO.cutZPercent == true) {
				if ((mPSA.mRSO.cutAwayFromTop == 50 & mPSA.mRSO.cutAwayFromBottom == 0) | (mPSA.mRSO.cutAwayFromTop == 0 & mPSA.mRSO.cutAwayFromBottom == 50)) {							
					iHeight = inCo.xmid.length / 2;
					if (mPSA.mRSO.cutAwayFromTop == 50) cFT = inCo.xmid.length / 2 - 1; 
				}			
				else {
					iHeight = inCo.xmid.length - (mPSA.mRSO.cutAwayFromTop + mPSA.mRSO.cutAwayFromBottom) * inCo.xmid.length / 100;
					cFT = inCo.xmid.length - mPSA.mRSO.cutAwayFromTop * inCo.xmid.length / 100;
				}
			}
			else {
				iHeight = inCo.xmid.length - mPSA.mRSO.cutAwayFromTop - mPSA.mRSO.cutAwayFromBottom;
			}
		}
		
		double[] innerMajorRadius = new double[iHeight];
		double[] innerMinorRadius = new double[iHeight];
		double[] xmid = new double[iHeight]; 
	    double[] ymid = new double[iHeight];
		double[] ixmid = new double[iHeight];
		double[] iymid = new double[iHeight];
		double[] theta = new double[iHeight];	
		double[] itheta = new double[iHeight]; 
		
		for (int i = 0 ; i < iHeight ; i ++) {
			
			innerMajorRadius[i] = scalingFactor * (inCo.innerMajorRadius[cFT + i] - mPSA.mRSO.cutAwayFromWall) - oddVoxelContribution[0];
			innerMinorRadius[i] = scalingFactor * (inCo.innerMinorRadius[cFT + i] - mPSA.mRSO.cutAwayFromWall) - oddVoxelContribution[1];
			xmid[i] = iDims[0] / 2 - oddVoxelContribution[0];
			ymid[i] = iDims[1] / 2 - oddVoxelContribution[0];
			ixmid[i] = iDims[0] / 2 - oddVoxelContribution[0];
			iymid[i] = iDims[1] / 2 - oddVoxelContribution[0];
			theta[i] = inCo.theta[cFT + i]; 	
			itheta[i] = inCo.itheta[cFT + i]; 	
			
		}
		
		outCo.innerMajorRadius = innerMajorRadius;
		outCo.innerMinorRadius = innerMinorRadius;
		outCo.xmid = xmid;
		outCo.ymid = ymid;
		outCo.ixmid = ixmid;	
		outCo.iymid = iymid;
		outCo.theta = theta;
		outCo.itheta = itheta; 
		outCo.topOfColumn = 0;
		outCo.heightOfColumn = iHeight;
		
		return outCo;
		
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
	
	public PolygonRoi makeMeAnIndependentRoi(int[] imageDimensions, MenuWaiter.ROISelectionOptions mRSO) {

		TailoredMaths math = new TailoredMaths();
		
		int maxAlpha = 360;
		int dAlpha = 1;	
		double[] myAngles = new double[maxAlpha/dAlpha];		
		double[][] xy = new double[maxAlpha/dAlpha][2];
		PolygonRoi pRoi = null;
		
		if (mRSO.choiceOfRoi.equals("Cuboid")) {	
			
			float x1 = mRSO.cubeX1;
			float y1 = mRSO.cubeY1;
			float x2 = mRSO.cubeX2;
			float y2 = mRSO.cubeY2;
			
			if (mRSO.cubeX1 == 0) {
				x1 = 1;
			}
			
			if (mRSO.cubeY1 == 0) {
				y1 = 1;
			}
			
			if (mRSO.cubeX2 == 0) {
				x2 = imageDimensions[0];
			}
			
			if (mRSO.cubeY2 == 0) {
				y2 = imageDimensions[1];
			}
			
			float[] x = new float[]{x1, x2, x2, x1};
			float[] y = new float[]{y1, y1, y2, y2};	
			
			pRoi = new PolygonRoi(x, y, Roi.POLYGON);
			
		} else if (mRSO.choiceOfRoi.equals("Cylinder")) {
			
			int cc = 0;
			for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {
				myAngles[cc] = angle; 
				cc++;
			}		
			xy = math.getXYOfEllipseFromAngle(myAngles, mRSO.cylX, mRSO.cylY, mRSO.cylRadius, mRSO.cylRadius, 0);
			
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

	public ColumnRoi prepareDesiredRoi(InputOutput.MyFileCollection mFC, ImagePlus nowTiff, MenuWaiter.ROISelectionOptions mRSO) {
		
		String pathSep = "/";
		
		//init units
		InputOutput jIO = new InputOutput();
		ObjectDetector jOD = new ObjectDetector();
		ImageManipulator jIM = new ImageManipulator();
		HistogramStuff hist = new HistogramStuff();
		RoiHandler roi = new RoiHandler();
		RollerCaster rC = new RollerCaster();
		ColumnRoi colRoi = new ColumnRoi();
		
		//init some important variables
		int startSlice = 0;
		int stopSlice = nowTiff.getNSlices();		
		ImagePlus[] outTiffs = new ImagePlus[2];		
		String imgName = mFC.colName;
		double area = 0;
				
		//and here we go..
		if (mRSO.choiceOfRoi.equals("RealSample")) {
			
			ImagePlus soilSurface = new ImagePlus();
						
			//select the correct gauge and surface files
			int[] myGandS = new int[2];
			if (mRSO.includeSurfaceTopography) myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);
			else myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);
			
			//read InnerCircle file			
			//IJ.error(mFC.myInnerCircleFolder);
			//IJ.error(mFC.myInnerCircleFiles[myGandS[0]]);
			
			String nowGaugePath = mFC.myInnerCircleFiles[myGandS[0]];			
			ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();
			int versio = jIO.checkInnerCircleFileVersion(nowGaugePath);			
			if (versio == 0) jCO = jIO.readInnerCircleVer0(nowGaugePath);	
			else jCO = jIO.readInnerCircleVer1(nowGaugePath);
			
			//also pass on column outline coordinates
			colRoi.jCO = jCO;
			
			//create Rois
			int averageRadius = (int)Math.round((StatUtils.mean(jCO.innerMajorRadius) + StatUtils.mean(jCO.innerMinorRadius)) / 2);
			PolygonRoi[] pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, -mRSO.cutAwayFromWall);
			PolygonRoi[] iRoi = null;
			if (mRSO.cutAwayFromWall == 0) iRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 10);
			if (mRSO.cutAwayFromCenter > 0) iRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, -(averageRadius - mRSO.cutAwayFromCenter));
			
			colRoi.pRoi = pRoi;
			colRoi.iRoi = iRoi;
			
			//calculate average area
			double[] avgRadius = new double[jCO.innerMajorRadius.length];
			for (int radiusFinder = 0 ; radiusFinder < jCO.innerMajorRadius.length ; radiusFinder++) {
				avgRadius[radiusFinder] = (jCO.innerMajorRadius[radiusFinder] + jCO.innerMinorRadius[radiusFinder]) / 2;
			}
			double myRadius = StatUtils.mean(avgRadius) - mRSO.cutAwayFromWall;			
			area = myRadius * myRadius * Math.PI;
			
			//load surfaces
			if (mRSO.includeSurfaceTopography) {
				String surfPath = mFC.mySurfaceFolder + pathSep + mFC.mySurfaceFileNames[myGandS[1]];
				soilSurface = jIO.openTiff3D(surfPath);
			}
			
			//extract image phase that is to be analyzed		
			//if (mPSA.imagePhase2BeAnalyzed != 255) nowTiff = jIM.extractPhaseOfInterest(nowTiff, mPSA.imagePhase2BeAnalyzed, mPSA.nameOfAnalyzedPhase);				
						
			//if neither soil surface should be included
			if ((mRSO.cutAwayFromTop > 0 & mRSO.cutAwayFromBottom > 0 ) | !mRSO.includeSurfaceTopography) {
				
				ImageStack[] tempStack = new ImageStack[2];
				ImageStack dimStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
				tempStack[0] = dimStack;tempStack[1] = dimStack;
				
				outTiffs = jIM.addColumnRegion2Stack(nowTiff, 0, nowTiff.getNSlices(), pRoi, iRoi);			
				outTiffs[1] = null;
				
				//also cut ROI
				for (int i = mFC.startSlice ; i < mFC.stopSlice - mFC.startSlice ; i++) {
					colRoi.pRoi[i - mFC.startSlice] = pRoi[i];
					if (iRoi != null) colRoi.iRoi[i - mFC.startSlice] = iRoi[i];
				}	
				
			}
			
			//if the ROI is defined by the choice of layer below the surface and the height of the ROI
			if ((mRSO.cutAwayFromTop > 0 & mRSO.cutAwayFromBottom == 0 ) & mRSO.includeSurfaceTopography & mRSO.heightOfRoi > 0 & mRSO.heightOfRoi < 4000) {
				
				//assign soil top and bottom surfaces		
				ImageStack surStack = soilSurface.getStack();						
				ImageProcessor botIP = surStack.getProcessor(2);
							
				int[] botSurfHist = botIP.getHistogram();
				int maxBotDepression = mFC.nOfSlices - hist.findMaxFromHistogram(botSurfHist);
				
				ImageStack[] tempStack = new ImageStack[2];
				ImageStack dimStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
				tempStack[0] = dimStack;tempStack[1] = dimStack;
				
				//extract image phase that is to be analyzed		
				//if (mPSA.imagePhase2BeAnalyzed != 255) nowTiff = jIM.extractPhaseOfInterest(nowTiff, mPSA.imagePhase2BeAnalyzed, mPSA.nameOfAnalyzedPhase);
				
				//assemble
				ImagePlus[] topHalf = jIM.addColumnRegion2Stack(nowTiff, 0, maxBotDepression, pRoi, iRoi);
							
				outTiffs = topHalf;
				//outTiffs[1] = null;
				
				//outTiffs[0].updateAndDraw();
				//outTiffs[0].show();				
				
				//also cut ROI
				for (int i = mFC.startSlice ; i < mFC.stopSlice - mFC.startSlice ; i++) {
					colRoi.pRoi[i - mFC.startSlice] = pRoi[i];
					if (iRoi != null) colRoi.iRoi[i - mFC.startSlice] = iRoi[i];
				}				
			}
			
			//if only bottom voxels need to be removed below surface
			if ((mRSO.cutAwayFromTop > 0 & mRSO.cutAwayFromBottom == 0 ) & mRSO.includeSurfaceTopography & (mRSO.heightOfRoi == 0 | mRSO.heightOfRoi > 4000)) {
				
				//assign soil top and bottom surfaces		
				ImageStack surStack = soilSurface.getStack();						
				ImageProcessor botIP = surStack.getProcessor(2);
							
				int[] botSurfHist = botIP.getHistogram();
				int maxBotDepression = mFC.nOfSlices - hist.findMaxFromHistogram(botSurfHist);
				
				ImageStack[] tempStack = new ImageStack[2];
				ImageStack dimStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
				tempStack[0] = dimStack;tempStack[1] = dimStack;
				
				//extract image phase that is to be analyzed		
				//if (mPSA.imagePhase2BeAnalyzed != 255) nowTiff = jIM.extractPhaseOfInterest(nowTiff, mPSA.imagePhase2BeAnalyzed, mPSA.nameOfAnalyzedPhase);
				
				//assemble
				ImagePlus[] topHalf = jIM.addColumnRegion2Stack(nowTiff, 0, maxBotDepression, pRoi, iRoi);
				ImagePlus[] botHalf = jIM.removePoresBelowSurface(nowTiff,maxBotDepression, botIP, pRoi, iRoi);
				
				outTiffs = jIM.stackImagePlusArrays(topHalf, botHalf);
				
				//also cut ROI
				for (int i = mFC.startSlice ; i < mFC.stopSlice - mFC.startSlice ; i++) {
					colRoi.pRoi[i - mFC.startSlice] = pRoi[i];
					if (iRoi != null) colRoi.iRoi[i - mFC.startSlice] = iRoi[i];
				}					
			}
			
			//if only top voxels need to be removed below surface
			if ((mRSO.cutAwayFromTop == 0 & mRSO.cutAwayFromBottom > 0 ) & mRSO.includeSurfaceTopography & (mRSO.heightOfRoi == 0 | mRSO.heightOfRoi > 4000)) {
				
				//assign soil top and bottom surfaces		
				ImageStack surStack = soilSurface.getStack();
				ImageProcessor topIP = surStack.getProcessor(1);		
				
				//get max surface depressions				
				int[] topSurfHist = topIP.getHistogram();
				int maxTopDepression = mFC.nOfSlices - hist.findMaxFromHistogram(topSurfHist);			
				
				//extract image phase that is to be analyzed		
				//if (mPSA.imagePhase2BeAnalyzed != 255) nowTiff = jIM.extractPhaseOfInterest(nowTiff, mPSA.imagePhase2BeAnalyzed, mPSA.nameOfAnalyzedPhase);
				
				//assemble
				ImagePlus[] topHalf = jIM.removePoresAboveSurface(nowTiff, maxTopDepression, topIP, pRoi, iRoi);
				ImagePlus[] botHalf = jIM.addColumnRegion2Stack(nowTiff,maxTopDepression, nowTiff.getNSlices(), pRoi, iRoi);
				
				outTiffs = jIM.stackImagePlusArrays(topHalf, botHalf);
				
				//also cut ROI
				for (int i = mFC.startSlice ; i < mFC.stopSlice - mFC.startSlice ; i++) {
					colRoi.pRoi[i - mFC.startSlice] = pRoi[i];
					if (iRoi != null) colRoi.iRoi[i - mFC.startSlice] = iRoi[i];
				}		
				
			}
			
			//if both soil surfaces should be included
			if (mRSO.cutAwayFromTop == 0 & mRSO.includeSurfaceTopography & mRSO.cutAwayFromBottom == 0 & (mRSO.heightOfRoi == 0 | mRSO.heightOfRoi > 4000)) {
				
				//assign soil top and bottom surfaces		
				ImageStack surStack = soilSurface.getStack();
				ImageProcessor topIP = surStack.getProcessor(1);		
				ImageProcessor botIP = surStack.getProcessor(2);
				
				//get max surface depressions				
				int[] topSurfHist = topIP.getHistogram();
				int maxTopDepression = hist.findMaxFromHistogram(topSurfHist);			
				int[] botSurfHist = botIP.getHistogram();
				int maxBotDepression = mFC.nOfSlices - hist.findMaxFromHistogram(botSurfHist);
				
				//extract image phase that is to be analyzed		
				//if (mPSA.imagePhase2BeAnalyzed != 255) nowTiff = jIM.extractPhaseOfInterest(nowTiff, mPSA.imagePhase2BeAnalyzed, mPSA.nameOfAnalyzedPhase);
				
				//assemble
				ImagePlus[] topOfC = jIM.removePoresAboveSurface(nowTiff, maxTopDepression, topIP, pRoi, iRoi);
				ImagePlus[] middleOfC= jIM.addColumnRegion2Stack(nowTiff,maxTopDepression, maxBotDepression, pRoi, iRoi);
				ImagePlus[] botOfC = jIM.removePoresBelowSurface(nowTiff,maxBotDepression, botIP, pRoi, iRoi);
			
				outTiffs = jIM.stackImagePlusArrays(topOfC, middleOfC);
				outTiffs = jIM.stackImagePlusArrays(outTiffs, botOfC);		
				
				colRoi.pRoi = pRoi;
				colRoi.iRoi = iRoi;
			
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
			mRSO.choiceOfRoi = "Cuboid";
			mRSO.cubeX1 = minX;
			mRSO.cubeY1 = minY;
			mRSO.cubeX2 = maxX;
			mRSO.cubeY2 = maxY;
			PolygonRoi cutRoi = makeMeAnIndependentRoi(imageDimensions, mRSO);			
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
			
			mRSO.choiceOfRoi = "RealSample";
			
			//if there is a surface File, also cut this one..
			if (mRSO.includeSurfaceTopography & (mRSO.cutAwayFromTop > 0 | mRSO.cutAwayFromBottom > 0)) {
				
				ImageStack cutStackS = new ImageStack(maxX-minX, maxY-minY);
					
				int roiEnds = 0;
				for (int i = 0 ; i < 2 ; i++) {
					
					soilSurface.setPosition(i + 1);					
					ImageProcessor nowIP = soilSurface.getProcessor().duplicate();
					
					int[] myHist = nowIP.getHistogram();					
					int surfMedian = hist.findMedianFromHistogram(myHist);									
					int roiStartsAt = surfMedian + mRSO.cutAwayFromTop;		
					
					if (i == 0)	roiEnds = roiStartsAt + mRSO.heightOfRoi;
					
					if (i == 1) {
						roiStartsAt = mFC.nOfSlices - roiEnds;
					}
					
					nowIP.subtract(roiStartsAt);
					nowIP.add(1);  //make sure that the new surface is not 0. 
					
					nowIP.setRoi(cutRoi);
					ImageProcessor cutIP = nowIP.crop();
					
					//ImagePlus test = new ImagePlus("",cutIP);
					//test.updateAndDraw();
					//test.show();
					
					cutStackS.addSlice(cutIP);
				}
				
				ImagePlus cutSurfacefile = new ImagePlus();
				cutSurfacefile.setStack(cutStackS);
						
				//save it
				jIO.tiffSaver(mFC.myCutSurfaceFolder, mFC.fileName, cutSurfacefile);				
			}
			
		}
		
		if (!mRSO.choiceOfRoi.equals("RealSample")) {
			
			startSlice = 0;
			stopSlice = nowTiff.getNSlices();
			
			if (mRSO.choiceOfRoi.equals("Cuboid")) {
				area = (mRSO.cubeX2 - mRSO.cubeX1) *  (mRSO.cubeY2 - mRSO.cubeY1);
			}
			
			if (mRSO.choiceOfRoi.equals("Cylinder")) {
				area = mRSO.cylRadius * mRSO.cylRadius * Math.PI;
			}
			
			if (mRSO.choiceOfRoi.equals("Everything!") | mRSO.choiceOfRoi.equals("TopOfEveryThing") | mRSO.choiceOfRoi.equals("BottomOfEveryThing")) {
				
				if (mRSO.areaOfInterest == 0) area = nowTiff.getWidth() * nowTiff.getHeight();
				else area = mRSO.areaOfInterest;
				
				startSlice = 0;
				stopSlice = nowTiff.getNSlices();
				
				if (mRSO.choiceOfRoi.equals("TopOfEveryThing")) stopSlice = (int)Math.round(stopSlice / 2);
				if (mRSO.choiceOfRoi.equals("BottomOfEveryThing")) startSlice = (int)Math.round(stopSlice / 2);
				
			}
			
			int[] imageDimensions = {nowTiff.getWidth(), nowTiff.getHeight()};
			
			PolygonRoi pRoi = makeMeAnIndependentRoi(imageDimensions, mRSO);
						
			ImageStack roiStack = new ImageStack(imageDimensions[0], imageDimensions[1]);
			if (mRSO.cutCanvas) roiStack = new ImageStack(pRoi.getBounds().width, pRoi.getBounds().height);
			
			//IJ.showMessage(imageDimensions[0] + " --- " + imageDimensions[1]);
			//IJ.showMessage(mPSA.cubeX1 + " --- " + mPSA.cubeX2);
			//IJ.showMessage(mPSA.cubeY1 + " --- " + mPSA.cubeY2);
			
			for (int i = startSlice ; i < stopSlice ; i++) {
				
				nowTiff.setPosition(i + 1);
				ImageProcessor nowIP = nowTiff.getProcessor();
				IJ.showStatus("Compiling cut-out slice #" + (i + 1 - startSlice));
			
				ImageProcessor modIP = nowIP.duplicate();
				modIP.setRoi(pRoi);
			
				if (!mRSO.choiceOfRoi.equals("TopOfEveryThing") & !mRSO.choiceOfRoi.equals("BottomOfEveryThing") & !mRSO.choiceOfRoi.equals("EveryThing!")) {					
					modIP.setColor(0);
					modIP.fillOutside(pRoi);
				}
				
				if (mRSO.cutCanvas) {
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
		if (mRSO.saveROI) {
			if (mRSO.cutAwayFromTop > 0 | mRSO.cutAwayFromWall > 0 | mRSO.choiceOfRoi.equals("RealSample")) {			
				jIO.tiffSaver(mFC.myPreOutFolder, imgName + ".tif", outTiffs[0]);						
			}
			if (mRSO.cutAwayFromTop == 0 & mRSO.choiceOfRoi.equals("RealSample") & !mRSO.includeSurfaceTopography) {
				jIO.tiffSaver(mFC.myPreOutFolder, imgName + ".tif", outTiffs[0]);						
			}
			if (mRSO.cutAwayFromTop == 0 & mRSO.choiceOfRoi.equals("RealSample") & mRSO.includeSurfaceTopography) {
				String thicknessTiffLocation = mFC.myPreOutFolder + pathSep + "Tiff4ThicknessAnalyses";
				new File(thicknessTiffLocation).mkdir();
				jIO.tiffSaver(thicknessTiffLocation, imgName + ".tif", outTiffs[1]);
			} 
			
			if (!mRSO.choiceOfRoi.equals("RealSample")) {
				jIO.tiffSaver(mFC.myPreOutFolder, imgName + ".tif", outTiffs[0]);
			}
		}
		
		//save also the cross-sectional area to allow for macroporosity calculations.. 
		mRSO.areaOfInterest = area;
		colRoi.area = area;
		String path4Area = mFC.myPreOutFolder + pathSep + "area.asc";
		jIO.writeStringIntoAsciiFile(path4Area, "" +  String.format("%8.0f\t",area));
	
		colRoi.nowTiff = outTiffs[0];
		colRoi.surfaceNotCut = outTiffs[1];
		
		return colRoi;
	}
	
	public ColumnRoi assembleInnerCircleROIs(InputOutput.MyFileCollection mFC, ImagePlus nowTiff, boolean useInnerCircle, int extraInnerCircle, boolean cutCanvas) {
		
		//init units
		InputOutput jIO = new InputOutput();
		ObjectDetector jOD = new ObjectDetector();
		ImageManipulator jIM = new ImageManipulator();
		RoiHandler roi = new RoiHandler();		
		ColumnRoi colRoi = new ColumnRoi();
		
		if (useInnerCircle) {				
			//select the correct gauge and surface files
			int[] myGandS = new int[2];
			myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);
			
			//read InnerCircle file
			String nowGaugePath = mFC.myInnerCircleFiles[myGandS[0]];			
			ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();
			int versio = jIO.checkInnerCircleFileVersion(nowGaugePath);			
			if (versio == 0) jCO = jIO.readInnerCircleVer0(nowGaugePath);	
			else jCO = jIO.readInnerCircleVer1(nowGaugePath);
			
			//create Rois	
			PolygonRoi[] pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, -extraInnerCircle);			
			colRoi.pRoi = pRoi;
			
			//calculate ROI area
			double[] avgRadius = new double[jCO.innerMajorRadius.length];
			for (int radiusFinder = 0 ; radiusFinder < jCO.innerMajorRadius.length ; radiusFinder++) {
				avgRadius[radiusFinder] = (jCO.innerMajorRadius[radiusFinder] + jCO.innerMinorRadius[radiusFinder]) / 2;
			}
			double myRadius = StatUtils.mean(avgRadius);			
			double area = myRadius * myRadius * Math.PI;
			colRoi.area = area;
					
			//Cut out image
			ImagePlus outTiff = jIM.cutImageInXYPlane(nowTiff, pRoi, cutCanvas);	
			colRoi.nowTiff = outTiff;	
			
		}		
		else {
			
			colRoi.nowTiff = nowTiff;
			colRoi.area = nowTiff.getWidth() * nowTiff.getHeight();
			colRoi.pRoi = null;
			colRoi.jCO = null;
			
		}
		
		return colRoi;
	}
}