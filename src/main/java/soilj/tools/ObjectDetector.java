package SoilJ.tools;

import java.awt.Color;
import java.util.ArrayList;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.measure.CurveFitter;
import ij.plugin.ContrastEnhancer;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.plugin.filter.RankFilters;
import ij.process.ByteProcessor;
import ij.process.EllipseFitter;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import SoilJ.copiedTools.JParticleCounter;
import process3d.Dilate_;
import process3d.Erode_;
import SoilJ.tools.HistogramStuff.IlluminationInfo;

public class ObjectDetector implements PlugIn {

	
	public void run(String arg) {
				//ok, this is not needed..
	}
	
	public static final long MilliSecondsOfOneDay=86400000L;
	
	public class InnerOuterWallPerimeter2D {
		
		public double xCenter;
		public double yCenter;
		public double zCenter;
		public double outerMinorRadius;
		public double outerMajorRadius;				
		public double theta;		
		public double outerR2;
		
		public double ixCenter;
		public double iyCenter;
		public double innerMinorRadius;
		public double innerMajorRadius;				
		public double itheta;		
		public double innerR2;
		
		public double wallThickness;
		
		public boolean columnIsAtThisDepth;		
		public int outerWallNotFound;
		public int innerWallNotFound;
		
	}
	
	public class BestWallFinds {
		
		public float[] xD;
		public float[] yD;
		public double[] angle;		
		
	}
	
	public class SampleBrightness {
		
		int[][] histogram;
		int[][] centralHistogram;
		int[][] referenceHistogram;
		double[] wall;
		int[] q01;
		int[] q50;
		int[] q60;
		int[] q80;
		int[] q95;
		int[] central_q01;
		int[] central_q50;
		int[] central_q60;
		int[] central_q80;
		int[] central_q95;
		int[] ref_q01;
		int[] ref_q50;
		int[] ref_q60;
		int[] ref_q80;
		int[] ref_q95;
		
	}
	
	public class ColumnContainer {
		
		public ImagePlus nowTiff;
		public ImagePlus cutTiff;
		public ImagePlus normalizedTiff;
		public ImagePlus binaryTiff;
		
		public int startslice;
		public int stopslice;
		public int increment;		
		public int dAlpha;
		
		public ColumnCoordinates prelimCC;
		public ColumnCoordinates preciseCC;
		public ColumnCoordinates mendedCC;
			
		public MenuWaiter.ColumnFinderMenuReturn jCFS;
		
	}
	
	public class ApproximateColumnIllumination {
		
		int outside;
		int wall;
		int inside;
		int globalContrast;
		
	}
	
	public class RadialModes {
		
		double[][] maskedRadialModes;
		double[][] maskedRadialMinima;
		double[] radius;
		int[] maskingThreshold;
		
	}
	
	public class ColumnCoordinates {
	
		public double[] xmid;			//x midpoint
		public double[] ymid;			//y midpoint
		public double[] zmid;			//y midpoint
		
		public double[] ixmid;			//x midpoint (inner circle)
		public double[] iymid;			//y midpoint (inner circle)
		
		public double[] outerMajorRadius;
		public double[] innerMajorRadius;
		public double[] outerMinorRadius;
		public double[] innerMinorRadius;
		public double[] wallThickness;
		public double[] theta;  //angle of major ellipse axis
		public double[] itheta;  //angle of major ellipse axis (inner circle)
			
		public double tiltInXZ; 
		public double tiltInYZ;
		public double tiltTotal;
		
		public int illuminationLevel;
		
		public boolean[] columnIsInThisDepth;
		public boolean[] innerColumnFound;
		
		public int topOfColumn;
		public int heightOfColumn;
		public int bottomOfColumn;
		public int numberOfImputedLayers;

		public int bevelStartsAt;
		public FitStuff.LinearFunction bevelFunction;
		
		public double[] outerR2;
		public double[] innerR2;
		
		public int anglesChecked;
		
	}
	
	/*public double[][][] getRadialSteelGreyValues(ImagePlus nowTiff, ObjectDetector.ColumnCoordinates jCO, float radialMappingFactor) {
		
		//init variables
		int standardRadius = (int)Math.round(radialMappingFactor * jCO.approximateRadius);
		double[][][] radialGreyValues = new double[jCO.heightOfColumn][jCO.anglesChecked][standardRadius];
		int maxAlpha = 360;
		int dAlpha = maxAlpha / jCO.anglesChecked;
		
		
		for (int i = 0 ; i < nowTiff.getNSlices() ; i++) {
			
			IJ.showStatus("Sampling radial illumination of slice #" + (i + 1) + "/" + nowTiff.getNSlices());
			
			//set image to next slice
			nowTiff.setPosition(i+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			//sweep through all checked angles and get radial grey values
			int angleCounter = 0;			
			for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {				
								
				float dx = jCO.xID[i][angleCounter] - (float)jCO.xmid[i];
				float dy = jCO.yID[i][angleCounter] - (float)jCO.ymid[i];
				double nowRadius = Math.sqrt((double)(dx*dx) + (double)(dy*dy));
				int[] x = new int[(int)Math.round(nowRadius)];
				int[] y = new int[(int)Math.round(nowRadius)];
				float[] greyAtThisAngle = new float[(int)Math.round(nowRadius)];
				int[] radius = new int[(int)Math.round(nowRadius)];
				
				int distanceCounter = 0;
				for (double checkRadius = 0 ; checkRadius < (int)Math.round(nowRadius); checkRadius++) {					
					x[distanceCounter] = (int)Math.round(Math.sin(angle) * checkRadius + jCO.xmid[i]);
					y[distanceCounter] = (int)Math.round(Math.cos(angle) * checkRadius + jCO.ymid[i]);				
					greyAtThisAngle[distanceCounter] = (int)Math.round(nowIP.getPixelValue(x[distanceCounter], y[distanceCounter]));	
					radius[distanceCounter]	= (int)Math.round(checkRadius);
					distanceCounter++;					
				}	
				
				//stretch grey values to standardized radius
				float[] newRadius = new float[distanceCounter];
				for (int j = 0 ; j < distanceCounter ; j++) newRadius[j] = j * standardRadius / distanceCounter;
				SplineFitter sF = new SplineFitter(newRadius, greyAtThisAngle, radius.length);
				for (int j = 0 ; j < standardRadius ; j++) {
					double newGrey = sF.evalSpline(j);
					radialGreyValues[i][angleCounter][j] = newGrey;
				}
				
				angleCounter++;
				 				
			}
			
		}
		
		return radialGreyValues;
	}*/
	
	public SampleBrightness getSampleBrightness(ImagePlus nowTiff, ObjectDetector.ColumnCoordinates jCO) {
		
		//init units
		SampleBrightness myBrightness = new SampleBrightness();
		RoiHandler roi = new RoiHandler();
		HistogramStuff hist = new HistogramStuff();
		
		//init variables
		int[][] outHist = new int[jCO.heightOfColumn][65536];
		int[][] centralHist = new int[jCO.heightOfColumn][65536];
		int[][] referenceHist = new int[jCO.heightOfColumn][65536];
		int[] q01 = new int[jCO.heightOfColumn];
		int[] q50 = new int[jCO.heightOfColumn];
		int[] q60 = new int[jCO.heightOfColumn];
		int[] q80 = new int[jCO.heightOfColumn];
		int[] q95 = new int[jCO.heightOfColumn];
		int[] c01 = new int[jCO.heightOfColumn];
		int[] c50 = new int[jCO.heightOfColumn];
		int[] c60 = new int[jCO.heightOfColumn];
		int[] c80 = new int[jCO.heightOfColumn];
		int[] c95 = new int[jCO.heightOfColumn];
		int[] r01 = new int[jCO.heightOfColumn];
		int[] r50 = new int[jCO.heightOfColumn];
		int[] r60 = new int[jCO.heightOfColumn];
		int[] r80 = new int[jCO.heightOfColumn];
		int[] r95 = new int[jCO.heightOfColumn];
		
		
		//load polygon roi of inner perimeter		
		PolygonRoi[] innerRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 2);
		PolygonRoi[] centralRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, (int)Math.round(StatUtils.percentile(jCO.innerMinorRadius,50) - 100));
		PolygonRoi[] referenceRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, (int)Math.round(StatUtils.percentile(jCO.innerMinorRadius,50) - 200));
		myBrightness.wall = findMedianWallGreyValues(nowTiff, jCO);
				
		for (int i = 0; i < jCO.heightOfColumn ; i++) {
			
			IJ.showStatus("Sampling 1-D brightness profiles ... reached slice " + (i + 1) + "/" + nowTiff.getNSlices());
			
			//set image to next slice
			nowTiff.setPosition(i+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			//get median grey value of slice
			nowIP.setRoi(innerRoi[i]);
			int[] myHist = nowIP.getHistogram();
			for (int j = 0 ; j < myHist.length ; j++) outHist[i][j] = myHist[j];			
			q01[i] = hist.findQuantileFromHistogram(myHist, 0.01);	
			q50[i] = hist.findQuantileFromHistogram(myHist, 0.50);			
			q80[i] = hist.findQuantileFromHistogram(myHist, 0.80);
			nowIP.resetRoi();
			
			//get median from small inner circle
			nowIP.setRoi(centralRoi[i]);
			myHist = nowIP.getHistogram();
			for (int j = 0 ; j < myHist.length ; j++) centralHist[i][j] = myHist[j];			
			c01[i] = hist.findQuantileFromHistogram(myHist, 0.01);	
			c50[i] = hist.findQuantileFromHistogram(myHist, 0.50);			
			c80[i] = hist.findQuantileFromHistogram(myHist, 0.80);
			nowIP.resetRoi();
			
			//get reference illumination
			ImageProcessor copyIP = nowIP.duplicate();
			copyIP.setRoi(centralRoi[i]);
			copyIP.setValue(0);			
			copyIP.fill(centralRoi[i]);
			copyIP.setRoi(referenceRoi[i]);
			copyIP.setValue(0);	
			copyIP.fillOutside(referenceRoi[i]);
			copyIP.resetRoi();
			
			myHist = copyIP.getHistogram();
			myHist[0] = 0;
			for (int j = 0 ; j < myHist.length ; j++) referenceHist[i][j] = myHist[j];			
			r01[i] = hist.findQuantileFromHistogram(myHist, 0.01);	
			r50[i] = hist.findQuantileFromHistogram(myHist, 0.50);			
			r80[i] = hist.findQuantileFromHistogram(myHist, 0.80);

		}
		
		myBrightness.histogram = outHist;
		myBrightness.centralHistogram = centralHist;
		myBrightness.referenceHistogram = referenceHist;
		myBrightness.q01 = q01;
		myBrightness.q50 = q50;		
		myBrightness.q80 = q80;		
		myBrightness.central_q01 = c01;
		myBrightness.central_q50 = c50;		
		myBrightness.central_q80 = c80;
		myBrightness.ref_q01 = c01;
		myBrightness.ref_q50 = c50;		
		myBrightness.ref_q80 = c80;		
		
		return myBrightness;
		
	}
	
	public RadialModes getRadialPVCIllumination(ImagePlus nowTiff, ColumnCoordinates jCO, MenuWaiter.BeamDeHardeningReturn mBDH) {
					
		//init units
		RoiHandler roi = new RoiHandler();
		HistogramStuff hist = new HistogramStuff();
		
		//init variables
		int cc=0;
		RadialModes myRM = new RadialModes(); 		
		int[] myThresh = new int[jCO.heightOfColumn];
		RankFilters myRF = new RankFilters();
				
		//prepare Rois..
		PolygonRoi[] nowRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 2);
		double averageRadius = StatUtils.percentile(jCO.innerMinorRadius,50);
		for (int frisbee = 0 ; frisbee < (1 - mBDH.cutoff) * averageRadius  ; frisbee += mBDH.stepsize * averageRadius) cc++;
		double[][] maskedModeBrightness = new double[nowTiff.getNSlices()][cc + 2];
		double[][] maskedMinimumBrightness = new double[nowTiff.getNSlices()][cc + 2];
		
		PolygonRoi[][] allRoi = new PolygonRoi[jCO.heightOfColumn][cc + 2];
		for (int frisbee = 0 ; frisbee < cc ; frisbee++) {
			int nowFris = (int)Math.round(frisbee * mBDH.stepsize * averageRadius); 
			PolygonRoi[] frisbeeRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, nowFris);
			for (int i = 0 ; i < frisbeeRoi.length ; i++) allRoi[i][frisbee] = frisbeeRoi[i];			
		}
		double[] r = new double[cc + 2];
		r[0] = 0;											// add a pseudo point at 0 
		r[cc + 1] = averageRadius;  // and as well one at the column wall...
		for (int i = 0 ; i < cc ; i++) r[cc - i] = averageRadius - (mBDH.stepsize * averageRadius / 2 + i * mBDH.stepsize * averageRadius);		
		
		//sample illumination
		//for (int i = 599; i < 800 ; i++) {
	    for (int i = 0; i < jCO.heightOfColumn ; i++) {
			
			IJ.showStatus("Sampling radial illumination of slice #" + (i + 1) + "/" + nowTiff.getNSlices());

			//set image to next slice
			nowTiff.setPosition(i+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			//apply threshold to mask out the low density phase
			ImageProcessor maskedIP = nowIP.duplicate();
			maskedIP.setRoi(nowRoi[i]);	
			if (mBDH.maskingMethod != "None") {
				maskedIP.setAutoThreshold(mBDH.maskingMethod);
				myThresh[i] = maskedIP.getAutoThreshold();
			}
			else myThresh[i] = 0;
			
			//get illumination
			for (int frisbee = 0 ; frisbee < cc - 1; frisbee++) {			
				ImageProcessor ceoIP = maskedIP.duplicate();
				
				//apply a median filter
				if (mBDH.radiusOfMedianFilter > 1) myRF.rank(ceoIP, mBDH.radiusOfMedianFilter, RankFilters.MEDIAN);
				
				//create donut shaped image segment
				ceoIP.setValue(0);	
				ceoIP.fillOutside(allRoi[i][frisbee]);
				ceoIP.fill(allRoi[i][frisbee+1]);
				ceoIP.resetRoi();
				
				//get histogram 
				int[] nowHist = ceoIP.getHistogram();
				
				//sample minimum grey values as a proxy for the air phase brightness
				maskedMinimumBrightness[i][cc - frisbee] = hist.findNonZeroMinimumFromHistogram(nowHist, (int)mBDH.airPhaseClassMemberMinimum);
				
				//set all values below the threshold to 0.. so that only the soil matrix is left in the images 
				for (int j = 0 ; j < myThresh[i] ; j++) nowHist[j] = 0;
				
				//pick out mode grey value
				maskedModeBrightness[i][cc - frisbee] = hist.findNonZeroModeFromHistogram(nowHist);
				
	/*			if (i == 100) {
					ImagePlus testImg = new ImagePlus("Slice",ceoIP);
					testImg.updateAndDraw();
					testImg.show();
				}*/
			}
			//and also do the same for the central area ..
			ImageProcessor ceoIP = maskedIP.duplicate();
			ceoIP.setValue(0);	
			ceoIP.fillOutside(allRoi[i][cc - 1]);
			ceoIP.resetRoi();
			int[] nowHist = ceoIP.getHistogram();	
			maskedMinimumBrightness[i][1] = hist.findNonZeroMinimumFromHistogram(nowHist, (int)mBDH.airPhaseClassMemberMinimum);
			for (int j = 0 ; j < myThresh[i] ; j++) nowHist[j] = 0;
			maskedModeBrightness[i][1] = hist.findNonZeroModeFromHistogram(nowHist);	
			
			//assign the neighboring grey values to the two pseudo points
			maskedMinimumBrightness[i][0] = maskedMinimumBrightness[i][1];
			maskedMinimumBrightness[i][cc + 1] = maskedMinimumBrightness[i][cc];
			maskedModeBrightness[i][0] = maskedModeBrightness[i][1];
			maskedModeBrightness[i][cc + 1] = maskedModeBrightness[i][cc];
			
			//fix skip
			double skipTarget = maskedModeBrightness[i][cc - mBDH.skip];
			for (int j = cc - mBDH.skip + 1 ; j < cc + 2 ; j++) maskedModeBrightness[i][j] = skipTarget;
		}
		
		myRM.maskedRadialModes = maskedModeBrightness;
		myRM.maskedRadialMinima = maskedMinimumBrightness;
		myRM.radius = r;
		myRM.maskingThreshold = myThresh;
		
		return myRM;
					
	}
	
	public ApproximateColumnIllumination probeApproximateColumnIllumination(ColumnContainer colCon) {
		
		Median jMed = new Median();
		RoiHandler roi = new RoiHandler();
		FitStuff fit = new FitStuff(); 
		HistogramStuff hist = new HistogramStuff();
		ApproximateColumnIllumination medACI = new ApproximateColumnIllumination();
		DisplayThings disp = new DisplayThings();
		
		//define in which part of the image to search for the wall
		int imageHeight = colCon.nowTiff.getNSlices();
		int startCheck = (int)Math.round((double)imageHeight/5);
		int stopCheck = (int)Math.round((double)imageHeight/5*4);
		ApproximateColumnIllumination[] myACI = new ApproximateColumnIllumination[stopCheck - startCheck];
		
		//calculate median column positions
		double oXmid = jMed.evaluate(colCon.prelimCC.xmid);
		double oYmid = jMed.evaluate(colCon.prelimCC.ymid);	
		double oMajR = jMed.evaluate(colCon.prelimCC.outerMajorRadius);
		double oMinR = jMed.evaluate(colCon.prelimCC.outerMinorRadius);
		double otheta = jMed.evaluate(colCon.prelimCC.theta);		
		FitStuff.FittedEllipse fE = fit.new FittedEllipse();
				
		//create rois
		int safetyNet = 10;
		float imageDims = (colCon.nowTiff.getWidth() + colCon.nowTiff.getHeight()) / 2;
		float thickestWallOnEarth = imageDims / 12;
		float thinnestWallOnEarth = imageDims / 30;
		
		//init grey levels
		double[] inside = new double[stopCheck - startCheck];
		double[] wall = new double[stopCheck - startCheck];
		double[] outside = new double[stopCheck - startCheck];
		double[] globalContrast = new double[stopCheck - startCheck];
		
		//mercury
		fE.xCenter = oXmid; fE.yCenter = oYmid; fE.majorRadius = oMajR - thickestWallOnEarth; fE.minorRadius = oMinR - thickestWallOnEarth; fE.theta = otheta;
		PolygonRoi mercury = roi.makeRoiFromFittedEllipse(fE);
		
		//venus
		fE.xCenter = oXmid; fE.yCenter = oYmid; fE.majorRadius = oMajR - thinnestWallOnEarth; fE.minorRadius = oMinR - thinnestWallOnEarth; fE.theta = otheta;
		PolygonRoi venus = roi.makeRoiFromFittedEllipse(fE);
		
		//earth
		fE.xCenter = oXmid; fE.yCenter = oYmid; fE.majorRadius = oMajR - safetyNet; fE.minorRadius = oMinR - safetyNet; fE.theta = otheta;
		PolygonRoi earth = roi.makeRoiFromFittedEllipse(fE);

		//mars
		fE.xCenter = oXmid; fE.yCenter = oYmid; fE.majorRadius = oMajR + safetyNet; fE.minorRadius = oMinR + safetyNet; fE.theta = otheta;
		PolygonRoi mars = roi.makeRoiFromFittedEllipse(fE);
		
		//universe
		Roi all = new Roi(1, 1, colCon.nowTiff.getWidth(), colCon.nowTiff.getHeight()); 
		
		//probe the illuminations..
		for (int i = startCheck ; i < stopCheck ; i++) {
			
			colCon.nowTiff.setPosition(i+1);
			ImageProcessor nowIP = colCon.nowTiff.getProcessor();
			nowIP.setValue(0);
			nowIP.setBackgroundValue(0);
			
			//get soil illumination
			ImageProcessor soilIP = nowIP.duplicate();			
			soilIP.fillOutside(mercury);			
			inside[i-startCheck] = hist.findMedianFromHistogram(soilIP.getHistogram());
			//disp.displayIP(soilIP, "soil");			
			
			//get wall illumination
			ImageProcessor wallIP = nowIP.duplicate();			
			wallIP.fillOutside(earth);
			wallIP.fill(venus);
			wall[i-startCheck] = hist.findMedianFromHistogram(wallIP.getHistogram());
			//disp.displayIP(wallIP, "wall");			
			
			//get outside illumination
			ImageProcessor outsiIP = nowIP.duplicate();			
			outsiIP.fill(mars);	
			outsiIP.setRoi(all);
			outside[i-startCheck] = hist.findMedianFromHistogram(outsiIP.getHistogram());			
			//disp.displayIP(outsiIP, "outside");
			
			//getGlobalContrast			
			IlluminationInfo jII = probeIlluminationLevel(nowIP, colCon, false);
			globalContrast[i-startCheck] = jII.quantile90 - jII.quantile10;
			
		}
		
		//calculate medians		
		medACI.inside = (int)Math.round(jMed.evaluate(inside));
		medACI.wall = (int)Math.round(jMed.evaluate(wall));
		medACI.outside = (int)Math.round(jMed.evaluate(outside));
		medACI.globalContrast = (int)Math.round(jMed.evaluate(globalContrast));
		
		return medACI;
	}
	
	public HistogramStuff.IlluminationInfo probeIlluminationLevel(ImageProcessor nowIP, ColumnContainer colCon, boolean cutAwaySoil) {
		
		//init objects
		HistogramStuff hist = new HistogramStuff();
		RoiHandler roi = new RoiHandler();
		FitStuff fit = new FitStuff();
		
		//init varis
		HistogramStuff.IlluminationInfo jII = hist.new IlluminationInfo();
		ImageProcessor iIP = nowIP.duplicate();
		
		//prepare roi on the wall .. make a roi that is approximately on the center of the column wall..
		if (cutAwaySoil == true) {									//if prelimCC is not empty
			FitStuff.FittedEllipse fE = fit.new FittedEllipse();		
			fE.xCenter = StatUtils.percentile(colCon.prelimCC.xmid,50);
			fE.yCenter = StatUtils.percentile(colCon.prelimCC.ymid,50);
			fE.majorRadius = colCon.jCFS.widthOfWallRelativeToOuterRadius * StatUtils.percentile(colCon.prelimCC.outerMajorRadius,50);
			fE.minorRadius = colCon.jCFS.widthOfWallRelativeToOuterRadius * StatUtils.percentile(colCon.prelimCC.outerMinorRadius,50);
			fE.theta = StatUtils.percentile(colCon.prelimCC.theta,50);
		
			PolygonRoi wallRoi = roi.makeRoiFromFittedEllipse(fE);
				
			//cut away soil to get grey values of wall and surrounding, only...			
			iIP.setRoi(wallRoi);
			iIP.setColor(Color.black);
			iIP.fill(wallRoi);
			iIP.resetRoi();
			
			//DisplayThings dT = new DisplayThings();
			//dT.displayIP(iIP, "");			
		}	
		
		//get histogram of image
		int[] myHist = iIP.getHistogram();
		myHist[0] = 0;		//set 0 entries to zero
				
		//assign values	
		jII.min = hist.findQuantileFromHistogram(myHist, 0.00001);
		jII.quantile01 = hist.findQuantileFromHistogram(myHist, 0.01);
		jII.quantile05 = hist.findQuantileFromHistogram(myHist, 0.05);
		jII.quantile10 = hist.findQuantileFromHistogram(myHist, 0.1);
		jII.lowerQuartile = hist.findQuantileFromHistogram(myHist, 0.25);
		jII.median = hist.findQuantileFromHistogram(myHist, 0.5);
		jII.upperQuartile = hist.findQuantileFromHistogram(myHist, 0.75);
		jII.quantile90 = hist.findQuantileFromHistogram(myHist, 0.9);
		jII.quantile95 = hist.findQuantileFromHistogram(myHist, 0.95);
		jII.quantile99 = hist.findQuantileFromHistogram(myHist, 0.99);
		jII.max = hist.findQuantileFromHistogram(myHist, 0.999999);
		
		jII.specialQuantile = hist.findQuantileFromHistogram(myHist, colCon.jCFS.minWallGreyValue);
		
		jII.mean = (int)Math.round(hist.findMeanFromHistogram(myHist));
		
		return jII;		
	}
	
	public ColumnCoordinates findOrientationOfPVCOrAluColumn(ImagePlus nowTiff, MenuWaiter.ColumnFinderMenuReturn jCFS) {

		ColumnCoordinates pCC = new ColumnCoordinates();
		FitStuff fit = new FitStuff();
		TailoredMaths math = new TailoredMaths();
		
		//init outline related variables
		int i, j, k;
		int maxAlpha = 360;
		int dAlpha = 5;		
		int checkRange = (nowTiff.getWidth()/2 + nowTiff.getHeight()/2) / 3 ;		
		double xmid = nowTiff.getWidth() / 2;
		double ymid = nowTiff.getHeight() / 2;
		double radius = nowTiff.getHeight() / 2;
		int footprintOfMedianFilter = jCFS.medianFilter;
		int iCounter = 0;
		
		//find slices corresponding to the 30 and 70 percentiles and 
		double imageHeight = nowTiff.getNSlices();
		int startSlice = (int)(imageHeight * 0.25) + 1;
		int stopSlice = (int)(imageHeight * 0.75) + 1;
		int samplingNumber = 50;
		int dSlice = (stopSlice - startSlice) / samplingNumber;
		
		//init outline related vectors
		double[] x = new double[checkRange];
		double[] y = new double[checkRange];	
		double[] greyAtThisAngle = new double[checkRange];
		double[] angleAtThisAngle = new double[maxAlpha/dAlpha];
		float[] xOD = new float[maxAlpha/dAlpha]; // x of outer diameter
		float[] yOD = new float[maxAlpha/dAlpha]; // y of outer diameter
				
		//init out vectors
		double confidenceLevel = 0.99;
		double[] xCenter = new double[samplingNumber];
		double[] yCenter = new double[samplingNumber];
		double[] zCenter = new double[samplingNumber];
		double[] majorRadius = new double[samplingNumber];
		double[] minorRadius = new double[samplingNumber];
		double[] xyAngle = new double[samplingNumber];
		double[] R2 = new double[samplingNumber];
		double[] tilt = new double[3];
		
		double[] myMeanIlluminationLevel = new double[samplingNumber];
								
		for (i = startSlice ; i < startSlice + samplingNumber * dSlice ; i = i + dSlice) {  
			
			IJ.showStatus("Searching for location and orientation of column..");
							
			int cc=0;
			nowTiff.setPosition(i);				
			ImageProcessor myIP = nowTiff.getProcessor();
			
			//probe illumination level of background
			ColumnContainer colCon = new ColumnContainer();
			colCon.jCFS = jCFS;
			HistogramStuff.IlluminationInfo jII = probeIlluminationLevel(myIP, colCon, false);
			int myContrast = jII.quantile90 - jII.quantile10;
			int minAirPVCGradient = (int)Math.round(jCFS.airWallContrast * (double)myContrast);
				
			//sample grey values around the column
			for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {				
				
				int cc0=1;
				angleAtThisAngle[cc] = angle;
				for (double checkRadius = radius - 1; checkRadius > radius - checkRange ; checkRadius--) {					
					x[cc0] = Math.sin(angle) * checkRadius + xmid;
					y[cc0] = Math.cos(angle) * checkRadius + ymid;				
					greyAtThisAngle[cc0] = myIP.getPixelValue((int)x[cc0], (int)y[cc0]);					
					cc0++;
				}
				greyAtThisAngle[0] = 65000; 	//set first entry to a very large value so that there is no gradient (it would be 0 otherwise) 
				
				//apply medianFilter
				double[] mGreyAtThisAngle = math.oneDMedianFilter(greyAtThisAngle, footprintOfMedianFilter);
				
				//calculate the first derivative, dMG
				double[] dMG = math.oneDFirstDerivative(mGreyAtThisAngle, footprintOfMedianFilter);
						
				//look for the outer edge
				for (j = 0 ; j < checkRange - 3; j++) {
					if (dMG[j] > minAirPVCGradient) {
						int myJ = j;
						for (k = j ; k < j + 3 ; k++) {
							if (dMG[k] > dMG[myJ]) myJ = k;
						}		
						xOD[cc] = (float)x[myJ];
						yOD[cc] = (float)y[myJ];
						break;
					}
				}
				
				//increase vector index
				cc++;
			}
			
			//kick out values where the wall was not found
			int zeroCounter = 0;
			cc = 0;
			for (j = 0 ; j < xOD.length ; j++) if (xOD[j] == 0 & yOD[j] == 0) zeroCounter++;
			float[] xOR = new float[xOD.length - zeroCounter]; // x of outer diameter and zeros removed
			float[] yOR = new float[yOD.length - zeroCounter]; // y of outer diameter and zeros removed
			double[] angleR = new double[yOD.length - zeroCounter];  //also make an angle that has equally many entries..
			for (j = 0 ; j < xOD.length ; j++) if (xOD[j] == 0 & yOD[j] == 0) {} // do nothing
			else {
				xOR[cc] = xOD[j];
				yOR[cc] = yOD[j];
				angleR[cc] = angleAtThisAngle[j];
				cc++;
			}		
			
			//fit an elliptic ROI to the x and y
			EllipseFitter jEF = new EllipseFitter();			
			PolygonRoi pRoi = new PolygonRoi(xOR, yOR, Roi.POLYLINE); // create pRoi with outer Wall coordinates
			myIP.setRoi(pRoi); 
			jEF.fit(myIP, myIP.getStatistics());	
			
			//also kick out outliers 
			FitStuff.GoodnessOfFit gOF = fit.calculateR2OfEllipsoidFit(angleR, xOR, yOR, jEF);
			double maxDeviationAllowed = 0.05;
			ArrayList<Integer> badOnes = new ArrayList<Integer>();
			for (j = 0 ; j < xOR.length ; j++) if (Math.abs((gOF.d2m[j] - gOF.median2d2m) / gOF.median2d2m) > maxDeviationAllowed) badOnes.add(j);
			float[] xOS = new float[xOR.length - badOnes.size()]; // x of outer diameter and zeros and outliers removed
			float[] yOS = new float[yOR.length - badOnes.size()]; // y of outer diameter and zeros and outliers removed
			double[] angleS = new double[yOR.length - badOnes.size()];  //also make an angle that has equally many entries..
			cc = 0;
			for (j = 0 ; j < xOR.length ; j++) if (Math.abs((gOF.d2m[j] - gOF.median2d2m) / gOF.median2d2m) < maxDeviationAllowed) {
				xOS[cc] = xOR[j];
				yOS[cc] = yOR[j];
				angleS[cc] = angleR[j];
				cc++;
			}
			
			//do the fit again but without outliers
			PolygonRoi sRoi = new PolygonRoi(xOS, yOS, Roi.POLYLINE); // create pRoi with outer Wall coordinates
			myIP.setRoi(sRoi); 
			EllipseFitter jEFS = new EllipseFitter();
			jEFS.fit(myIP, myIP.getStatistics());			
			gOF = fit.calculateR2OfEllipsoidFit(angleS, xOS, yOS, jEFS);
			
			//transfer fitting results to vectors for export
			xCenter[iCounter] = jEFS.xCenter;
			yCenter[iCounter] = jEFS.yCenter;
			zCenter[iCounter] = i;
			minorRadius[iCounter] = jEFS.minor / 2;
			majorRadius[iCounter] = jEFS.major / 2;				
			xyAngle[iCounter] = jEFS.theta;			
			R2[iCounter] = gOF.R2;
			
			iCounter++;
		}
		
		//kick out fits that did not work
		int cc = 0;
		for (i = 0 ; i < R2.length ; i++) if (R2[i] > confidenceLevel) cc++;
		double[] xFiltered = new double[cc];
		double[] yFiltered = new double[cc];
		double[] zFiltered = new double[cc];
		double[] oMajRF = new double[cc];
		double[] oMinRF = new double[cc];
		double[] thetaF = new double[cc];
		double[] R2F = new double[cc];
		
		cc = 0;
		for (i = 0 ; i < R2.length ; i++) {
			if (R2[i] > confidenceLevel) {
				xFiltered[cc]=xCenter[i];
				yFiltered[cc]=yCenter[i];
				zFiltered[cc]=zCenter[i];
				oMajRF[cc]=majorRadius[i];
				oMinRF[cc]=minorRadius[i];
				thetaF[cc]=xyAngle[i];
				R2F[cc]=R2[i];
				cc++;
			}
		}
		
		//find column Tilt!
		tilt = findColumnTilt(xFiltered, yFiltered, zFiltered);
		
		//save results in struct
		pCC.xmid = xFiltered;
		pCC.ymid = yFiltered;
		pCC.zmid = zFiltered;
		pCC.outerMinorRadius = oMinRF;
		pCC.outerMajorRadius = oMajRF;
		pCC.theta = thetaF;
		pCC.outerR2 = R2F;
		pCC.tiltInXZ = tilt[0];
		pCC.tiltInYZ = tilt[1];
		pCC.tiltTotal = tilt[2];
		pCC.illuminationLevel = (int)StatUtils.mean(myMeanIlluminationLevel);
		
		return pCC;

	}
		
	public ColumnCoordinates findOrientationSteelColumn(ImagePlus nowTiff) {

		//init objects
		ColumnCoordinates pCC = new ColumnCoordinates();
		HistogramStuff hist = new HistogramStuff();
		TailoredMaths math = new TailoredMaths();
		FitStuff fit = new FitStuff();
		
		//init outline related variables		
		int i, j;
		int maxAlpha = 360;
		int dAlpha = 5;		
		int checkRange = nowTiff.getWidth() / 4;
		double xmid = nowTiff.getWidth() / 2;
		double ymid = nowTiff.getHeight() / 2;
		double radius = nowTiff.getHeight() / 2;				
		int footprintOfMedianFilter = 5;
		int iCounter = 0;
		int gradientSearchWindow = 5;		
		double steelGreyThreshold;
		
		//find slices corresponding to the 30 and 70 percentiles and 
		double imageHeight = nowTiff.getNSlices();
		int startSlice = (int)(imageHeight * 0.3) + 1;
		int stopSlice = (int)(imageHeight * 0.7) + 1;
		int samplingNumber = 30;
		int dSlice = (stopSlice - startSlice) / samplingNumber;		
		
		//init outline related vectors		
		double[] x = new double[checkRange];
		double[] y = new double[checkRange];	
		double[] greyAtThisAngle = new double[checkRange];
		double[] angleAtThisAngle = new double[maxAlpha/dAlpha];
		float[] xOD = new float[maxAlpha/dAlpha]; // x of outer diameter
		float[] yOD = new float[maxAlpha/dAlpha]; // y of outer diameter
		
		//init out vectors
		double[] xCenter = new double[samplingNumber];
		double[] yCenter = new double[samplingNumber];
		double[] zCenter = new double[samplingNumber];
		double[] majorRadius = new double[samplingNumber];
		double[] minorRadius = new double[samplingNumber];
		double[] medianOuterRadius = new double[samplingNumber];				
		double[] xyAngle = new double[samplingNumber];
		double[] R2 = new double[samplingNumber];
		double[] tilt = new double[3];	
	
								
		//loop over the some of the slices in the image..
		for (i = startSlice ; i < startSlice + samplingNumber * dSlice ; i = i + dSlice) {  
			
			IJ.showStatus("Searching for location and orientation of column..");

			int cc=0;
			nowTiff.setPosition(i);				
			ImageProcessor myIP = nowTiff.getProcessor();
			
			//determine steelGrey Threshold
			int[] myHist = myIP.getHistogram();
			double[] cumHist = hist.calcCumulativeHistogram(myHist);
			steelGreyThreshold = hist.findPercentileFromCumHist(cumHist, 0.95);
			
			//nowTiff.draw();
			//nowTiff.show();
			
			//sample grey values around the column
			for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {				
				
				int cc0=0;
				angleAtThisAngle[cc] = angle;
				for (double checkRadius = radius ; checkRadius > radius - checkRange ; checkRadius--) {					
					x[cc0] = Math.sin(angle) * checkRadius + xmid;
					y[cc0] = Math.cos(angle) * checkRadius + ymid;				
					greyAtThisAngle[cc0] = myIP.getPixelValue((int)x[cc0], (int)y[cc0]);					
					cc0++;
				}
				
				//apply medianFilter
				double[] mGreyAtThisAngle = math.oneDMedianFilter(greyAtThisAngle, footprintOfMedianFilter);
							
				//calculate the first derivative, dMG
				double[] dMG = math.oneDFirstDerivative(mGreyAtThisAngle, footprintOfMedianFilter);
								
				//find first steel
				int myWallPosition = 0;
				for (j = 0 ; j < mGreyAtThisAngle.length ; j++) if (mGreyAtThisAngle[j] > steelGreyThreshold) {
					myWallPosition = j;
					break;
				}
				
				//if first steel was not found..
				if (myWallPosition == 0) for (j = 0 ; j < mGreyAtThisAngle.length ; j++) if (mGreyAtThisAngle[j] > 0.95 * steelGreyThreshold) {
					myWallPosition = j;
					break;
				}
								
				//look for the outer edge
				double maxGradient = -50000;
				int maxGradPosition = 0;
				if (myWallPosition > gradientSearchWindow & myWallPosition < dMG.length - gradientSearchWindow) {
					for (j = myWallPosition - gradientSearchWindow ; j < myWallPosition + gradientSearchWindow; j++) {
						if (dMG[j] > maxGradient) {	
							maxGradient = dMG[j];
							maxGradPosition = j;
						}
					}
				} else {
					maxGradPosition = 0;
				}
				
				xOD[cc] = (float)x[maxGradPosition];
				yOD[cc] = (float)y[maxGradPosition];
				
				//increase vector index
				cc++;
			}
			
			//fit an elliptic ROI to the x and y
			EllipseFitter jEF = new EllipseFitter();			
			PolygonRoi pRoi = new PolygonRoi(xOD, yOD, Roi.POLYLINE); // create pRoi with outer Wall coordinates
			myIP.setRoi(pRoi); 
			jEF.fit(myIP, myIP.getStatistics());	
			
			//transfer fitting results to vectors for export
			xCenter[iCounter] = jEF.xCenter;
			yCenter[iCounter] = jEF.yCenter;
			zCenter[iCounter] = i;
			minorRadius[iCounter] = jEF.minor / 2;
			majorRadius[iCounter] = jEF.major / 2;
			medianOuterRadius[iCounter] = (minorRadius[iCounter] + majorRadius[iCounter]) / 2;
			xyAngle[iCounter] = jEF.theta;
			FitStuff.GoodnessOfFit gOF = fit.calculateR2OfEllipsoidFit(angleAtThisAngle, xOD, yOD, jEF);
			R2[iCounter] = gOF.R2;
			
			iCounter++;		
 
		}
		
		//find column Tilt!
		tilt = findColumnTilt(xCenter, yCenter, zCenter);
		
		//save results in struct
		pCC.xmid = xCenter;
		pCC.ymid = yCenter;
		pCC.zmid = zCenter;
		pCC.tiltInXZ = tilt[0];
		pCC.tiltInYZ = tilt[1];
		pCC.tiltTotal = tilt[2];
		
		return pCC;		
		
	}
	
	/*public ColumnCoordinates findAllSteelWalls(ImagePlus nowTiff, ColumnCoordinates prelimCC, double wallThickness) {

		//init objects
		Median jMed = new Median();
		ColumnCoordinates preciseCC = new ColumnCoordinates();
		TailoredMaths math = new TailoredMaths();
		
		//init outline related variables
		int i, j;
		int maxAlpha = 360;
		int dAlpha = 5;
		double xmid = nowTiff.getWidth() / 2;
		double ymid = nowTiff.getHeight() / 2;
		double radius = nowTiff.getWidth() / 2;
		int checkRange = (int)Math.round(radius / 2) ;		
		int footprintOfMedianFilter = 5;
		int imageHeight = nowTiff.getNSlices();
		double steelGreyThreshold = 25000;
		double appliedSteelGreyThreshold = steelGreyThreshold;
		int gradientSearchWindow = 5;
		
		//create angle vector
		double[] myAngle = new double[maxAlpha/dAlpha];
		int cc = 0;	
		for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) 
			{myAngle[cc] = angle; cc++;}
		
		//set starting point for outerColumn search
		double[] x = new double[checkRange];
		double[] y = new double[checkRange];	
		double[] greyAtThisAngle = new double[checkRange];
		double[] angleAtThisAngle = new double[maxAlpha/dAlpha];
		float[] xOD = new float[maxAlpha/dAlpha + 1]; // x of outer diameter
		float[] yOD = new float[maxAlpha/dAlpha + 1]; // y of outer diameter
		float[] xID = new float[maxAlpha/dAlpha + 1]; // x of inner diameter
		float[] yID = new float[maxAlpha/dAlpha + 1]; // y of inner diameter
		float[][] aXOD = new float[nowTiff.getNSlices()][maxAlpha/dAlpha + 1]; // x of outer diameter
		float[][] aYOD = new float[nowTiff.getNSlices()][maxAlpha/dAlpha + 1]; // y of outer diameter
		float[][] aXID = new float[nowTiff.getNSlices()][maxAlpha/dAlpha + 1]; // x of inner diameter
		float[][] aYID = new float[nowTiff.getNSlices()][maxAlpha/dAlpha + 1]; // y of inner diameter
		
		//init pre-out vectors
		double[] xCenter = new double[imageHeight];
		double[] yCenter = new double[imageHeight];
		double[] zCenter = new double[imageHeight];
		double[] tilt = new double[3];	
		
		//init out vectors
		PolygonRoi[] oRoi = new PolygonRoi[imageHeight]; 	
		PolygonRoi[] iRoi = new PolygonRoi[imageHeight]; 	
		
		//nowTiff.draw();
		//nowTiff.show();
								
		for (i = 0 ; i < nowTiff.getNSlices() ; i++) {  
		
			IJ.showStatus("Searching for column's inner wall of slice #" + (i + 1) + "/" + nowTiff.getNSlices());
							
			int coco=0;
			nowTiff.setPosition(i+1);				
			ImageProcessor myIP = nowTiff.getProcessor();
			
			//nowTiff.draw();
			//nowTiff.show();
			
			//sample grey values around the column
			for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {				
				
				int cc0=0;
				angleAtThisAngle[coco] = angle;
				for (double checkRadius = radius ; checkRadius > radius - checkRange ; checkRadius--) {					
					x[cc0] = Math.sin(angle) * checkRadius + xmid;
					y[cc0] = Math.cos(angle) * checkRadius + ymid;				
					greyAtThisAngle[cc0] = myIP.getPixelValue((int)x[cc0], (int)y[cc0]);					
					cc0++;
				}
				
				//apply medianFilter
				double[] mGreyAtThisAngle = math.oneDMedianFilter(greyAtThisAngle, footprintOfMedianFilter);
							
				//calculate the first derivative, dMG
				double[] dMG = math.oneDFirstDerivative(mGreyAtThisAngle, footprintOfMedianFilter);
		
				//find first steel
				int myWallPosition = 0;
				appliedSteelGreyThreshold = steelGreyThreshold;
				for (j = 5 ; j < mGreyAtThisAngle.length ; j++) if (mGreyAtThisAngle[j] > appliedSteelGreyThreshold) {
					myWallPosition = j;
					break;
				}
				
				//if first steel was not found..
				appliedSteelGreyThreshold = steelGreyThreshold * 3 / 4;
				if (myWallPosition == 0) for (j = 0 ; j < mGreyAtThisAngle.length ; j++) if (mGreyAtThisAngle[j] > appliedSteelGreyThreshold) {
					myWallPosition = j;
					break;
				}
						
				//if first steel was not found..
				appliedSteelGreyThreshold = steelGreyThreshold * 1 / 2;
				if (myWallPosition == 0) for (j = 0 ; j < mGreyAtThisAngle.length ; j++) if (mGreyAtThisAngle[j] > appliedSteelGreyThreshold) {
					myWallPosition = j;
					break;
				}
				
				//look for the outer edge
				double maxGradient = -500000;
				int maxGradPosition = 0;
				if (myWallPosition > gradientSearchWindow & myWallPosition < dMG.length - gradientSearchWindow) {
					for (j = myWallPosition - gradientSearchWindow ; j < myWallPosition + gradientSearchWindow; j++) {
						if (dMG[j] > maxGradient) {	
							maxGradient = dMG[j];
							maxGradPosition = j;
						}
					}
				} else {
					maxGradPosition = 0;
				}
				
				xOD[coco] = (float)x[maxGradPosition];
				yOD[coco] = (float)y[maxGradPosition];
				if (j + (int)Math.round(wallThickness) < x.length) {
					xID[coco] = (float)x[maxGradPosition + (int)Math.round(wallThickness)];
					yID[coco] = (float)y[maxGradPosition + (int)Math.round(wallThickness)];
				}
				
				//test approach
				//myIP.setValue(60000);
				//myIP.drawDot((int)xID[coco], (int)yID[coco]);
				//nowTiff.updateAndDraw();
				//nowTiff.show();
				
				//increase vector index
				coco++;
			}
			
			xOD[coco]=xOD[0];
			yOD[coco]=yOD[0];
			xID[coco]=xID[0];
			yID[coco]=yID[0];
			
			//fit an elliptic ROI to the x and y						
			oRoi[i] = new PolygonRoi(xOD, yOD, Roi.POLYLINE); // create pRoi with outer Wall coordinates
			iRoi[i] = new PolygonRoi(xID, yID, Roi.POLYLINE); // create pRoi with inner Wall coordinates
			
			//test if ROI was found properly..			 
			myIP.setValue(50000);	
			myIP.draw(oRoi[i]);
			myIP.draw(iRoi[i]);			
			nowTiff.updateAndDraw();
			nowTiff.show();
			
			double xODsum = 0; for(j = 0 ; j < xOD.length - 1 ; j++) xODsum += xOD[j];
			double yODsum = 0; for(j = 0 ; j < yOD.length - 1 ; j++) yODsum += yOD[j];
			xCenter[i] = xODsum / (xOD.length - 1);
			yCenter[i] = yODsum / (yOD.length - 1);
			zCenter[i] = i;
		
			//write ROI coordinates into array.. 
			for (j = 0 ; j < xOD.length ; j++) aXOD[i][j] = xOD[j]; 
			for (j = 0 ; j < yOD.length ; j++) aYOD[i][j] = yOD[j];
			for (j = 0 ; j < xID.length ; j++) aXID[i][j] = xID[j];
			for (j = 0 ; j < yID.length ; j++) aYID[i][j] = yID[j];
		 
		}
		
		//find column Tilt!
		tilt = findColumnTilt(xCenter, yCenter, zCenter);

		preciseCC.xmid = xCenter; 
		preciseCC.ymid = yCenter;
		preciseCC.zmid = zCenter;
		
		preciseCC.wallThickness = wallThickness;
		
	    //re-check the found values..
		j = 0;
		double[] myMinPos = new double[xOD.length - 1];
		for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {  				//loop over the angles..
			
			double[] thisRadius = new double[nowTiff.getNSlices()];	//vector for the inner radii..
			double[] thisOuter = new double[nowTiff.getNSlices()];	//vector for the outer radii..			
			
			for (i = 0 ; i < nowTiff.getNSlices() ; i++) {  	//loop over slices
				
				//check inner radius
				double x0 = xCenter[i];
				double y0 = yCenter[i];
				double x1 = aXID[i][j]; 
				double y1 = aYID[i][j];
				double dx = x1 - x0;
				double dy = y1 - y0;
				
				thisRadius[i] = Math.sqrt(dx * dx + dy * dy);				
				
				//also check outerRadius
				double ox1 = aXOD[i][j]; 
				double oy1 = aYOD[i][j];
				double ox = ox1 - x0;
				double oy = oy1 - y0;
				
				thisOuter[i] = Math.sqrt(ox * ox + oy * oy);
				
			}
			
			//find the minimum radii of the bottom of the column.. 
			double[] bottomRadii = new double[thisOuter.length - 500];
			for (int k = 500 ; k < thisOuter.length ; k++) bottomRadii[k-500] = thisOuter[k];
			double myMin = StatUtils.min(bottomRadii);			
			for (int k = 0 ; k < thisOuter.length ; k++) if (thisOuter[k] == myMin) {
				myMinPos[j] = k;			//remember where the minimum radius was .. == bottom of column..
				break;
			}
			
			//handle column's bottom..
			double[] middleRadii = new double[3];
			int sampleLastInnerDiameter = (int)myMinPos[j] - 3 - (int)preciseCC.coningStartsBeforeEnd;
			for (int k = sampleLastInnerDiameter; k < sampleLastInnerDiameter + 3; k++) middleRadii[k - sampleLastInnerDiameter] = thisRadius[k];
			double mRadius = jMed.evaluate(middleRadii);
			double aRadius = jMed.evaluate(thisRadius);
			
			//correct inner wall coordinates where necessary
			for (i = 0 ; i < nowTiff.getNSlices() ; i++) {
				
				if (Math.sqrt(Math.pow(thisRadius[i] - mRadius,2)) > wallThickness / 2) {
				
					float xnew = (float)xCenter[i] + (float)(Math.sin(angle) * aRadius);
					float ynew = (float)yCenter[i] + (float)(Math.cos(angle) * aRadius);					
					//float xnow = aXID[i][j]; 
					//float ynow = aYID[i][j];
					
					aXID[i][j] = xnew;
					aYID[i][j] = ynew;
					
				}
				
				if (i + preciseCC.coningStartsBeforeEnd > myMinPos[j]) {
					
					float xnew = (float)xCenter[i] + (float)(Math.sin(angle) * mRadius);
					float ynew = (float)yCenter[i] + (float)(Math.cos(angle) * mRadius);					
					//float xnow = aXID[i][j]; 
					//float ynow = aYID[i][j];
					
					aXID[i][j] = xnew;
					aYID[i][j] = ynew;
					
				}
				
				//add the additional point to close the circle..
				if (j == 0) {
					aXID[i][xID.length - 1] = aXID[i][j];
					aYID[i][yID.length - 1] = aYID[i][j];
				}
								
			}	
			
			j++; 			//set angle counter + 1
		}
			
		preciseCC.xOD = aXOD;
		preciseCC.yOD = aYOD;
		preciseCC.xID = aXID;
		preciseCC.yID = aYID;		
		
		preciseCC.bottomOfColumn = (int)jMed.evaluate(myMinPos) - 5;		
		
		preciseCC.anglesChecked = xOD.length - 1;

		preciseCC.tiltInXZ = tilt[0];
		preciseCC.tiltInYZ = tilt[1];
		preciseCC.tiltTotal = tilt[2];
						
		for(i = 1 ; i < nowTiff.getStackSize() ; i++) {
			
			Overlay mO = new Overlay();
			mO.add(oRoi[i]);
			
			nowTiff.setPosition(i);				
			ImageProcessor myIP = nowTiff.getProcessor();
			
			myIP.setRoi(oRoi[i]);
			myIP.drawOverlay(mO);
			
			mO.remove(oRoi[i]);
			
		}
		
		return preciseCC;		
		
	}*/
	
	public int[] findTopAndBottomOfPVCColumn(ImagePlus nowTiff, ColumnCoordinates pCC) {
		
		//init objects		
		Median jMed = new Median();
		HistogramStuff hist = new HistogramStuff();

		//init outline related variables		
		int i;
		int[] topAndBot = new int[2];
		double imageHeight = nowTiff.getNSlices();
		double[] myP = new double[(int)imageHeight];
		double[] dP = new double[(int)imageHeight - 1];
		double threshold = 0.99;

		for (i = 0 ; i < imageHeight ; i++) {
		
			IJ.showStatus("Searching slice for column's top and bottom #" + (i + 1) + "/" + imageHeight);
		
			nowTiff.setPosition(i+1);				
			ImageProcessor myIP = nowTiff.getProcessor();
			
			int[] myHist = myIP.getHistogram();						
			double[] cumHist = hist.calcCumulativeHistogram(myHist);
		
			//calculate 99 percentile of grey values at this depth
			myP[i] = hist.findPercentileFromCumHist(cumHist, threshold);
			
			//calculate 1st derivative of 99 percentile profile
			if (i > 0) dP[i-1] = myP[i] - myP[i-1];
			
		}
		
		//search maximum and minimum derivatives
		double minDP = 0;
		double maxDP = 0;
		
		for (i = 0 ; i < imageHeight - 1 ; i++) {
			if (dP[i] < minDP & dP[i] < -0.1 * jMed.evaluate(myP)) {
				minDP = dP[i];
				topAndBot[0] = i;
			}
			if (dP[i] > maxDP & dP[i] > 0.1 * jMed.evaluate(myP)) {
				maxDP = dP[i];
				topAndBot[1] = i;
			}
			
		}
		
		return topAndBot;
		
	}
	
	public int findMedianSoilSurfacePosition(ImagePlus surfTiff) {
		
		HistogramStuff hist = new HistogramStuff();
		
		surfTiff.setPosition(1);
		ImageProcessor topSoilIP = surfTiff.getProcessor();
		int[] topSurfHist = topSoilIP.getHistogram();
		topSurfHist[0] = 0;
		return hist.findMedianFromHistogram(topSurfHist);		
	
	}
	
	/*public ColumnCoordinates findTopAndBottomOfSteelColumn(ImagePlus nowTiff, ColumnCoordinates pCC) {
		
		//init objects
		Median jMed = new Median();
		HistogramStuff hist = new HistogramStuff();
				
		//init outline related variables
		int i;
		double imageHeight = nowTiff.getNSlices();
		double[] myP = new double[(int)imageHeight];
		double[] dP = new double[(int)imageHeight - 1];
		double threshold = 0.99;

		for (i = 0 ; i < imageHeight ; i++) {
		
			IJ.showStatus("Searching slice for column's top and bottom #" + (i + 1) + "/" + imageHeight);
		
			nowTiff.setPosition(i+1);				
			ImageProcessor myIP = nowTiff.getProcessor();
			
			int[] myHist = myIP.getHistogram();						
			double[] cumHist = hist.calcCumulativeHistogram(myHist);
					
			//calculate 99 percentile of grey values at this depth
			myP[i] = hist.findPercentileFromCumHist(cumHist, threshold);
			
			//calculate 1st derivative of 99 percentile profile
			if (i > 0) dP[i-1] = myP[i] - myP[i-1];
			
		}
		
		//search maximum and minimum derivatives
		double minDP = 0;
		double maxDP = 0;
		pCC.bottomOfColumn = myP.length; 
		pCC.topOfColumn = 0;
		
		for (i = 0 ; i < imageHeight - 1 ; i++) {
			if (dP[i] < minDP & myP[i] < -0.1 * jMed.evaluate(myP)) {
				minDP = dP[i];
				pCC.bottomOfColumn = i;
				
				if (pCC.bottomOfColumn - pCC.topOfColumn < 0.1 * pCC.approximateHeight) {
					maxDP = 0;
					minDP = 0;
				}
				
			}
			if (dP[i] > maxDP & myP[i] > 0.1 * jMed.evaluate(myP)) {
				maxDP = dP[i];
				pCC.topOfColumn = i;
			}
			
		}
		
		return pCC;
		
	}*/
	
	public ImagePlus extractFreshOM(ImagePlus nowTiff, MenuWaiter.OMFinderSettings oMF) {
		
		ImageCalculator myIC =  new ImageCalculator();
		
		ImagePlus outTiff = new ImagePlus();	

		int i, j;
		int numberOfWindows = (oMF.maxGreyValue - oMF.minGreyValue - oMF.overlap) / (oMF.windowSize - oMF.overlap);
		int[] lowerWindowBound = new int[numberOfWindows]; 
		for (i = 0 ; i < numberOfWindows ; i++) lowerWindowBound[i] = oMF.minGreyValue + i * (oMF.windowSize - oMF.overlap);
		
		//init sub-images to merge
		ImagePlus[] mTiff = new ImagePlus[numberOfWindows];
		
		for (j = 0 ; j < numberOfWindows ; j++) {
			
			ImageStack outStack1 = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
			ImageStack outStack2 = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
			
			for (i = 0 ; i < nowTiff.getNSlices() ; i++) {			
					
				IJ.showStatus("Searching for fresh organic material stage " + (j + 1) + "/" + numberOfWindows  + " in slice " + (i + 1) + "/" + nowTiff.getNSlices());
		
				//set tiff to the correct position and get Processor etc..
				nowTiff.setPosition(i+1);
				ImageProcessor myIP = nowTiff.getProcessor();
				ByteProcessor modIP1 = myIP.duplicate().convertToByteProcessor(false);
				ByteProcessor modIP2 = myIP.duplicate().convertToByteProcessor(false);
											
				//apply threshold	
				modIP1.threshold(lowerWindowBound[j]);
				//modIP2.invert();
				modIP2.threshold(lowerWindowBound[j] + oMF.windowSize);
				//modIP2.invert();
				
				outStack1.addSlice(modIP1);
				outStack2.addSlice(modIP2);
			}
					
			ImagePlus zTiff1 = new ImagePlus();
			zTiff1.setStack(outStack1);
			ImagePlus zTiff2 = new ImagePlus();
			zTiff2.setStack(outStack2);
			
			//merge the two images			
			outTiff = myIC.run("xor create stack", zTiff1, zTiff2);
			zTiff1.killStack();
			zTiff2.killStack();
			System.gc();
			
			//erode 3D
			Erode_ jE = new Erode_();
			ImagePlus outTiff2 = jE.erode(outTiff, 255, true);
			outTiff.killStack();
			System.gc();
			
			//kick out clusters smaller than 1000
			JParticleCounter jJPA = new JParticleCounter();
			int slicesPerChunk = 2; //input parameter 4 particle analyzer
			double minVol = 1000; //
			double maxVol = Double.POSITIVE_INFINITY; //input parameter 4 particle analyzer
			Boolean doExclude = false;
			int FORE = -1;
			
			Object[] result = jJPA.getParticles(outTiff2, slicesPerChunk, minVol, maxVol, FORE, doExclude);
			int[][] particleLabels = (int[][]) result[1];
			outTiff = jJPA.displayParticleLabels(particleLabels, outTiff2);			
			
			//re-convert to binary image
			ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
			for (i = 0 ; i < nowTiff.getNSlices() ; i++) {			
				
				//set tiff to the correct position and get Processor etc..			
				outTiff.setPosition(i+1);
				ImageProcessor myIP = outTiff.getProcessor().convertToByteProcessor(false);
				ImageProcessor modIP = myIP.duplicate();
				modIP.threshold(0);
				
				outStack.addSlice(modIP);
			}
			
			outTiff2.setStack(outStack);				
			outTiff.killStack();
			System.gc();
			
			//dilate 3D
			Dilate_ jD = new Dilate_();
			mTiff[j] = jD.dilate(outTiff2, 255, true);
			outTiff2.killStack();
			System.gc();
		}
		
		//merge the sub-images
		ImagePlus zwischi = new ImagePlus();		
		if (numberOfWindows > 1) zwischi = myIC.run("or create stack",mTiff[0], mTiff[1]);
		else outTiff = mTiff[0];
		if (numberOfWindows > 2) {
			for (i = 2 ; i < numberOfWindows ; i++) {
				zwischi = myIC.run("or create stack",zwischi, mTiff[i]);
			}
		}
		
		outTiff = zwischi;
		
		return outTiff;
	}
	
	public BestWallFinds findBestWallDetectionPerSector(int numOfSectors, int i, float[] xOD, float[] yOD, double[] myAngle, ColumnCoordinates prelimCC) {
				
		BestWallFinds bWF = new BestWallFinds();
		
		int numOfAngles = myAngle.length;
		int increment = numOfAngles / numOfSectors;
		
		double xMid = StatUtils.percentile(prelimCC.xmid, 50);
		double yMid = StatUtils.percentile(prelimCC.ymid, 50);
		
		ArrayList<Integer> bestOnes = new ArrayList<Integer>();
		
		for (int j = 0 ; j < numOfAngles ; j += increment) {
			
			int bestOne = -1;
			double farthestDist = -1;
			
			for (int k = 0 ; k < increment ; k++) {
			
				double xDist = xOD[j + k] - xMid;
				double yDist = yOD[j + k] - yMid;				
				double nowDist = Math.sqrt(xDist * xDist + yDist * yDist);
				
				if (nowDist > farthestDist) {
					
					farthestDist = nowDist;
					bestOne = j + k;
					
				}		
			}
			
			if (bestOne > 1) {
				bestOnes.add(bestOne);
			}
			
		}
		
		float[] xD = new float[bestOnes.size()];
		float[] yD = new float[bestOnes.size()];
		double[] angle = new double[bestOnes.size()];
		
		for (int j = 0 ; j < bestOnes.size() ; j++) {
			
			xD[j] = xOD[bestOnes.get(j)];
			yD[j] = yOD[bestOnes.get(j)];
			angle[j] = myAngle[bestOnes.get(j)];
			
		}
		
		bWF.xD = xD;
		bWF.yD = yD;
		bWF.angle = angle;
		
		return bWF;
	}
	
	public ColumnContainer findColumnsWalls(ColumnContainer colCon, String countMsg) {
				
		//import objects
		Median jMed = new Median();
		TailoredMaths maths = new TailoredMaths();
		FitStuff fit = new FitStuff();		
		RoiHandler roi = new RoiHandler();
				
		//init variables
		int maxAlpha = 360;
		double[] myAngle = new double[maxAlpha/colCon.dAlpha];		
		int cc = 0;	
		for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/colCon.dAlpha)) 
			{myAngle[cc] = angle; cc++;}
		
		//set starting point for outerColumn search		
		int footprintOfMedianFilter = colCon.jCFS.medianFilter;
		if ((Math.round((float)footprintOfMedianFilter/2) - Math.floor((float)footprintOfMedianFilter/2)) == 0) footprintOfMedianFilter += 1;

		double xmid = jMed.evaluate(colCon.prelimCC.xmid);
		double ymid = jMed.evaluate(colCon.prelimCC.ymid);
		double radius = (jMed.evaluate(colCon.prelimCC.outerMajorRadius) + jMed.evaluate(colCon.prelimCC.outerMinorRadius)) / 2;
		double startingPoint = Math.round(radius) + 30;
		double stoppingPoint = Math.round(radius) - 200;	
		int myWindowSize = (int)Math.ceil(0.003 * radius);
		
		//allow start, stop and stepsize in case of debug mode
		if (colCon.jCFS.debug) {
			//show precision of wall findings
			GenericDialog gd = new GenericDialog("Slice Stepper");
			
			gd.addNumericField("start slice", colCon.startslice, 0);
			gd.addNumericField("stop slice", colCon.stopslice, 0);
			gd.addNumericField("increment", 100, 0);
			
			gd.showDialog();
		    if (gd.wasCanceled()) return null;
		    else {		    	
		    	colCon.startslice = (int)gd.getNextNumber();
		    	colCon.stopslice = (int)gd.getNextNumber();
		    	colCon.increment = (int)gd.getNextNumber();	
		    }
		    
		    //probe illumination levels of image
		    //ApproximateColumnIllumination aCI = probeApproximateColumnIllumination(colCon);
		    //MenuWaiter.ColumnFinderMenuReturn mW = menu.showTuneThreshold4ColumnDetectionMenu(aCI);		    
		}
	
		//init column wall coordinates
		int colHeight = colCon.nowTiff.getNSlices();
		double[] xCenter = new double[colHeight];
		double[] yCenter = new double[colHeight];
		double[] zCenter = new double[colHeight];
		double[] outerMinorRadius = new double[colHeight];
		double[] outerMajorRadius = new double[colHeight];				
		double[] theta = new double[colHeight];		
		double[] outerR2 = new double[colHeight];
		double[] ixCenter = new double[colHeight];
		double[] iyCenter = new double[colHeight];			
		double[] innerMinorRadius = new double[colHeight];
		double[] innerMajorRadius = new double[colHeight];				
		double[] itheta = new double[colHeight];		
		double[] innerR2 = new double[colHeight];		
		boolean[] columnIsAtThisDepth = new boolean[colHeight];		
		double[] wallThickness = new double[colHeight];		
				
		//re-find the column's outer wall		
		for (int i = colCon.startslice ; i < colCon.stopslice ; i+=colCon.increment) {  
		
			IJ.showStatus(countMsg + (i + 1) + "/" + colCon.stopslice);
			
			//init vectors
			float[] xOD = new float[maxAlpha/colCon.dAlpha]; // x of outer diameter
			float[] yOD = new float[maxAlpha/colCon.dAlpha]; // y of outer diameter
			double[] locationOD = new double[maxAlpha/colCon.dAlpha]; // y of outer diameter
			float[] xID = new float[maxAlpha/colCon.dAlpha]; // x of inner diameter
			float[] yID = new float[maxAlpha/colCon.dAlpha]; // y of inner diameter
			double[] locationID = new double[maxAlpha/colCon.dAlpha]; // y of inner diameter
		
			//set tiff to the correct position and get Processor etc..
			colCon.nowTiff.setPosition(i+1);
			ImageProcessor myIP = colCon.nowTiff.getProcessor();			
			
			//find illumination level 			
			HistogramStuff.IlluminationInfo jII = probeIlluminationLevel(myIP, colCon, true);
			int myContrast = jII.quantile99 - jII.lowerQuartile;
			int minAirColGradient = (int)Math.round(colCon.jCFS.airWallContrast * (double)myContrast);			
			int minColSoilGradient = (int)Math.round(colCon.jCFS.wallSoilContrast * ((double)jII.quantile99 - (double)jII.lowerQuartile));
			double wallBrightnessRange = 4 * colCon.jCFS.coeffVarOfWallBrightness * (double)jII.quantile99;
			int minWallBrightness = (int)Math.round((double)(1 - 4 * colCon.jCFS.coeffVarOfWallBrightness) * (double)jII.quantile99);
			
			
			double[][] x = new double[myAngle.length][(int)Math.round(startingPoint + 1)];
			double[][] y = new double[myAngle.length][(int)Math.round(startingPoint + 1)];
			
			boolean[] innerPerimeterFound = new boolean[myAngle.length];
			boolean[] outerPerimeterFound = new boolean[myAngle.length];
			
			for (int j = 0 ; j < myAngle.length ; j++) {
	
				double[] greyAtThisAngle = new double[(int)Math.round(startingPoint + 1)];
				locationOD[j] = 0; //nil it!
				
				//re-find outer column wall
				cc=0;
				for (double checkRadius = startingPoint ; checkRadius > stoppingPoint ; checkRadius--) {					
					x[j][cc] = Math.sin(myAngle[j]) * checkRadius + xmid;
					y[j][cc] = Math.cos(myAngle[j]) * checkRadius + ymid;				
					greyAtThisAngle[cc] = myIP.getPixelValue((int)x[j][cc], (int)y[j][cc]);					
					cc++;
				}
				
				//apply medianFilter
				double[] mGreyAtThisAngle = maths.oneDMedianFilter(greyAtThisAngle, footprintOfMedianFilter);
				
				//calculate the first derivative, dMG
				double[] dMG = maths.oneDFirstDerivative(mGreyAtThisAngle, footprintOfMedianFilter);
				
				//set wall found trackers to false
				innerPerimeterFound[j] = false;
				outerPerimeterFound[j] = false;
				
				//look for the outer and then inner edge				
				for (int k = 0 ; k < (int)Math.round(startingPoint) - myWindowSize - 2; k++) {
					double[] theNextInTheWindow = new double[myWindowSize];
					for(int l = 0 ; l < theNextInTheWindow.length ; l++) theNextInTheWindow[l] = mGreyAtThisAngle[k+l+2];
					if (dMG[k] > minAirColGradient & StatUtils.min(theNextInTheWindow) >= minWallBrightness & locationOD[j] == 0) {
						int myJ = k;
						for (int l = 1 ; l < 3 ; l++) {		//find last gradient over threshold..
							if (dMG[k+l] > minAirColGradient) myJ = k + l; 						
						}
						xOD[j] = (float)x[j][myJ];
						yOD[j] = (float)y[j][myJ];
						locationOD[j] = myJ;
						outerPerimeterFound[j] = true;
						if (!colCon.jCFS.isAlu) break;
					}
					if (k >= myWindowSize + 1 & locationOD[j] > 0 & locationID[j] == 0) {
						double[] thePervInTheWindow = new double[myWindowSize];
						for(int l = 0 ; l < thePervInTheWindow.length ; l++) thePervInTheWindow[l] = mGreyAtThisAngle[k-l];
						double[] twoInARowForTheWin = {greyAtThisAngle[k+1], greyAtThisAngle[k+2]};						
						if (dMG[k] < -minColSoilGradient & 
								StatUtils.min(thePervInTheWindow) > minWallBrightness & 
								(StatUtils.mean(twoInARowForTheWin) > minWallBrightness + 2 * wallBrightnessRange | 
								StatUtils.mean(twoInARowForTheWin) < minWallBrightness) &
								!innerPerimeterFound[j]) {
							int myJ = k;
							int myWindow = k + 3;
							if ((int)Math.round(startingPoint + 1) <= myWindow) myWindow = (int)Math.round(startingPoint + 1); 
							for (int l = k ; l < myWindow ; l++) {
								if (dMG[l] < dMG[myJ]) myJ = l;
							}		
							xID[j] = (float)x[j][myJ];
							yID[j] = (float)y[j][myJ];
							locationID[j] = myJ;
							innerPerimeterFound[j] = true; 
							break;
						}
					}	
				}
				
			/*	//display radial profile
				DisplayThings dT = new DisplayThings();
				double[] X = new double[dMG.length];
				for (int k = 0 ; k < X.length ; k++) X[k] = (double)k;
				double[] foundPoints = new double[2];
				foundPoints[0] = locationID[j];foundPoints[1] = locationOD[j];
				double[] plotAtZero = {0d, 0d};
				
				dT.plotXYXY(X, dMG, foundPoints, plotAtZero, "angle #" + j, "radial distance (vx)", "dMG");			*/	
				
			}
			
			//filter out the best wall detection per pi/4 sector for outer perimeter
			//BestWallFinds bWF = findBestWallDetectionPerSector(numOfSectors, i, xOD, yOD, myAngle, colCon.prelimCC); 			
						
			//fit an elliptic ROI to the x and y of outer edge
			FitStuff.FittedEllipse fO = fit.doRobustEllipseFit(i, xOD, yOD, myAngle, myIP, colCon.jCFS);
								
			//transfer fitting results to vectors for export
			xCenter[i] = fO.xCenter;
			yCenter[i] = fO.yCenter;
			zCenter[i] = i;
			outerMinorRadius[i] = fO.minorRadius;
			outerMajorRadius[i] = fO.majorRadius;				
			theta[i] = fO.theta;		
			outerR2[i] = fO.R2;
			
			//first assume it is a PVC column and infer to inner diameter by subtracting wall thicknes
			FitStuff.FittedEllipse fI = fit.doRobustEllipseFit(i, xOD, yOD, myAngle, myIP, colCon.jCFS);
			fI.minorRadius -= colCon.jCFS.wallThickness;	
			fI.majorRadius -= colCon.jCFS.wallThickness;				

			//filter out the best wall detection per pi/4 sector for inner perimeter
			//BestWallFinds ibWF = findBestWallDetectionPerSector(numOfSectors, i, xID, yID, myAngle, colCon.prelimCC); 
			
			//in case it is an alu column, fit fI								
			if (colCon.jCFS.isAlu) fI = fit.doRobustEllipseFit(i, xID, yID, myAngle, myIP, colCon.jCFS);
			
			//transfer results to array
			ixCenter[i] = fI.xCenter;
			iyCenter[i] = fI.yCenter;			
			innerMinorRadius[i] = fI.minorRadius;
			innerMajorRadius[i] = fI.majorRadius;				
			itheta[i] = fI.theta;		
			innerR2[i] = fI.R2;	
			
			//calculate wall thickness...			
			wallThickness[i] = 0.5 * (outerMajorRadius[i] - innerMajorRadius[i] + outerMinorRadius[i] - innerMinorRadius[i]) ;
			
			//compare new fits with formerly found properties and decide whether the column was at this depth or not..			
			boolean colDetected = true; //first set to true			
			double czechX_IX2 = Math.sqrt(Math.pow(xCenter[i] - ixCenter[i],2));
			double czechY_IY2 = Math.sqrt(Math.pow(yCenter[i] - iyCenter[i],2));
			if (xCenter[i] == 0 | ixCenter[i] == 0 | yCenter[i] == 0 | iyCenter[i] == 0) colDetected = false;
			if (Double.isNaN(outerR2[i])) colDetected = false;
			if (outerR2[i] < colCon.jCFS.minR24PerimeterFit) colDetected = false;	
			if (colCon.jCFS.isAlu & innerR2[i] < colCon.jCFS.minR24PerimeterFit) colDetected = false;	
			if (calculateCVOfWallThickness(fI, fO) > colCon.jCFS.maxCVOfWallThickness) colDetected = false;	
			if (myContrast < colCon.jCFS.minAbsoluteGradientbetweenAirAndWall) colDetected = false;	
			if (wallThickness[i] < myWindowSize + 3) colDetected = false;
			columnIsAtThisDepth[i] = colDetected;
			
			if (colCon.jCFS.debug) {
				
				//add Overlays
				PolygonRoi outerPerimeter = roi.makeRoiFromFittedEllipse(fO);
				PolygonRoi innerPerimeter = roi.makeRoiFromFittedEllipse(fI);

				Overlay myO = new Overlay(outerPerimeter);
				myO.add(innerPerimeter);
				myO.setStrokeColor(Color.YELLOW);
				colCon.nowTiff.setOverlay(myO);
				
				//also adjust image brightness and contrast		
				ContrastEnhancer myCE = new ContrastEnhancer();
				myCE.stretchHistogram(colCon.nowTiff, 0.35);
				colCon.nowTiff.updateAndDraw();
				colCon.nowTiff.show();		
				
				//show precision of wall findings
				GenericDialog gd = new GenericDialog("ColRecogDebugger");
				
				gd.addMessage("czech X: " + String.format("%3.6f\t",(float)czechX_IX2));
				gd.addMessage("czech Y: " + String.format("%3.6f\t",(float)czechY_IY2));
				if (!colCon.jCFS.isPVC) gd.addMessage("inner R2: " + String.format("%1.6f\t",(float)innerR2[i]));
				gd.addMessage("outer R2: " + String.format("%1.6f\t",(float)outerR2[i]));
				gd.addMessage("wall thickness: " + String.format("%3.1f\t",(float)wallThickness[i]));		
				gd.addMessage("CV of wall thickness: " + String.format("%3.2f\t",(float)calculateCVOfWallThickness(fI, fO)));
				gd.addMessage("column found: " + columnIsAtThisDepth[i]);
				
				gd.showDialog();
			    if (gd.wasCanceled()) return null;
			    else new WaitForUserDialog("press any key");
				
			}	
			
			int u = 0;
		}
		
		//assign results to colCon
		colCon.preciseCC.xmid = xCenter;		
		colCon.preciseCC.ymid = yCenter;
		colCon.preciseCC.zmid = zCenter;
		colCon.preciseCC.outerMinorRadius = outerMinorRadius;
		colCon.preciseCC.outerMajorRadius = outerMajorRadius;				
		colCon.preciseCC.theta = theta;		
		colCon.preciseCC.outerR2 = outerR2;
		colCon.preciseCC.ixmid = ixCenter;
		colCon.preciseCC.iymid = iyCenter;			
		colCon.preciseCC.innerMinorRadius = innerMinorRadius;
		colCon.preciseCC.innerMajorRadius = innerMajorRadius;				
		colCon.preciseCC.itheta = itheta;
		colCon.preciseCC.innerR2 = innerR2;
		colCon.preciseCC.columnIsInThisDepth = columnIsAtThisDepth;
		colCon.preciseCC.wallThickness = wallThickness;
				
		return colCon;
		
	}
		
	public ColumnContainer findPVCOrAluWalls(ImagePlus nowTiff, ColumnCoordinates prelimCC, MenuWaiter.ColumnFinderMenuReturn jCFS) {
		
		//init objects	
		ColumnContainer colCon = new ColumnContainer();
		
		//setup colCon
		colCon.dAlpha = 5;			
		colCon.startslice = 0;
		colCon.stopslice = nowTiff.getNSlices();
		colCon.increment = 1;	
		colCon.prelimCC = prelimCC;		
		colCon.preciseCC = prelimCC;
		colCon.nowTiff = nowTiff;
		colCon.jCFS = jCFS;		
		
		//find the column wall
		colCon = findColumnsWalls(colCon, "Re-finding column wall in slice #");
		
		//flag doubtworthy results
		colCon = flagDoubtworthyResults(colCon);
		
		//look if the detected top and bottom of the columns are likely correct		
		colCon = findtopAndBottomOfSoilColumn(colCon);
		
		//find bevel of column
		colCon = findColumnsBevel(colCon);
		
		//impute layer outline on top if R2 < chosen confidence criterion
		colCon = imputeMissingLayers(colCon);
		
		//detect rubber band
		//colCon = rubberBandDetection(colCon);
		
		return colCon;
		
	}
	
	public ColumnContainer findtopAndBottomOfSoilColumn(ColumnContainer colCon) {
		
		colCon.preciseCC.topOfColumn = 0;		
		colCon.preciseCC.bottomOfColumn = 0;
		if (colCon.jCFS.try2FindColumnTopAndBottom) {
			for (int i = 0 ; i < colCon.nowTiff.getNSlices() ; i++) {
				if ((colCon.preciseCC.topOfColumn == 0) && colCon.preciseCC.columnIsInThisDepth[i]) colCon.preciseCC.topOfColumn = i;
				if ((colCon.preciseCC.topOfColumn > 0) && colCon.preciseCC.columnIsInThisDepth[i]) colCon.preciseCC.bottomOfColumn = i;
			}
		}
		else {
			colCon.preciseCC.topOfColumn =	colCon.jCFS.topOfColumn; 
			if (colCon.jCFS.bottomOfColumn == 1) colCon.preciseCC.bottomOfColumn = colCon.nowTiff.getNSlices();
			else colCon.preciseCC.bottomOfColumn = colCon.jCFS.bottomOfColumn;
		}		
		colCon.preciseCC.heightOfColumn = colCon.preciseCC.bottomOfColumn - colCon.preciseCC.topOfColumn;
		
		return colCon;
	}
	
	public ColumnContainer flagDoubtworthyResults(ColumnContainer colCon) {
		
		//set R2s for columnNotInThisDepth to 0;
		for (int z = colCon.preciseCC.topOfColumn ; z < colCon.preciseCC.bottomOfColumn ; z++) {
			if (!colCon.preciseCC.columnIsInThisDepth[z]) {
				colCon.preciseCC.innerR2[z] = 0;
				colCon.preciseCC.outerR2[z] = 0;
			}
			if (colCon.preciseCC.outerR2[z] == 0) colCon.preciseCC.columnIsInThisDepth[z] = false;
		}
		
/*		//find layers with complete values
		ArrayList<Integer> hasOuterPerimeter = new ArrayList<Integer>();
		ArrayList<Integer> hasInnerPerimeter = new ArrayList<Integer>();
		for (int z = colCon.preciseCC.topOfColumn ; z < colCon.preciseCC.bottomOfColumn ; z++) {			
			if (colCon.preciseCC.innerR2[z] >= colCon.jCFS.minR24PerimeterFit) hasInnerPerimeter.add(z);
			if (colCon.preciseCC.outerR2[z] >= colCon.jCFS.minR24PerimeterFit) hasOuterPerimeter.add(z);
		}
		
		double[] xmid = new double[hasOuterPerimeter.size()];
		double[] ymid = new double[hasOuterPerimeter.size()];
		double[] ixmid = new double[hasInnerPerimeter.size()];*/
		
		return colCon;
	}
	
	public ColumnContainer findColumnsBevel(ColumnContainer colCon) {
		
		FitStuff fit = new FitStuff();
		
		//find entries with outer perimeter data 
		ArrayList<Integer> outerPerimeter = new ArrayList<Integer>();
		for (int i = colCon.preciseCC.topOfColumn ; i < colCon.preciseCC.bottomOfColumn ; i++) {
			if (colCon.preciseCC.xmid[i] > 0 & colCon.preciseCC.ymid[i] > 0 & colCon.preciseCC.ixmid[i] > 0 & colCon.preciseCC.iymid[i] > 0) {
				outerPerimeter.add(i);
			}
		}
		
		//define the median wall thickness above the bevel..
		ArrayList<Double> wallThicknessWithoutBevel = new ArrayList<Double>();
		for (int i = 0 ; i < outerPerimeter.size() ; i++) {
			int j = outerPerimeter.get(i);
			if (colCon.preciseCC.zmid[j] > colCon.preciseCC.topOfColumn + 0.5 * colCon.preciseCC.heightOfColumn) {
				wallThicknessWithoutBevel.add(colCon.preciseCC.wallThickness[j]);
			}
		}
		double[] wallThickness = new double[wallThicknessWithoutBevel.size()]; 
		for (int i = 0 ; i < wallThickness.length ; i++) wallThickness[i] = wallThicknessWithoutBevel.get(i);
		double medianWallThickness = StatUtils.percentile(wallThickness, 50);
		
		//find all parts of the column with wallthickness smaller than 0.66 * median WallThickness
		ArrayList<Double> bevelWall = new ArrayList<Double>();
		ArrayList<Double> z = new ArrayList<Double>();
		for (int i = 0 ; i < outerPerimeter.size() ; i++) {
			int j = outerPerimeter.get(i);
			if (colCon.preciseCC.wallThickness[j] < 0.9 * medianWallThickness) {
				bevelWall.add(colCon.preciseCC.wallThickness[j]);
				z.add(colCon.preciseCC.zmid[j]);
			}
		}
		
		//fit straight to bevel 
		double[] Z = new double[z.size()];for (int i = 0 ; i < Z.length ; i++) Z[i] = z.get(i);
		double[] bw = new double[bevelWall.size()];for (int i = 0 ; i < bw.length ; i++) bw[i] = bevelWall.get(i);
		colCon.preciseCC.bevelFunction = fit.fitLinearFunction(Z, bw);
		
		//evaluate start of bevel
		double nowColumnWall = bevelWall.get(0);
		double zz = z.get(0);
		while (nowColumnWall < medianWallThickness ) {
			zz--;
			nowColumnWall = fit.evalLinearFunction(colCon.preciseCC.bevelFunction, zz);
		}
		colCon.preciseCC.bevelStartsAt = (int)Math.round(zz++);
		
		return colCon;
	}
	
	public ColumnContainer imputeMissingLayers(ColumnContainer colCon) {
		
		FitStuff fit = new FitStuff();
		
		//init mendedCC
		colCon.mendedCC = colCon.preciseCC;
		
		//find layers with complete values
		ArrayList<Integer> hasOuterPerimeter = new ArrayList<Integer>();
		ArrayList<Integer> hasInnerPerimeter = new ArrayList<Integer>();
		for (int z = colCon.preciseCC.topOfColumn ; z < colCon.preciseCC.bottomOfColumn ; z++) {
			if (colCon.preciseCC.outerR2[z] >= colCon.jCFS.minR24PerimeterFit) hasOuterPerimeter.add(z);
			if (colCon.preciseCC.innerR2[z] >= colCon.jCFS.minR24PerimeterFit) hasInnerPerimeter.add(z);
		}
		
		//impute the missing values for inner perimeter...
		for (int z = colCon.preciseCC.topOfColumn ; z < colCon.preciseCC.bottomOfColumn ; z++) {
			
			if (!hasInnerPerimeter.contains(z)) {			
			
				colCon.mendedCC.numberOfImputedLayers++;
				
				//find existing entries above and below the missing value
				int data4Imputation = -1;
				for (int j = 0 ; j < hasInnerPerimeter.size() ; j++) {
					if (hasInnerPerimeter.get(j) > z) {
						data4Imputation = j;
						break;					
					}
				}
									
				if (data4Imputation > 0){
					
					//find available data above and below
					int jabove = hasInnerPerimeter.get(data4Imputation - 1);
					int jbelow = hasInnerPerimeter.get(data4Imputation);
					
					//calculate weights according to distance to available finds
					double totalDistance = jbelow - jabove;
					double wAbove = (double)(jbelow - z) / totalDistance;
					double wBelow = (double)(z - jabove) / totalDistance;
					
					colCon.mendedCC.ixmid[z] = wAbove * colCon.preciseCC.ixmid[jabove] + wBelow * colCon.preciseCC.ixmid[jbelow];
					colCon.mendedCC.iymid[z] = wAbove * colCon.preciseCC.iymid[jabove] + wBelow * colCon.preciseCC.iymid[jbelow];
					colCon.mendedCC.innerMinorRadius[z] = wAbove * colCon.preciseCC.innerMinorRadius[jabove] + wBelow * colCon.preciseCC.innerMinorRadius[jbelow];
					colCon.mendedCC.innerMajorRadius[z] = wAbove * colCon.preciseCC.innerMajorRadius[jabove] + wBelow * colCon.preciseCC.innerMajorRadius[jbelow];
					colCon.mendedCC.itheta[z] = wAbove * colCon.preciseCC.itheta[jabove] + wBelow * colCon.preciseCC.itheta[jbelow];
					
				}
				if (data4Imputation == 0) {
					
					colCon.mendedCC.ixmid[z] = colCon.preciseCC.ixmid[data4Imputation];
					colCon.mendedCC.iymid[z] = colCon.preciseCC.iymid[data4Imputation];						
					colCon.mendedCC.innerMinorRadius[z] = colCon.preciseCC.innerMinorRadius[data4Imputation];
					colCon.mendedCC.innerMajorRadius[z] = colCon.preciseCC.innerMajorRadius[data4Imputation];
					colCon.mendedCC.itheta[z] = colCon.preciseCC.itheta[data4Imputation];	
					
					
				}
				if (data4Imputation < 0) {
					
					colCon.mendedCC.ixmid[z] = colCon.preciseCC.ixmid[hasInnerPerimeter.size() - 1];
					colCon.mendedCC.iymid[z] = colCon.preciseCC.iymid[hasInnerPerimeter.size() - 1];
					colCon.mendedCC.innerMinorRadius[z] = colCon.preciseCC.innerMinorRadius[hasInnerPerimeter.size() - 1];
					colCon.mendedCC.innerMajorRadius[z] = colCon.preciseCC.innerMajorRadius[hasInnerPerimeter.size() - 1];
					colCon.mendedCC.itheta[z] = colCon.preciseCC.itheta[hasInnerPerimeter.size() - 1];
					
				}
			}	
		}
		
		//impute the missing values for inner perimeter...
		for (int z = colCon.preciseCC.topOfColumn ; z < colCon.preciseCC.bottomOfColumn ; z++) {
			
			if (!hasOuterPerimeter.contains(z)) {			
			
				//find existing entries above and below the missing value
				int data4Imputation = -1;
				for (int j = 0 ; j < hasOuterPerimeter.size() ; j++) {
					if (hasOuterPerimeter.get(j) > z) {
						data4Imputation = j;
						break;					
					}
				}
				
				//find available data above	
				if (data4Imputation > 0){
					
					//find available data above and below
					int jabove = hasOuterPerimeter.get(data4Imputation - 1);
					int jbelow = hasOuterPerimeter.get(data4Imputation);
					
					//calculate weights according to distance to available finds
					double totalDistance = jbelow - jabove;
					double wAbove = (double)(jbelow - z) / totalDistance;
					double wBelow = (double)(z - jabove) / totalDistance;
					
					colCon.mendedCC.xmid[z] = wAbove * colCon.preciseCC.xmid[jabove] + wBelow * colCon.preciseCC.xmid[jbelow];
					colCon.mendedCC.ymid[z] = wAbove * colCon.preciseCC.ymid[jabove] + wBelow * colCon.preciseCC.ymid[jbelow];
					
					colCon.mendedCC.theta[z] = wAbove * colCon.preciseCC.theta[jabove] + wBelow * colCon.preciseCC.theta[jbelow];
					
					//take care of bevel
					if (z > colCon.preciseCC.bevelStartsAt) {
						colCon.mendedCC.outerMinorRadius[z] = colCon.mendedCC.innerMinorRadius[z] + fit.evalLinearFunction(colCon.preciseCC.bevelFunction, z);
						colCon.mendedCC.outerMajorRadius[z] = colCon.mendedCC.innerMajorRadius[z] + fit.evalLinearFunction(colCon.preciseCC.bevelFunction, z);
					}
					else {
						colCon.mendedCC.outerMinorRadius[z] = wAbove * colCon.preciseCC.outerMinorRadius[jabove] + wBelow * colCon.preciseCC.outerMinorRadius[jbelow];
						colCon.mendedCC.outerMajorRadius[z] = wAbove * colCon.preciseCC.outerMajorRadius[jabove] + wBelow * colCon.preciseCC.outerMajorRadius[jbelow];
					}
					
				}
				if (data4Imputation == 0) {
					
					colCon.mendedCC.xmid[z] = colCon.preciseCC.xmid[data4Imputation];
					colCon.mendedCC.ymid[z] = colCon.preciseCC.ymid[data4Imputation];
					colCon.mendedCC.theta[z] = colCon.preciseCC.theta[data4Imputation];
					
					//take care of bevel
					if (z > colCon.preciseCC.bevelStartsAt) {
						colCon.mendedCC.outerMinorRadius[z] = colCon.mendedCC.innerMinorRadius[z] + fit.evalLinearFunction(colCon.preciseCC.bevelFunction, z);
						colCon.mendedCC.outerMajorRadius[z] = colCon.mendedCC.innerMajorRadius[z] + fit.evalLinearFunction(colCon.preciseCC.bevelFunction, z);
					}
					else {
						colCon.mendedCC.outerMinorRadius[z] = colCon.preciseCC.outerMinorRadius[data4Imputation];
						colCon.mendedCC.outerMajorRadius[z] = colCon.preciseCC.outerMajorRadius[data4Imputation];						
					}
					
					
				}
				if (data4Imputation < 0) {
					
					colCon.mendedCC.xmid[z] = colCon.preciseCC.xmid[hasOuterPerimeter.size() - 1];
					colCon.mendedCC.ymid[z] = colCon.preciseCC.ymid[hasOuterPerimeter.size() - 1];					
					colCon.mendedCC.theta[z] = colCon.preciseCC.theta[hasOuterPerimeter.size() - 1];
					
					//take care of bevel
					if (z > colCon.preciseCC.bevelStartsAt) {
						colCon.mendedCC.outerMinorRadius[z] = colCon.mendedCC.innerMinorRadius[z] + fit.evalLinearFunction(colCon.preciseCC.bevelFunction, z);
						colCon.mendedCC.outerMajorRadius[z] = colCon.mendedCC.innerMajorRadius[z] + fit.evalLinearFunction(colCon.preciseCC.bevelFunction, z);
					}
					else{
						colCon.mendedCC.outerMinorRadius[z] = colCon.preciseCC.outerMinorRadius[hasOuterPerimeter.size() - 1];
						colCon.mendedCC.outerMajorRadius[z] = colCon.preciseCC.outerMajorRadius[hasOuterPerimeter.size() - 1];
					}
				
				}
			}	
		}
		
		//re-calculate wall thickness
		for (int z = colCon.preciseCC.topOfColumn ; z < colCon.preciseCC.bottomOfColumn ; z++) {
			colCon.mendedCC.wallThickness[z] = 0.5 * (colCon.mendedCC.outerMajorRadius[z] - colCon.mendedCC.innerMajorRadius[z] + colCon.mendedCC.outerMinorRadius[z] - colCon.mendedCC.innerMinorRadius[z]);
		}
		
		return colCon;
	}
	
	public ColumnContainer rubberBandDetection(ColumnContainer colCon) {
		
		//check if there are layers in between where the wall was not properly detected (because of dirt or rubber band)
		double[] checker = new double[colCon.preciseCC.zmid.length];
		int diffCheckWindowSize = 1;
			
		double[] xf2 = makeChecker(colCon.preciseCC.xmid, diffCheckWindowSize);
		double[] yf2 = makeChecker(colCon.preciseCC.ymid, diffCheckWindowSize);
		double[] mif2 = makeChecker(colCon.preciseCC.outerMinorRadius, diffCheckWindowSize);
		double[] maf2 = makeChecker(colCon.preciseCC.outerMajorRadius, diffCheckWindowSize);		
		for (int i = 0 ; i < checker.length ; i++) checker[i] = Math.sqrt(xf2[i] + yf2[i] + mif2[i] + maf2[i]) / 4;				
		
		//find high and low probability outlines
		int filterSize = 16;
		int[] goodones = findSmallCheckers(checker, filterSize);
		
		//kick out values above and below column
		//for (int i = 0 ; i < goodones.length ; i++) if (i < colCon.preciseCC.topOfColumn | i > colCon.preciseCC.bottomOfColumn) goodones[i] = 0;
		
		//apply correction		
	    int MaximumGoodValueInVicinityNumber = 5;
	    colCon.preciseCC = correct4RubberBandArtifacts(colCon.preciseCC, goodones, MaximumGoodValueInVicinityNumber);		
	    
	    return colCon;
	}
	
	public double[] makeChecker(double[] values, int diffCheckWindowSize) {
		
		double[] diff = new double[values.length];
		diff[values.length - 1] = 0;
		
		for (int i = 0 ; i < values.length - 1; i++) {
			
			ArrayList<Double> diffcheck = new ArrayList<Double>();			
			
			//downstream
			for (int j = i - 1 ; j > i - 2 ; j--) if (j > 0) diffcheck.add((values[j] - values[i]) * (values[j] - values[i]));
			for (int j = i + 1 ; j < i + 2 ; j++) if (j < values.length) diffcheck.add((values[j] - values[i]) * (values[j] - values[i]));
			
			double[] diffCheckVector = new double[diffcheck.size()];
			for (int j = 0 ; j < diffcheck.size() ; j++) diffCheckVector[j] = diffcheck.get(j);
			
			diff[i] = StatUtils.max(diffCheckVector);
			
		}
		
		return diff;
		
	}
	
	public ColumnCoordinates correct4RubberBandArtifacts(ColumnCoordinates preciseCC, int[] goodone, int MaximumGoodValueInVicinityNumber) {
		
		int j,k,cc;
		ColumnCoordinates superPreciseCC = new ColumnCoordinates();
		
		double[] xmid = preciseCC.xmid;
		double[] ymid = preciseCC.ymid;				
		double[] zmid = preciseCC.zmid;
		double[] ixmid = preciseCC.ixmid; 
		double[] iymid = preciseCC.iymid;	
		double[] outerMinorRadius = preciseCC.outerMinorRadius;
		double[] outerMajorRadius = preciseCC.outerMajorRadius;
		double[] innerMinorRadius = preciseCC.innerMinorRadius; 
		double[] innerMajorRadius = preciseCC.innerMajorRadius;
		double[] wallThickness = preciseCC.wallThickness;
		double[] theta = preciseCC.theta;
		double[] itheta = preciseCC.itheta;
		double[] outerR2 = preciseCC.outerR2;
		double[] innerR2 = preciseCC.innerR2;
		boolean[] columnIsInThisDepth = preciseCC.columnIsInThisDepth;
		int topOfColumn = preciseCC.topOfColumn;
		int heightOfColumn = preciseCC.heightOfColumn;
		int numberOfImputedLayers = preciseCC.numberOfImputedLayers;
		
		for (j = 0 ; j < goodone.length ; j++) {
			if (goodone[j] == 0 & j >= topOfColumn & j <= heightOfColumn + topOfColumn) {
				
				numberOfImputedLayers++;
				
				//search downstream
				k = j - 1; cc = 1;				
				ArrayList<Integer> downStream = new ArrayList<Integer>();
				while (k >= topOfColumn & cc <= MaximumGoodValueInVicinityNumber) {
					if (goodone[k] == 1) {
						downStream.add(k);
						cc++;
					}
					k--;
				}
				
				//search upstream
				k = j + 1; cc = 1;				
				ArrayList<Integer> upStream = new ArrayList<Integer>();
				while (k <= heightOfColumn + topOfColumn & cc <= MaximumGoodValueInVicinityNumber & k < goodone.length) {
					if (goodone[k] == 1) {
						upStream.add(k);
						cc++;
					}
					k++;
				}
				
				//calculate corrected value
				xmid[j] = plugItTogether(downStream, upStream, preciseCC.xmid);
				ymid[j] = plugItTogether(downStream, upStream, preciseCC.ymid);				
				ixmid[j] = plugItTogether(downStream, upStream, preciseCC.ixmid);
				iymid[j] = plugItTogether(downStream, upStream, preciseCC.iymid);
				outerMinorRadius[j] = plugItTogether(downStream, upStream, preciseCC.outerMinorRadius);
				outerMajorRadius[j] = plugItTogether(downStream, upStream, preciseCC.outerMajorRadius);
				innerMinorRadius[j] = plugItTogether(downStream, upStream, preciseCC.innerMinorRadius); 
				innerMajorRadius[j] = plugItTogether(downStream, upStream, preciseCC.innerMajorRadius);
				wallThickness[j] = plugItTogether(downStream, upStream, preciseCC.wallThickness);
				theta[j] = plugItTogether(downStream, upStream, preciseCC.theta);
				itheta[j] = plugItTogether(downStream, upStream, preciseCC.theta);
				outerR2[j] = plugItTogether(downStream, upStream, preciseCC.outerR2);
				innerR2[j] = plugItTogether(downStream, upStream, preciseCC.outerR2);
				
			}
		}
		
		//assign corrected values to output struct
		superPreciseCC.xmid = xmid; 
		superPreciseCC.ymid = ymid;				
		superPreciseCC.zmid = zmid;
		superPreciseCC.ixmid = ixmid; 
		superPreciseCC.iymid = iymid;	
		superPreciseCC.outerMinorRadius = outerMinorRadius;
		superPreciseCC.outerMajorRadius = outerMajorRadius;
		superPreciseCC.innerMinorRadius = innerMinorRadius; 
		superPreciseCC.innerMajorRadius = innerMajorRadius;
		superPreciseCC.wallThickness = wallThickness;
		superPreciseCC.theta = theta;
		superPreciseCC.itheta = theta;
		superPreciseCC.outerR2 = outerR2;
		superPreciseCC.innerR2 = outerR2;
		superPreciseCC.columnIsInThisDepth = columnIsInThisDepth;
		superPreciseCC.topOfColumn = topOfColumn;
		superPreciseCC.heightOfColumn = heightOfColumn;
		superPreciseCC.numberOfImputedLayers = numberOfImputedLayers;
		
		return superPreciseCC;
		
	}
           
    public double plugItTogether(ArrayList<Integer> downStream, ArrayList<Integer> upStream, double[] ccFeature) {
    	
    	Median jMed = new Median();
    	double outValue;
    	double[] imputors =  new double[downStream.size() + upStream.size()];
    	
    	for (int k = 0 ; k < downStream.size() ; k++) imputors[k] = ccFeature[downStream.get(k)];
    	for (int k = downStream.size() ; k < imputors.length; k++) imputors[k] = ccFeature[upStream.get(k - downStream.size())];
    	
    	outValue = jMed.evaluate(imputors);
    	
    	return outValue;
    }

	public int[] findSmallCheckers(double[] checker, int filterSize) {
		
	      int[] goodone = new int[checker.length];
	   
	      for (int j = filterSize ; j < checker.length - filterSize ; j++) {
	    	  
	    	  double[] consider = new double[2 * filterSize + 1];
	    	  for (int i = 0 ; i < consider.length ; i++) consider[i] = checker[j - filterSize + i];
	    	  if (StatUtils.max(consider) > 0.1) goodone[j] = 0;
	    	  else goodone[j] = 1;
	      }

	      //also assign values to top and bottom bits of column
	      for (int j = 0 ; j < filterSize ; j++) goodone[j] = 0;
	      for (int j = checker.length - filterSize ; j < checker.length ; j++) goodone[j] = 0;
	      
	      return goodone;
	}
	
	
	public double[] findColumnTilt(double[] xCenter, double[] yCenter, double[] zCenter) {
		
		double[] xyzAngle = new double[3];
		double[] xzLine = new double[3];
		double[] yzLine = new double[3];
		double[] initialGuess = {2, 555};
			
		// first fit in the xz-plane
		CurveFitter xz = new CurveFitter(zCenter, xCenter);
		xz.setInitialParameters(initialGuess);
		xz.doFit(0);
		xzLine = xz.getParams();
				
		// first fit in the yz-plane				
		CurveFitter yz = new CurveFitter(zCenter, yCenter);
		yz.setInitialParameters(initialGuess);
		yz.doFit(0);				
		yzLine = yz.getParams();
			
		//calculate angle of tilt in xz-plane and yz-plane
		xyzAngle[0] = Math.atan(xzLine[1]);
		xyzAngle[1] = Math.atan(yzLine[1]);
		xyzAngle[2] = Math.atan(Math.sqrt(Math.tan(xyzAngle[0]) * Math.tan(xyzAngle[0]) + Math.tan(xyzAngle[1]) * Math.tan(xyzAngle[1])));
		
		return xyzAngle;
	}
	
	public ImageProcessor findSurfaceSubroutine(ImagePlus nowTiff0, PolygonRoi[] pRoi0, int minimumObjectThickness, String topOrBottom) {
		
		JParticleCounter dPC = new JParticleCounter();
		ImageManipulator jIM = new ImageManipulator();
		AndAllTheRest aa = new AndAllTheRest();
		DisplayThings disp = new DisplayThings();
		
		Dilate_ mD = new Dilate_();
		Erode_ mE = new Erode_();
		
		int i;
		int startSlice = 0;
		int stopSlice = 300;
		int slicesPerChunk = 2;		
		double volumeOfChunks2Exclude = 100000;
		PolygonRoi[] pRoi = new PolygonRoi[pRoi0.length];
		
		ImageStack inverseStack = new ImageStack(nowTiff0.getWidth(), nowTiff0.getHeight());
		ImagePlus inverseTop = new ImagePlus();
		ImagePlus nowTiff = new ImagePlus();
	
		if (topOrBottom.contentEquals("bottom")) {
			nowTiff = jIM.flipTiff(nowTiff0);
			for (i = pRoi0.length - 1; i >= 0 ; i--) pRoi[pRoi0.length - 1 - i] = pRoi0[i];			
		} else {
			pRoi = pRoi0;
			nowTiff = nowTiff0;
		}
		
		//find approximate location of upper soil surface
		int approximateTopSurface = 9999;		
		for (i = 0 ; i < nowTiff.getNSlices() ; i++) {
			
			IJ.showStatus("Finding approximate location of soil top surface #" + (i + 1) + "/" + nowTiff.getNSlices());
			
			nowTiff.setPosition(i + 1);			
			ImageProcessor myIP = nowTiff.getProcessor();		
			
			myIP.setRoi(pRoi[i]);
			int[] myHist = myIP.getHistogram();
			float solids = myHist[0];
			float all = myHist[0] + myHist[255];
			float myFractionOfSoil = solids / all;
			if (approximateTopSurface == 9999 & myFractionOfSoil > 0.5) {
				approximateTopSurface = i;
				break;
			}
						
		}
		if (approximateTopSurface > 400) startSlice = approximateTopSurface - 400;
		stopSlice = approximateTopSurface + 300;
		
		//create inverse of binary 
		for (i = startSlice ; i < stopSlice ; i++) {
			
			IJ.showStatus("Inverting slice #" + (i + 1 - startSlice) + "/" + (stopSlice - startSlice));
			
			nowTiff.setPosition(i + 1);			
			ImageProcessor myIP = nowTiff.getProcessor();		
			
			ImageProcessor invIP = jIM.invertSoilBinary(myIP, pRoi[i]);	
			inverseStack.addSlice(invIP);			
		}
		inverseTop.setStack(inverseStack);
		
		IJ.freeMemory();IJ.freeMemory();
		
		//dilate 3-D
		IJ.showStatus("Eroding ...");
		ImagePlus erodedTop = mE.erode(inverseTop, 255, false);
		
		//find the real soil matrix..  (get rid of smaller aggregates on top);
		IJ.showStatus("Removing loose soil aggregates lying on top of the surface using the ParticleAnalyzer  ...");
		Object[] result = dPC.getParticles(erodedTop, slicesPerChunk, volumeOfChunks2Exclude, Double.POSITIVE_INFINITY,-1, false);
		int[][] particleLabels = (int[][]) result[1];
		ImagePlus topMatrix = disp.jDisplayParticleLabels(particleLabels, erodedTop);
				
		//binarize the image again		
		ImageStack purifiedSurfStack = new ImageStack(topMatrix.getWidth(), topMatrix.getHeight());
		for (i = 0 ; i < topMatrix.getNSlices() ; i++) {
			
			topMatrix.setPosition(i + 1);			
			
			ImageProcessor myIP = topMatrix.getProcessor();
			ImageProcessor hisIP = myIP.convertToByte(false);
			hisIP.max(1);
			hisIP.multiply(255);
			purifiedSurfStack.addSlice(hisIP);
		}
		
		ImagePlus purifiedTop = new ImagePlus();
		purifiedTop.setStack(purifiedSurfStack);
	
		//dilate again		
		IJ.showStatus("Dilating ...");
		ImagePlus reDilated = mD.dilate(purifiedTop, 255, false);
			
		IJ.freeMemory();IJ.freeMemory();
		
		//create xy grid
		ImageProcessor outIP = new ShortProcessor(nowTiff.getWidth(), nowTiff.getHeight());
		int[][] xy = new int[outIP.getWidth()][outIP.getHeight()];
		for (int x = 0 ; x < outIP.getWidth() ; x++) {
			for (int y = 0 ; y < outIP.getHeight() ; y++) {
				xy[x][y] = 0;
			}
		}
	
		//find the soil surface
		for (i = 0 ; i < reDilated.getNSlices() ; i++) {
			
			IJ.showStatus("Searching for soil surface in slice #" + (i + 1) + "/" + (stopSlice - startSlice));
			
			reDilated.setPosition(i + 1);			
			
			ImageProcessor myIP = reDilated.getProcessor();
					
			for (int x = 0 ; x < myIP.getWidth() ; x++) {
				for (int y = 0 ; y < myIP.getHeight() ; y++) {					
					if (xy[x][y] == 0) {					//if there was not yet another particle above at this position 
						int thisImgPixel = myIP.getPixel(x, y);
				
						if (thisImgPixel > 0) {
							int cc = 0;
							int theValue2Put = i + startSlice;
							if (i==0) theValue2Put = 9999; //mark the pixels located on the PVC column as zero..
							for (int j = i + 1 ; j < reDilated.getNSlices() ; j++) {
								reDilated.setPosition(j + 1);
								thisImgPixel = myIP.getPixel(x, y);
								if (thisImgPixel > 0) cc++;
								else break;
								if (cc > minimumObjectThickness) {
									xy[x][y] = theValue2Put;									
									break;
								}
							}
						}					
					}				
				}
			}
		}

		//transfer results to image		
		for (int x = 0 ; x < outIP.getWidth() ; x++) {
			for (int y = 0 ; y < outIP.getHeight() ; y++) {
				if (xy[x][y] == 9999) xy[x][y] = 0;
				outIP.putPixelValue(x, y, xy[x][y]);
			}			
		}
		
		return outIP;				
	}
	
	public ImagePlus findSoilSurface(ImagePlus nowTiff, String nowGaugePath, MenuWaiter.SurfaceFinderReturn mSFR) {
				
		//init objects
		ImageManipulator jIM = new ImageManipulator();
		InputOutput jIO = new InputOutput();
		MenuWaiter menu = new MenuWaiter();
		MenuWaiter.ThresholderMenuReturn mTMR = menu.new ThresholderMenuReturn();
		RoiHandler roi = new RoiHandler();
				
		//init varis
		ImagePlus mySurface = new ImagePlus();
		ImagePlus[] rightTiff = {null, null, null};
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		ImageProcessor topIP = new ShortProcessor(nowTiff.getWidth(), nowTiff.getHeight());
		ImageProcessor bottomIP = new ShortProcessor(nowTiff.getWidth(), nowTiff.getHeight());
	
		int minimumObjectThickness = mSFR.neglectCrumbsOfLessThan;
		
		//read gauge file
		ColumnCoordinates jCO = jIO.readGaugeFile(nowGaugePath);	
		PolygonRoi[] pRoi = roi.makeMeAPolygonRoiStack("inner", "wide", jCO, 0);
		
		if (mSFR.greyscale) { 

			//set up thresholding 
			mTMR.airThreshold = mSFR.threshold;		//ok, this naming is misleading.. apologies..
			mTMR.save4Evaluation = false;
			mTMR.useConstantThreshold = true;
			
			rightTiff = jIM.applyAChosenThresholdingMethod(nowTiff, nowGaugePath, mTMR, "");
		
		}
		else rightTiff[0] = nowTiff;
		
		//find top and bottom surface considering the air-phase / non-air-phase boundary
		ImageProcessor prelimTopIP = findSurfaceSubroutine(rightTiff[0], pRoi, minimumObjectThickness, "top");  
		ImageProcessor prelimBottomIP = findSurfaceSubroutine(rightTiff[0], pRoi, minimumObjectThickness, "bottom");
		
		if (mSFR.croppedSoil) {
			
			topIP = jIM.fillHolesWithNeighboringValues(prelimTopIP, pRoi, "top");
			bottomIP = jIM.fillHolesWithNeighboringValues(prelimBottomIP, pRoi, "bottom");	
			
		}
		else {
			topIP = prelimTopIP;
			bottomIP = prelimBottomIP;
		}
			
		//smooth pores that open to the soil surface.. 
		ImageProcessor filledTop = smoothAwaySmallPores(topIP, mSFR.radius4PoreFilter, pRoi[0]);
		ImageProcessor filledBottom = smoothAwaySmallPores(bottomIP, mSFR.radius4PoreFilter, pRoi[0]);
		
		//transfer surface IPs to imagePlus
		outStack.addSlice(filledTop);
		outStack.addSlice(filledBottom);		
		mySurface.setStack(outStack);
		
		//mySurface.draw();mySurface.show();
		
		return mySurface;		
	}
	
	private ImageProcessor smoothAwaySmallPores(ImageProcessor origIP, double smootherRadius, PolygonRoi pRoi) {

		JParticleCounter jPC = new JParticleCounter();
		DisplayThings disp = new DisplayThings();
		AndAllTheRest aa = new AndAllTheRest();
		
		if (smootherRadius > 0) {
		
			for (int i = 0 ; i < (int)Math.round(origIP.getMax()) ; i++) {
				
				IJ.showStatus("Filling pores opening to the surface in layer " + (i + 1));
				
				ImageProcessor smoothIP = origIP.duplicate();	
				smoothIP.min(i);
				smoothIP.max(i + 1);
				smoothIP.subtract(i);
				smoothIP.multiply(255);
				
				ImageProcessor binIP = smoothIP.convertToByte(false);	
				ImagePlus binTiff = new ImagePlus("", binIP);
				
				//binTiff.draw();binTiff.show();
	
				//apply the particle counter to measure the different pore sizes;
				Object[] result = jPC.getParticles(binTiff, 1, 0d, Double.POSITIVE_INFINITY,-1, false);			
				int[][] particleLabels = (int[][]) result[1];
				long[] particleSizes = jPC.getParticleSizes(particleLabels);
				final int nParticles = particleSizes.length;
				int[][] limits = jPC.getParticleLimits(binTiff, particleLabels, nParticles);  //xmin, xmax, ymin, ymax, zmin, zmax
				double[] volumes = jPC.getVolumes(binTiff, particleSizes);
				ImagePlus myHoles = disp.jDisplayParticleLabels(particleLabels, binTiff);	
				
				//spot the holes that should be filled
				ArrayList<Integer> nowOut = new ArrayList<Integer>();
				for (int j = 0 ; j < nParticles ; j++) if (volumes[j] < smootherRadius * smootherRadius * Math.PI) nowOut.add(j);
				
				//fill holes that were spotted
				ImageProcessor holeIP = myHoles.getProcessor();				
				for (int j = 0 ; j < nowOut.size() ; j++) {
					
					int xmin = limits[j][0] - 1;
					int xmax = limits[j][1] + 1;
					int ymin = limits[j][2] - 1;
					int ymax = limits[j][3] + 1;
					
					for (int x = xmin ; x < xmax + 1 ; x++) {
						for (int y = ymin ; y < ymax + 1 ; y++) {
							int holeCheck = (int)Math.round(holeIP.getPixelValue(x, y));
							if (holeCheck == nowOut.get(j)) {
								origIP.putPixel(x, y, i);	//if is hole than fill it..
							}
						}
					}
				}
				
				IJ.freeMemory();IJ.freeMemory();				
			}
		
		}
		
		return origIP;
		
	}
	
	public double[] calcHeightAndDiameterOfColumn(ImagePlus soilSurface, int numberOfSlices) {
		
		double[] hod = new double[2];
		int iW = soilSurface.getWidth();
		int iH = soilSurface.getHeight();
		
		//calculate top of soil
		soilSurface.setPosition(1);
		ImageProcessor surIP = soilSurface.getProcessor();
		double depthSum = 0;
		double area = 0;
		for (int x = 0 ; x < iW ; x++) {
			for (int y = 0 ; y < iH ; y++) {
				int nowPixel = (int)surIP.getPixelValue(x, y);
				if (nowPixel > 0) {
					area++;
					depthSum += nowPixel;
				}
			}
		}
		double topMinus = depthSum / area;
		hod[1] = 2 * Math.sqrt(area / Math.PI);
		
		//calculate top of soil
		soilSurface.setPosition(2);
		surIP = soilSurface.getProcessor();
		depthSum = 0;
		area = 0;
		for (int x = 0 ; x < iW ; x++) {
			for (int y = 0 ; y < iH ; y++) {
				int nowPixel = (int)surIP.getPixelValue(x, y);
				if (nowPixel > 0) {
					area++;
					depthSum += nowPixel;
				}
			}
		}
		
		double botMinus = depthSum / area;
		
		hod[0] = numberOfSlices - topMinus - botMinus;		
		
		return hod;
	}
	
	public double[] findMedianWallGreyValues(ImagePlus nowTiff, ObjectDetector.ColumnCoordinates preciseCC) {

		RoiHandler roi = new RoiHandler();
		
		int i;		
		
		ImagePlus origTiff = nowTiff.duplicate();
		
		PolygonRoi[] iRoi = roi.makeMeAPolygonRoiStack4RawImage("inner", "wide", preciseCC, 0);
		PolygonRoi[] oRoi = roi.makeMeAPolygonRoiStack4RawImage("outer", "tight", preciseCC, 0);
		
		double[] wallMedian = new double[preciseCC.heightOfColumn];
		
		//sample column wall grey
		for (i = preciseCC.topOfColumn ; i < preciseCC.topOfColumn + preciseCC.heightOfColumn ; i++) { 
			
			IJ.showStatus("Sampling illumination of slice " + (i - preciseCC.topOfColumn + 1) + "/" + (preciseCC.heightOfColumn));	
			
			origTiff.setPosition(i+1);
			ImageProcessor myIP = origTiff.getProcessor();		
			
			//origTiff.draw();
			//origTiff.show();
					
			myIP.setRoi(oRoi[i - preciseCC.topOfColumn]);
			myIP.setValue(0);
			//myIP.setColor(Color.white);
			//myIP.draw(oRoi[i]);
			myIP.fillOutside(oRoi[i - preciseCC.topOfColumn]);		
			
			myIP.setRoi(iRoi[i - preciseCC.topOfColumn]);
			myIP.setValue(0);
			//myIP.setColor(Color.white);
			//myIP.draw(iRoi[i]);
			myIP.fill(iRoi[i - preciseCC.topOfColumn]);	
			
			myIP.resetRoi();
			
			//origTiff.draw();
			//origTiff.show();
				
			int[] myHist = myIP.getHistogram();
			float[] cHist = new float[myHist.length];	
			cHist[0] = 0;
			for (int j = 1 ; j < myHist.length ; j++) cHist[j] = cHist[j-1] + (float)myHist[j];
			for (int j = 1 ; j < myHist.length ; j++) cHist[j] = cHist[j] / cHist[cHist.length - 1];			
			for (int j = 1 ; j < myHist.length ; j++) if (cHist[j] > 0.5) {
				wallMedian[i - preciseCC.topOfColumn] = j;
				break;				
			}
	
		}
		
/*		double[] xax =  new double[preciseCC.heightOfColumn];
		for (i = 0 ; i < preciseCC.heightOfColumn ; i++) xax[i] = i; 		
		Plot nP = new Plot("","","",xax, wallMedian);
		nP.setLimits(0, xax.length, 6500, 7500);
		nP.draw();
		nP.show();*/
		
		return wallMedian;
		
	}	
	
	public int[] postProcessSteelGreyValues(double[] medianOfWallGreyValue) {
		
		int[] topBotHeight = new int[3];
		Median jMed = new Median();
		
		int steelMedian = (int)Math.round(jMed.evaluate(medianOfWallGreyValue));	
		
		int firstGood = -1;
		int lastGood = medianOfWallGreyValue.length - 1;
		for (int i = 0 ; i < medianOfWallGreyValue.length ; i++) {
			if (firstGood == -1) if (medianOfWallGreyValue[i] > 0.9 * steelMedian) firstGood = i;
			if (firstGood > 0 & lastGood == medianOfWallGreyValue.length - 1) if (medianOfWallGreyValue[i] < 0.9 * steelMedian) lastGood = i - 1;
		}
				
		topBotHeight[0] = firstGood;
		topBotHeight[1] = lastGood; 
		topBotHeight[2] = lastGood - firstGood + 1;
		
		return topBotHeight;
		
	}
	
	/*public double[] findMedianSteelGreyValues(ImagePlus nowTiff, ObjectDetector.ColumnCoordinates preciseCC) {
		
		int i;		

		ImagePlus origTiff = nowTiff.duplicate();
	
		double[] wallQ = new double[nowTiff.getNSlices()];
		double[] q1 = new double[nowTiff.getNSlices()];		
		
		//sample column wall grey
		for (i = 0 ; i < nowTiff.getNSlices() ; i++) { 
			
			IJ.showStatus("Sampling illumination of slice " + (i + 1) + "/" + nowTiff.getNSlices());	
			
			origTiff.setPosition(i+1);
			ImageProcessor myIP = origTiff.getProcessor();		
			
			float[] myX = new float[preciseCC.xID[i].length];for (int j = 0 ; j < myX.length ; j++) myX[j] = preciseCC.xID[i][j];
			float[] myY = new float[preciseCC.yID[i].length];for (int j = 0 ; j < myY.length ; j++) myY[j] = preciseCC.yID[i][j];
			PolygonRoi iRoi = new PolygonRoi(myX,myY,Roi.POLYLINE);
			float[] moX = new float[preciseCC.xOD[i].length];for (int j = 0 ; j < moX.length ; j++) moX[j] = preciseCC.xOD[i][j];
			float[] moY = new float[preciseCC.yOD[i].length];for (int j = 0 ; j < moY.length ; j++) moY[j] = preciseCC.yOD[i][j];			
			PolygonRoi oRoi = new PolygonRoi(moX,moY,Roi.POLYLINE);
					
			//fit a spline to it..
			//iRoi.fitSpline();
			//oRoi.fitSpline();
			
			//trick imagej to fill the outside of oRoi
			ImageProcessor copyIP = myIP.duplicate();
			copyIP.setValue(0);
			copyIP.fill();
			
			copyIP.setRoi(oRoi);
			copyIP.setValue(1);
			copyIP.fill(oRoi);
			
			copyIP.setRoi(iRoi);
			copyIP.setValue(0);
			copyIP.fill(iRoi);
			
			//multiply both image processors
			for (int xx = 0 ; xx < myIP.getWidth() ; xx++) {
				for (int yy = 0 ; yy < myIP.getHeight() ; yy++) {
					if (copyIP.getPixel(xx, yy) == 0) myIP.putPixel(xx, yy, 0);
				}
			}
			
			//origTiff.updateAndDraw();			
			//origTiff.show();
				
			int[] myHist = myIP.getHistogram();
			float[] cHist = new float[myHist.length];	
			cHist[0] = 0;
			for (int j = 1 ; j < myHist.length ; j++) cHist[j] = cHist[j-1] + (float)myHist[j];
			for (int j = 1 ; j < myHist.length ; j++) cHist[j] = cHist[j] / cHist[cHist.length - 1];			
			
			//evaluate thickness of wall
			//double wallThickness = Math.sqrt((preciseCC.xOD[i][0] - preciseCC.xID[i][0]) * (preciseCC.xOD[i][0] -preciseCC.xID[i][0]) + (preciseCC.yOD[i][0] -preciseCC.yID[i][0]) * (preciseCC.yOD[i][0] -preciseCC.yID[i][0])); 			 
			
			for (int j = 1 ; j < myHist.length ; j++) if (cHist[j] > 0.5) {
				q1[i] = j;
				break;				
			}
			
			for (int j = 1 ; j < myHist.length ; j++) if (cHist[j] > 0.9) {
				wallQ[i] = j;
				break;			
			}
			
		}
		
		//filter samples with dark q1
		//double mq1 = jMed.evaluate(q1);
		//for (int j = 0 ; j < q1.length ; j++) if (q1[j] < 0.5 * mq1) wallQ[j] = 0;
		
		//origTiff.updateAndDraw();			
		//origTiff.show();
		
		return wallQ;
		
	}	*/

	public MorphologyAnalyzer.SurfaceStatistics extractSurfaceStatistics(String myOutPath, String myTiff, ImagePlus nowTiff, String nowGaugePath, int numberOfSlices) {

		MorphologyAnalyzer morph = new MorphologyAnalyzer();
		MorphologyAnalyzer.SurfaceStatistics mSS = morph.new SurfaceStatistics();
		RoiHandler jRH = new RoiHandler();	
		InputOutput jIO = new InputOutput();
		HistogramStuff hist = new HistogramStuff();
		
		//read gauge file
		ColumnCoordinates jCO = jIO.readGaugeFile(nowGaugePath);	
		PolygonRoi[] pRoi = jRH.makeMeAPolygonRoiStack("inner", "wide", jCO, 0);
		
		//get top stats
		nowTiff.setPosition(1);
		ImageProcessor topIP = nowTiff.getProcessor();
		topIP.setRoi(pRoi[0]);		
		int[] topHist = topIP.getHistogram();
		topHist[0] = 0;
		
		mSS.highestElevation = hist.findMinFromHistogram(topHist);
		mSS.medianElevation = hist.findMedianFromHistogram(topHist);
		mSS.meanElevation = (int)Math.round(hist.findMeanFromHistogram(topHist));
		mSS.lowestElevation = hist.findMaxFromHistogram(topHist);
		
		//get bottom stats
		nowTiff.setPosition(2);
		ImageProcessor bottomIP = nowTiff.getProcessor();
		bottomIP.setRoi(pRoi[pRoi.length - 1]);
	
		int[] bottomHist = bottomIP.getHistogram();
		mSS.highestIntrusion = numberOfSlices - hist.findMaxFromHistogram(bottomHist);
		mSS.medianIntrusion = numberOfSlices - hist.findMedianFromHistogram(bottomHist);
		mSS.meanIntrusion = numberOfSlices - (int)Math.round(hist.findMeanFromHistogram(bottomHist));
		mSS.lowestIntrusion = numberOfSlices - hist.findMinFromHistogram(bottomHist);
		
		return mSS;
				
	}
	
	public double calculateCVOfWallThickness(FitStuff.FittedEllipse fI, FitStuff.FittedEllipse fO) {
		
		RoiHandler roi = new RoiHandler();
		
		PolygonRoi iRoi = roi.makeRoiFromFittedEllipse(fI);
		PolygonRoi oRoi = roi.makeRoiFromFittedEllipse(fO);
		
		int[] iX = iRoi.getXCoordinates();
		int[] iY = iRoi.getYCoordinates();
		int[] oX = oRoi.getXCoordinates();
		int[] oY = oRoi.getYCoordinates();
		
		double ixb = iRoi.getXBase();
		double iyb = iRoi.getYBase();
		double oxb = oRoi.getXBase();
		double oyb = oRoi.getYBase();
		
		double[] wallThickness = new double[iX.length];
		
		for (int i = 0 ; i < iX.length ; i++) {
			
			double dX = ((double)oX[i] + oxb) - ((double)iX[i] + ixb);
			double dY = ((double)oY[i] + oyb) - ((double)iY[i] + iyb);
			
			wallThickness[i] = Math.sqrt(dX * dX + dY * dY);
			
		}
		
		double stdev = Math.sqrt(StatUtils.variance(wallThickness));
		double avg = StatUtils.mean(wallThickness);
		
		double CV = stdev / avg;
		
		return CV;
		
	}

}



