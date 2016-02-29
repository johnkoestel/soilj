package soilj.tools;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.Filters3D;
import ij.plugin.GaussianBlur3D;
import ij.plugin.PlugIn;
import ij.plugin.Slicer;
import ij.plugin.filter.GaussianBlur;
import ij.process.AutoThresholder;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import soilj.tools.ObjectDetector;
import soilj.tools.ObjectDetector.ColumnContainer;

import java.awt.Polygon;
import java.util.ArrayList;
import java.util.Arrays;

import java.util.Random;

import org.apache.commons.math3.analysis.function.Logistic;
import org.apache.commons.math3.analysis.function.Logistic.Parametric;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import soilj.copiedTools.JParticleCounter;

public class ImageManipulator implements PlugIn {
	
	// carries out all sorts of image manipulations (illumination correction, changes tilt of object, applying filters..)
	
	public void run(String arg) {
		//nothing will be run here..		
	}
	
	///////////////////////////////////////////
	// classes
	///////////////////////////////////////////
	
	public class StackCalculator {
		
		public ImagePlus add(ImagePlus a, ImagePlus b) {
			
			ImageStack outStack = new ImageStack(a.getWidth(), a.getHeight()); 
			ImagePlus outTiff = new ImagePlus();
			
			for (int z = 0 ; z < a.getNSlices() ; z++) {
				
				a.setPosition(z + 1);
				b.setPosition(z + 1);
				
				ImageProcessor aIP = a.getProcessor();
				ImageProcessor bIP = b.getProcessor();
				
				ImageProcessor outIP = aIP.duplicate();
				
				for (int x = 0 ; x < a.getWidth() ; x++) {
					for (int y = 0 ; y < a.getHeight() ; y++) {
						int avox = aIP.getPixel(x, y);
						int bvox = bIP.getPixel(x, y);
						outIP.putPixel(x, y, avox + bvox);
					}
				}
				
				outStack.addSlice(outIP);
			}
			
			outTiff.setStack(outStack);
			
			return outTiff;
		}
		
		public ImagePlus subtract(ImagePlus a, ImagePlus b) {
			
			ImageStack outStack = new ImageStack(a.getWidth(), a.getHeight()); 
			ImagePlus outTiff = new ImagePlus();
			
			for (int z = 0 ; z < a.getNSlices() ; z++) {
				
				a.setPosition(z + 1);
				b.setPosition(z + 1);
				
				ImageProcessor aIP = a.getProcessor();
				ImageProcessor bIP = b.getProcessor();
				
				ImageProcessor outIP = aIP.duplicate();
				
				for (int x = 0 ; x < a.getWidth() ; x++) {
					for (int y = 0 ; y < a.getHeight() ; y++) {
						int avox = aIP.getPixel(x, y);
						int bvox = bIP.getPixel(x, y);
						outIP.putPixel(x, y, avox - bvox);
					}
				}
				
				outStack.addSlice(outIP);
			}
			
			outTiff.setStack(outStack);
			
			return outTiff;
		}
		
	}
		
	///////////////////////////////////////////////
	// Melitta
	//////////////////////////////////////////////
	
	public ImagePlus applyMedianFilterAndUnsharpMask(ImagePlus nowTiff, MenuWaiter.MedianFilterAndUnsharpMaskReturn mMUS) {				
		StackCalculator mSC = new StackCalculator();
		
		//construct some objects
		ImageStack zStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImageStack readyStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());		
		
		//define filter codes			
		int median3DFilter = 11; // 11 is the code for the 3D median filter..	
	
		//apply 3-D median filter		
		readyStack = nowTiff.getStack();
		if (mMUS.medianFilterSizeXDir > 1 & mMUS.medianFilterSizeYDir > 1 & mMUS.medianFilterSizeZDir > 1 ) {
			zStack = Filters3D.filter(readyStack, median3DFilter, mMUS.medianFilterSizeXDir, mMUS.medianFilterSizeYDir, mMUS.medianFilterSizeZDir);
		}
		else {
			zStack = readyStack;
		}
		ImagePlus zTiff = new ImagePlus();
		zTiff.setStack(zStack);		
	
		//apply 3-D unsharp mask
		ImagePlus blurTiff = zTiff.duplicate();
		GaussianBlur3D.blur(blurTiff, mMUS.uMaskStandardDeviationXDir, mMUS.uMaskStandardDeviationYDir, mMUS.uMaskStandardDeviationZDir);
		
		//calculate difference between blurTiff and original		
		ImagePlus diffTiff = mSC.subtract(zTiff,blurTiff);
		
		//and weighted mask
		ImageStack weightedStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		for (int i = 0 ; i < diffTiff.getNSlices() ; i++) {
			diffTiff.setPosition(i+1);
			ImageProcessor nowIP = diffTiff.getProcessor();
			nowIP.multiply(mMUS.uMaskSharpeningWeight);
			weightedStack.addSlice(nowIP);
		}
		ImagePlus unsharpMask = new ImagePlus();
		unsharpMask.setStack(weightedStack);
		
		//sharpened Image
		ImagePlus filtTiff = mSC.add(zTiff, unsharpMask);
		
		//return filtered 3D image
		return filtTiff;
	}
		
	///////////////////////////////////////////////
	// geo and illu correct
	//////////////////////////////////////////////

	public ColumnContainer innerCircleAndTrimAluOrPVCColumns(ImagePlus nowTiff, MenuWaiter.ColumnFinderMenuReturn jCFS) {
		
		//construct some objects
		ObjectDetector jOD = new ObjectDetector();
		ObjectDetector.ColumnCoordinates prelimCC;		
		
		//re-determine the column outlines
		prelimCC = jOD.findOrientationOfPVCOrAluColumn(nowTiff, jCFS);
		
		IJ.freeMemory();IJ.freeMemory();
		
		//find precise column Outlines				
		ObjectDetector.ColumnContainer colCon = jOD.findPVCOrAluWalls(nowTiff, prelimCC, jCFS);
		
		IJ.freeMemory();IJ.freeMemory();
				
		//cut column according to found top and bottom
		colCon.cutTiff = cutColumnsEndsOff(colCon);

		colCon.nowTiff.flush();
		IJ.freeMemory();IJ.freeMemory();
		
		return colCon;
		
	}
	
	/*public ColumnContainer carveOutSteelColumn(ImagePlus nowTiff, double wallThickness) {
	
		//construct some objects
		Median jMed = new Median();
		
		ObjectDetector jOD = new ObjectDetector();
		ObjectDetector.ColumnCoordinates prelimCC = null;
		ObjectDetector.ColumnCoordinates preciseCC = null;
		ObjectDetector.ColumnContainer colCon = jOD.new ColumnContainer();
		
		ImagePlus upTiff = new ImagePlus();
	
		//find the inner edge of the PVC column
		prelimCC = jOD.findOrientationSteelColumn(nowTiff);
		
		//init colCon
		colCon.nowTiff = nowTiff;
		colCon.prelimCC = prelimCC;
		
		//check if column is already upright					
		upTiff = putSteelColumnUprightInCenter(nowTiff, prelimCC);
		
		//upTiff.draw();
		//upTiff.show();
		
		//re-determine the column outlines
		colCon.preciseCC = jOD.findAllSteelWalls(upTiff, prelimCC, wallThickness);
	
		//set top of the soil column to 0
		preciseCC.topOfColumn = 0;
	
		//sample wall illumination and save it	
		double[] medianOfWallGreyValue0 = jOD.findMedianSteelGreyValues(upTiff, preciseCC);
		
		//evaluate grey values..
		int[] topBotHeight = jOD.postProcessSteelGreyValues(medianOfWallGreyValue0);
		
		//shift top and bottom accordingly..
		preciseCC.topOfColumn = topBotHeight[0];
		preciseCC.bottomOfColumn = topBotHeight[1];
		preciseCC.heightOfColumn = topBotHeight[2];
		double[] medianOfWallGreyValue = new double[preciseCC.heightOfColumn];
		for (int i = preciseCC.topOfColumn ; i < preciseCC.bottomOfColumn ; i++) medianOfWallGreyValue[i - preciseCC.topOfColumn] = medianOfWallGreyValue0[i]; 		
		preciseCC.steelGreyValue = jMed.evaluate(medianOfWallGreyValue);
		preciseCC.wallGreyValues = medianOfWallGreyValue0;
		
		//do illumination correction and cut away top and bottom
		colCon.cutTiff = doSteelIlluCorrAndCutOffEnds(upTiff, preciseCC);
			
		return colCon;
	
	}*/
	
	public ImagePlus beamDeHardenThis(ImagePlus nowTiff, String nowGaugePath, MenuWaiter.BeamDeHardeningReturn mBDH) {
		
		//init objects
		InputOutput jIO = new InputOutput();
		ObjectDetector jOD = new ObjectDetector(); 
		FitStuff fittings = new FitStuff();
		ImageManipulator jIM = new ImageManipulator();
		
		//init variables
		ImagePlus outTiff = new ImagePlus();	
		ObjectDetector.ColumnCoordinates jPCO = null;		

		//read steelGauge file
		/*if (mBDH.isSteelColumn) {
			jSCO = jIO.readSteelGaugeFile(nowGaugePath);
			int anglesChecked = jSCO.anglesChecked;
			double[][][] radialGreyValues = jOD.getRadialSteelGreyValues(nowTiff, jSCO, radialMappingFactor);			
		}*/
		//else {
			
			//load gauge file
			jPCO = jIO.readGaugeFile(nowGaugePath);		
			
			//get 1D sample's brightness 
			//ObjectDetector.SampleBrightness myBri = jOD.getSampleBrightness(nowTiff, jPCO);
			
			//get radial illumination
			//radialGreyValues = jOD.getRadialPVCGreyValues(nowTiff, jPCO, radialMappingFactor, mBDH);
			
			//get modes of radial histograms
			ObjectDetector.RadialModes myModes = jOD.getRadialPVCIllumination(nowTiff, jPCO, mBDH);
			
			//fit a curve to the matrix brightness data..
			FitStuff.FittingResults myLogisticFit = fittings.fitGeneralizedLogisticFunction(myModes, mBDH.maxEval);
						
			//fill gaps and bad fits
			double[][] filteredGLFParams = jIM.fillGapsInFittedBeamHardeningCorrectionMap(myLogisticFit, myModes, mBDH);
			
			//assemble beam hardening correction map
			ImageProcessor blurIP = assembleCorrectionMapForBeamHardening(myModes, filteredGLFParams);
			
			//apply brightness correction to collected minima..
			myModes = applyBeamHardeningCorrection2AirPhaseGammaData(myModes, blurIP);
			
			//fit a curve to the air-phase brightness data..
			FitStuff.FittingResults myHyperbolicFit = fittings.fitHyperbolicFunction2Params(myModes, mBDH.maxEval);
			
			//fill gaps and bad fits
			double[][] filteredHFParams = jIM.fillGapsInFittedAirPhaseGammaCorrectionMap(myLogisticFit, myHyperbolicFit, myModes, mBDH);
			
			//assemble air phase gamma correction map
			ImageProcessor gammaIP = assembleCorrectionMapForAirPhaseGamma(myModes, filteredHFParams);

			//correct image for beam hardening and air-phase gamma
			outTiff = correct4PVCBeamHardening(nowTiff, blurIP, gammaIP, jPCO, mBDH);
			
		//}
		
		//apply correction
		/*if (mBDH.isSteelColumn) {
			double[][][] bhc = jIO.readSteelBeamHardeningCorrectionParameters(nowGaugePath, nowTiff.getNSlices());
			outTiff = correct4SteelBeamHardening(nowTiff, jSCO, bhc, standardRadius);
		}
		else {			
			//
		} */
		
		//return corrected tiffs..
		return outTiff;
		
	}
	
	
	public ObjectDetector.RadialModes applyBeamHardeningCorrection2AirPhaseGammaData(ObjectDetector.RadialModes myModes, ImageProcessor blurIP) { 
		
		for (int i = 0 ; i < myModes.maskingThreshold.length ; i++) {
			for (int j = 0 ; j < myModes.radius.length - 1 ; j++) {
				double nowCorrFac = blurIP.getPixelValue((int)Math.floor(myModes.radius[j]) + 1, i);
				double nowRef = blurIP.getPixelValue(blurIP.getWidth() - 1, i);
				myModes.maskedRadialMinima[i][j] = myModes.maskedRadialMinima[i][j] * nowRef / nowCorrFac;
			}
		}
		
		return myModes;
	}	
	
	public ImageProcessor assembleCorrectionMapForBeamHardening(ObjectDetector.RadialModes myModes, double[][] filteredGLFParams) {
		
		TailoredMaths maths = new TailoredMaths();
		
		//create a correction function map
		double[] radius = new double[(int)maths.max(myModes.radius)];
		for (int i = 0 ; i < radius.length ; i++) radius[i] = i;
		Parametric myLogL = new Logistic.Parametric();
		FloatProcessor corrIP = new FloatProcessor(radius.length, myModes.maskingThreshold.length);
		for (int y = 0 ; y < myModes.maskingThreshold.length ; y++) {
			double[] nowGLFP = new double[6];
			for (int i = 0 ; i < 6 ; i++) nowGLFP[i] = filteredGLFParams[y][i];
			for (int x = 0 ; x < radius.length ; x++) {	
				double myValue2Put = myModes.maskedRadialModes[y][myModes.radius.length - 1];
				try {myValue2Put = myLogL.value(x, nowGLFP);} 
				catch (NotStrictlyPositiveException e) {} 
				corrIP.putPixelValue(x, y, myValue2Put);
			}
		}
		
		//ImagePlus letsSee = new ImagePlus("nonFiltered", corrIP);
		//letsSee.updateAndDraw();letsSee.show();
		
		GaussianBlur myGB = new GaussianBlur();
		double accuracy = 0.01d;
		FloatProcessor blurIP = (FloatProcessor)corrIP.duplicate();
		myGB.blurFloat(blurIP, 20, 20, accuracy);			
		//ImagePlus letsSeeM = new ImagePlus("Filtered",blurIP);
		//letsSeeM.updateAndDraw();letsSeeM.show();
		
		return blurIP;
	}
	
	public ImageProcessor assembleCorrectionMapForAirPhaseGamma(ObjectDetector.RadialModes myModes, double[][] myGapFilledHFParameters) {
		
		FitStuff fittings = new FitStuff();
		TailoredMaths maths = new TailoredMaths();
		
		//create a correction function map
		double[] radius = new double[(int)maths.max(myModes.radius)];
		for (int i = 0 ; i < radius.length ; i++) radius[i] = i;
		FitStuff.HyperbolicFunction2Params myHF = fittings.new HyperbolicFunction2Params();
		FloatProcessor corrIP = new FloatProcessor(radius.length, myModes.maskingThreshold.length);
		for (int y = 0 ; y < myModes.maskingThreshold.length ; y++) {
			double[] nowHFparams = new double[2];
			myHF.setup(StatUtils.max(myModes.radius), myModes.maskedRadialMinima[y][myModes.radius.length - 1]);
			for (int i = 0 ; i < 2 ; i++) nowHFparams[i] = myGapFilledHFParameters[y][i];
			for (int x = 0 ; x < radius.length ; x++) {	
				double myValue2Put = myModes.maskedRadialModes[y][myModes.radius.length - 1];
				try {myValue2Put = myHF.value(x, nowHFparams);} 
				catch (NotStrictlyPositiveException e) {} 
				corrIP.putPixelValue(x, y, myValue2Put);
			}
		}
		
		//ImagePlus letsSee = new ImagePlus("nonFiltered", corrIP);
		//letsSee.updateAndDraw();letsSee.show();
		
		GaussianBlur myGB = new GaussianBlur();
		double accuracy = 0.01d;
		FloatProcessor blurIP = (FloatProcessor)corrIP.duplicate();
		myGB.blurFloat(blurIP, 20, 20, accuracy);			
		//ImagePlus letsSeeM = new ImagePlus("Filtered",blurIP);
		//letsSeeM.updateAndDraw();letsSeeM.show();
		
		return blurIP;
	}
	
	public double[][] fillGapsInFittedBeamHardeningCorrectionMap(FitStuff.FittingResults myFit, ObjectDetector.RadialModes myModes, MenuWaiter.BeamDeHardeningReturn mBDH) {
		
		RollerCaster cast = new RollerCaster();
		
		double[] R2 = myFit.R2;
		int standardRadius = (int)Math.round(StatUtils.max(myModes.radius));
		double[] realRadius = new double[standardRadius];
		ArrayList<Integer> testOnes = new ArrayList<Integer>();
		ArrayList<Integer> goodOnes = new ArrayList<Integer>();		
		
		double[][] gapFilledParameters = new double[R2.length][myFit.numberOfParams];
		
		//sample the maximal beamhardening artifact in the soil
		for (int i = R2.length/3 ; i < 2 * R2.length / 3 ; i++) if (R2[i] > mBDH.severeGoodnessCriterion) testOnes.add(i);		
		double[] A = new double[testOnes.size()];			//minimum of radial brightness profile
		double[] K = new double[testOnes.size()];			//maximum of radial brightness profile
		for (int i = 0 ; i < testOnes.size() ; i++) {
			A[i] = myFit.params[testOnes.get(i)][4];
			K[i] = myFit.params[testOnes.get(i)][0];
		}
		double expectedBrightnessDifference = StatUtils.percentile(K, 50) - StatUtils.percentile(A, 50);
		
		//set all R2 of fits with more than the 2 expected brightness difference to 0;
		for (int i = 0 ; i < R2.length ; i++) if (myFit.params[i][0] - myFit.params[i][4] > expectedBrightnessDifference | myFit.params[i][0] - myFit.params[i][4] < 0) {
			R2[i] = 0;
		}
		
		//find first and last good fit
		int firstGood = 0;
		int lastGood = 0;
		for (int i = 0 ; i < R2.length ; i++) {
			if (firstGood == 0 & R2[i] > mBDH.severeGoodnessCriterion) firstGood = i;
			if (R2[i] > mBDH.severeGoodnessCriterion) lastGood = i;
		}
		
		//plug in close-to-wall values before firstGood and after lastGood
		for (int i = 0 ; i < firstGood ; i++) {
			myFit.params[i][0] = myModes.maskedRadialModes[i][myModes.radius.length - 1];
			myFit.params[i][1] = myFit.params[firstGood][1];
			myFit.params[i][2] = myFit.params[firstGood][2];
			myFit.params[i][3] = myFit.params[firstGood][3];
			myFit.params[i][4] = myModes.maskedRadialModes[i][myModes.radius.length - 1] - 1;
			myFit.params[i][5] = myFit.params[firstGood][5];
			R2[i] = mBDH.goodnessCriterion + 0.01;
		}
		for (int i = lastGood + 1 ; i < R2.length ; i++) {
			myFit.params[i][0] = myModes.maskedRadialModes[i][myModes.radius.length - 1];
			myFit.params[i][1] = myFit.params[lastGood][1];
			myFit.params[i][2] = myFit.params[lastGood][2];
			myFit.params[i][3] = myFit.params[lastGood][3];
			myFit.params[i][4] = myModes.maskedRadialModes[i][myModes.radius.length - 1] - 1;
			myFit.params[i][5] = myFit.params[lastGood][5];
			R2[i] = mBDH.goodnessCriterion + 0.01;
		}
		
		//pick out the good ones and imputed ones
		for (int i = 0 ; i < R2.length ; i++) if (R2[i] > mBDH.goodnessCriterion) goodOnes.add(i);
		
		//re-mould them in vectors and matrices
		int[] depth = new int[goodOnes.size()];
		double[][] data = new double[goodOnes.size()][myFit.numberOfParams];		
		for (int i = 0 ; i < depth.length ; i++) {
			depth[i] = goodOnes.get(i);
			for (int j = 0 ; j < myFit.numberOfParams ; j++) data[i][j] = myFit.params[depth[i]][j];
		}
		
		//init realRadius
		for (int i = 0 ; i < realRadius.length ; i++) realRadius[i] = i;
		
		//interpolate the gaps for all parameters
		for (int i = 0 ; i < myFit.numberOfParams ; i++) {			
			double[] nowData = new double[depth.length];
			for (int j = 0 ; j < nowData.length ; j++) nowData[j] = data[j][i];
			LinearInterpolator myLI = new LinearInterpolator();		
			PolynomialSplineFunction mySF = myLI.interpolate(cast.castInt2Double(depth), nowData);			 
			for (int j = 0 ; j < R2.length ; j++) gapFilledParameters[j][i] = mySF.value(j);						
		}
		
		return gapFilledParameters;
		
	}
	
	public double[][] fillGapsInFittedAirPhaseGammaCorrectionMap(FitStuff.FittingResults myMatrixBrightnessFit, FitStuff.FittingResults myHFFit, ObjectDetector.RadialModes myModes, MenuWaiter.BeamDeHardeningReturn mBDH) {
		
		FitStuff fittings = new FitStuff();
		RollerCaster cast = new RollerCaster();
		
		
		double[] R2 = myHFFit.R2;
		double[] matrixR2 = myMatrixBrightnessFit.R2;
		int standardRadius = (int)Math.round(StatUtils.max(myModes.radius));
		double[] realRadius = new double[standardRadius];
		ArrayList<Integer> testOnes = new ArrayList<Integer>();
		ArrayList<Integer> goodOnes = new ArrayList<Integer>();		
		
		double[][] gapFilledParameters = new double[R2.length][myHFFit.numberOfParams];
		
		//sample the maximal beam-hardening artifact in the soil ... using the matrix fits here..
		for (int i = matrixR2.length/3 ; i < 2 * matrixR2.length / 3 ; i++) if (matrixR2[i] > mBDH.severeGoodnessCriterion) testOnes.add(i);		
		double[] A = new double[testOnes.size()];			//minimum of radial brightness profile
		double[] K = new double[testOnes.size()];			//maximum of radial brightness profile
		for (int i = 0 ; i < testOnes.size() ; i++) {
			A[i] = myMatrixBrightnessFit.params[testOnes.get(i)][4];
			K[i] = myMatrixBrightnessFit.params[testOnes.get(i)][0];
		}
		double expectedBrightnessDifference = StatUtils.percentile(K, 50) - StatUtils.percentile(A, 50);
		
		//set all R2 of fits with more than the 2 expected brightness difference to 0;
		for (int i = 0 ; i < R2.length ; i++) {
			if (myMatrixBrightnessFit.params[i][0] - myMatrixBrightnessFit.params[i][4] > expectedBrightnessDifference 
					| myMatrixBrightnessFit.params[i][0] - myMatrixBrightnessFit.params[i][4] < 0) {		
				R2[i] = 0;
			}
		}
		
		//perform a similar check for the air-phase gamma correction function
		FitStuff.HyperbolicFunction2Params myHF = fittings.new HyperbolicFunction2Params();		
		for (int i = 0 ; i < R2.length ; i++) {
			myHF.setup(StatUtils.max(myModes.radius), myModes.maskedRadialMinima[i][myModes.radius.length - 1]);
			double refMatrixBrightness = myModes.maskedRadialModes[i][myModes.radius.length - 1]; 
			double refAirPhaseBrightness = myModes.maskedRadialMinima[i][myModes.radius.length - 1];
			double[] nowParams = new double[2];			
			for (int j = 0 ; j < 2 ; j++) nowParams[j] = myHFFit.params[i][j];
			double fittedAirPhaseBrightness = myHF.value(0, nowParams);
			if ( (refMatrixBrightness - refAirPhaseBrightness) / (refMatrixBrightness - fittedAirPhaseBrightness) > 3) {		
				R2[i] = 0;
			}
		}		
		
		//find first and last good fit .. using the matrix fits here..
		int firstGood = 0;
		int lastGood = 0;
		for (int i = 0 ; i < matrixR2.length ; i++) {
			if (firstGood == 0 & matrixR2[i] > mBDH.severeGoodnessCriterion) firstGood = i;
			if (matrixR2[i] > mBDH.severeGoodnessCriterion) lastGood = i;
		}
		
		//plug in close-to-wall values before firstGood and after lastGood
		for (int i = 0 ; i < firstGood ; i++) {			
			myHFFit.params[i][1] = 0;
			R2[i] = mBDH.gammaGoodnessCriterion + 0.01;
		}
		for (int i = lastGood + 1 ; i < matrixR2.length ; i++) {			
			myHFFit.params[i][1] = 0;
			R2[i] = mBDH.gammaGoodnessCriterion + 0.01;
		}
		
		//pick out the good ones and imputed ones
		for (int i = 0 ; i < R2.length ; i++) if (R2[i] > mBDH.gammaGoodnessCriterion) goodOnes.add(i);
		
		//re-mould them in vectors and matrices
		int[] depth = new int[goodOnes.size()];
		double[][] data = new double[goodOnes.size()][myHFFit.numberOfParams];		
		for (int i = 0 ; i < depth.length ; i++) {
			depth[i] = goodOnes.get(i);
			for (int j = 0 ; j < myHFFit.numberOfParams ; j++) data[i][j] = myHFFit.params[depth[i]][j];
		}
		
		//init realRadius
		for (int i = 0 ; i < realRadius.length ; i++) realRadius[i] = i;
		
		//interpolate the gaps for all parameters
		for (int i = 0 ; i < myHFFit.numberOfParams ; i++) {			
			double[] nowData = new double[depth.length];
			for (int j = 0 ; j < nowData.length ; j++) nowData[j] = data[j][i];
			LinearInterpolator myLI = new LinearInterpolator();		
			PolynomialSplineFunction mySF = myLI.interpolate(cast.castInt2Double(depth), nowData);			 
			for (int j = 0 ; j < R2.length ; j++) gapFilledParameters[j][i] = mySF.value(j);						
		}
		
		return gapFilledParameters;
		
	}
	
	public ImagePlus clipImage(int i, ImagePlus nowTiff, InputOutput.MyFolderCollection mFC, MenuWaiter.ClipperMenuReturn mSCM) {
		
		InputOutput jIO = new InputOutput(); 
		ObjectDetector jOD = new ObjectDetector(); 
		RoiHandler roi = new RoiHandler();
		
		ImagePlus outTiff = new ImagePlus();		
		
		String nowGaugePath = null;
		String nowSurfPath = null;
		int[] myGandS = new int[2];
		ObjectDetector.ColumnCoordinates jCO = jOD.new ColumnCoordinates();		
		PolygonRoi[] pRoi = new PolygonRoi[nowTiff.getNSlices()];		
		int xmin = nowTiff.getWidth();
		int xmax = 0;
		int ymin = nowTiff.getHeight();
		int ymax = 0;
		
		int referenceSlice = 1;		
		
		//load gauge and surface files
		if (mSCM.isSoilColumn == true) {			
			if (mSCM.referenceIsTopMostSlice) myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC.myTiffs[i], mFC.myGaugePaths, null);			
			if (mSCM.referenceIsSoilSurface) {
				myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC.myTiffs[i], mFC.myGaugePaths, mFC.mySurfaceFileNames);
				nowSurfPath = mFC.mySurfaceFolder + "//" + mFC.mySurfaceFileNames[myGandS[1]];
			}
			nowGaugePath = mFC.myGaugePaths[myGandS[0]];
			
			jCO = jIO.readGaugeFile(nowGaugePath);
			pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, mSCM.clipFromInnerPerimeter);
					
			if (mSCM.referenceIsSoilSurface == true) {
				ImagePlus mySurfs = jIO.openTiff3D(nowSurfPath);	
				referenceSlice = jOD.findMedianSoilSurfacePosition(mySurfs);	
			}
			
			//find min and max values for clip
			for (int j = referenceSlice + mSCM.startAtSlice ; j < referenceSlice + mSCM.stopAtSlice ; j++) {
				
				Polygon p = pRoi[j-2].getNonSplineCoordinates();
				
				// for x			
				int[] nowX = new int[pRoi[j].getNCoordinates()]; 				
				nowX = p.xpoints;
				Arrays.sort(nowX, 0, nowX.length); 
				int xmi = nowX[0] + (int)Math.round(pRoi[j].getXBase());
				int xma = nowX[nowX.length - 1] + (int)Math.round(pRoi[j].getXBase());
				
				// and for y
				int[] nowY = new int[pRoi[j].getNCoordinates()]; 				
				nowY = p.ypoints;
				Arrays.sort(nowY, 0, nowY.length); 
				int ymi = nowY[0] + (int)Math.round(pRoi[j].getYBase());
				int yma = nowY[nowY.length - 1] + (int)Math.round(pRoi[j].getYBase());
				
				// remember the mins and maxes
				if (xmi < xmin) xmin = xmi - mSCM.canvasExceedsBy;
				if (xma > xmax) xmax = xma + mSCM.canvasExceedsBy;
				if (ymi < ymin) ymin = ymi - mSCM.canvasExceedsBy;
				if (yma > ymax) ymax = yma + mSCM.canvasExceedsBy;				
			}			
		} 
		else {			// if the sample is not a soil column..			
			xmin = 0 + mSCM.clipFromCanvasEdge - mSCM.canvasExceedsBy;
			xmax = nowTiff.getWidth() - mSCM.clipFromCanvasEdge + mSCM.canvasExceedsBy;
			ymin = 0 + mSCM.clipFromCanvasEdge - mSCM.canvasExceedsBy;
			ymax = nowTiff.getHeight() - mSCM.clipFromCanvasEdge + mSCM.canvasExceedsBy;			
		}			
		
		//create clip-out-roi
		float[] xpf = {xmin, xmax, xmax, xmin, xmin};
		float[] ypf = {ymin, ymin, ymax, ymax, ymin};	
		PolygonRoi coRoi = new PolygonRoi(xpf, ypf, Roi.POLYGON);
		
		//create simple Roi in case that the image does not contain a soil column
		float[] xps = {xmin + mSCM.canvasExceedsBy, xmax - mSCM.canvasExceedsBy, xmax - mSCM.canvasExceedsBy, xmin + mSCM.canvasExceedsBy, xmin + mSCM.canvasExceedsBy};
		float[] yps = {ymin + mSCM.canvasExceedsBy, ymin + mSCM.canvasExceedsBy, ymax - mSCM.canvasExceedsBy, ymax - mSCM.canvasExceedsBy, ymin + mSCM.canvasExceedsBy};
		PolygonRoi simpleRoi = new PolygonRoi(xps, yps, Roi.POLYGON);
		
		//clip and cut
		ImageStack outStack = null;
		if (!mSCM.preserveOriginialCanvasSize) {
		
			outStack = new ImageStack(xmax - xmin, ymax - ymin);
		
			for (int j = referenceSlice + mSCM.startAtSlice ; j < referenceSlice + mSCM.stopAtSlice ; j++) {
			
				nowTiff.setPosition(j);
				ImageProcessor nowIP = nowTiff.getProcessor();
			
				if (mSCM.isSoilColumn == true) simpleRoi = pRoi[j];
			
				nowIP.setRoi(simpleRoi);
				nowIP.setBackgroundValue(0);
				nowIP.fillOutside(simpleRoi);
			
				nowIP.setRoi(coRoi);
				ImageProcessor cropIP = nowIP.crop();
			
				outStack.addSlice(cropIP);
			}			
		} 
		else {
			
			outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
			
			for (int j = referenceSlice + mSCM.startAtSlice ; j < referenceSlice + mSCM.stopAtSlice ; j++) {
			
				nowTiff.setPosition(j);
				ImageProcessor nowIP = nowTiff.getProcessor();
			
				if (mSCM.isSoilColumn == true) simpleRoi = pRoi[j];
			
				nowIP.setRoi(simpleRoi);
				nowIP.setBackgroundValue(0);
				nowIP.fillOutside(simpleRoi);
			
				outStack.addSlice(nowIP);
			}
		}			 
		
		//add blank canvas to bottom if it is wished..
		if (mSCM.addCanvasExccedance2Bottom == true) {			
			ImageProcessor blankIP = null;
			if (nowTiff.getBitDepth() == 8) blankIP = new ByteProcessor(outStack.getWidth(), outStack.getHeight());
			if (nowTiff.getBitDepth() == 16) blankIP = new ShortProcessor(outStack.getWidth(), outStack.getHeight());
			for (int j = 0 ; j < mSCM.canvasExceedsBy ; j++) {
				outStack.addSlice(blankIP);	
			}
		}
		
		outTiff.setStack(outStack);
		
		//outTiff.updateAndDraw();
		//outTiff.show();
		
		return outTiff;
		
	}
	
	/*public ImagePlus correct4SteelBeamHardening(ImagePlus nowTiff, ObjectDetector.ColumnCoordinates jCO, double[][][] bhc, int standardRadius) {
		
		//init output variables
		ImagePlus outTiff = new ImagePlus();
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		//calculate number of angles needed to paint the whole canvas..
		double spacing = 2 * Math.PI / jCO.anglesChecked;		
		float[] sR = new float[standardRadius]; 
		float[] sA = new float[jCO.anglesChecked + 1];
		
		//create array for standardradius
		for (int i = 0 ; i < standardRadius ; i++) sR[i] = i + 1;
		for (int i = 0 ; i < jCO.anglesChecked + 1 ; i++) sA[i] = i;		 

		//sweep over slices
		for (int i = 0 ; i < nowTiff.getNSlices() ; i++) {
			
			IJ.showStatus("Correcting for beam hardening in slice #" + (i + 1) + "/" + nowTiff.getNSlices());
			
			//set image to next slice
			nowTiff.setPosition(i+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			//calculate normalized correction functions
			double r01 = bhc[i][3][0];
			float[] ftot01 = correctionFunctionCalculator(i, bhc, 0, sR);
			
			//double r40 = bhc[i][3][1];
			//float[] ftot40 = correctionFunctionCalculator(i, bhc, 1, sR);
			
			//double r60 = bhc[i][3][2];
			//float[] ftot60 = correctionFunctionCalculator(i, bhc, 2, sR);
			
			double r80 = bhc[i][3][3];
			float[] ftot80 = correctionFunctionCalculator(i, bhc, 3, sR);
			
			//do the correction
			for (int x = 0 ; x < nowIP.getWidth() ; x++) {
				for (int y = 0 ; y < nowIP.getHeight() ; y++) {			
					
					//get distance of pixel to center
					float dx = x - (float)jCO.xmid[i];
					float dy = y - (float)jCO.ymid[i];
					float nowRadius = (float)Math.sqrt((double)(dx*dx) + (double)(dy*dy));
					
					//get angle
					double alpha = Math.atan(dy/dx);
					if (dx < 0 & dy >= 0) alpha = 2 * Math.PI + alpha;
					if (dy < 0) alpha = Math.PI + alpha;
															
					//get radius at this angle
					double dx0 = jCO.xID[i][0] - jCO.xmid[i];
					double dy0 = jCO.yID[i][0] - jCO.ymid[i];
					double radiusAtThisAngle = Math.sqrt((dx0*dx0) + (dy0*dy0));;
					for (double j = 0 ; j < jCO.anglesChecked ; j++) {
						double nowAngle  = 2 * j / jCO.anglesChecked * Math.PI;
						if (nowAngle > alpha) {
							
							double weight = (nowAngle - alpha) / spacing;							
							
							double dx1 = jCO.xID[i][(int)j-1] - jCO.xmid[i];
							double dx2 = jCO.xID[i][(int)j] - jCO.xmid[i];
							double dy1 = jCO.yID[i][(int)j-1] - jCO.ymid[i];
							double dy2 = jCO.yID[i][(int)j] - jCO.ymid[i];
							double r1 = Math.sqrt((dx1*dx1) + (dy1*dy1));
							double r2 = Math.sqrt((dx2*dx2) + (dy2*dy2));
							
							radiusAtThisAngle = weight * r2 + (1 - weight) * r1;
							
							break;
						}
					}
					
					//check if pixel is within column and if yes apply correction
					if (nowRadius < radiusAtThisAngle) {
						
						int renormalizedRadius = (int)Math.floor(nowRadius / radiusAtThisAngle * standardRadius);
						int C01 = (int)Math.round(ftot01[renormalizedRadius]);
						//int C40 = (int)Math.round(ftot40[renormalizedRadius]);
						//int C60 = (int)Math.round(ftot60[renormalizedRadius]);
						int C80 = (int)Math.round(ftot80[renormalizedRadius]);
						
						float myGrey = nowIP.getPixelValue(x, y);
						//int newGrey = (int)Math.round(r60 + myGrey - C60);
						//int newGrey = (int)Math.round(r60 + ((r80-r60)/(C80-C60)) * (myGrey - C60));
						int newGrey = (int)Math.round(r01 + ((r80-r01)/(C80-C01)) * (myGrey - C01));
						
						nowIP.putPixel(x, y, newGrey);
						
					}			
				}	
		
			}
			
			//nowTiff.updateAndDraw();
			//nowTiff.show();

			outStack.addSlice(nowIP);
		}
		
		outTiff.setStack(outStack);
		
		//outTiff.updateAndDraw();
		//outTiff.show();
		
		return outTiff;		
		
	}*/
	
	public ImagePlus correct4PVCBeamHardening(ImagePlus nowTiff, ImageProcessor blurIP, ImageProcessor gammaIP, ObjectDetector.ColumnCoordinates jCO, MenuWaiter.BeamDeHardeningReturn mBDH) {
		
		RoiHandler roi = new RoiHandler();
		
		//init output variables
		ImagePlus outTiff = new ImagePlus();
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		//calculate number of angles needed to paint the whole canvas..
		double spacing = 2 * Math.PI / mBDH.anglesChecked;
		int standardRadius = blurIP.getWidth();
		float[] sR = new float[standardRadius]; 
		float[] sA = new float[mBDH.anglesChecked + 1];
		
		//create array for standardradius
		for (int i = 0 ; i < standardRadius ; i++) sR[i] = i + 1;
		for (int i = 0 ; i < mBDH.anglesChecked + 1 ; i++) sA[i] = i;	
		
		//load polygon roi of inner perimeter
		PolygonRoi[] nowRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 2);

		//sweep over slices
		for (int i = 0 ; i < nowTiff.getNSlices() ; i++) {
			
			IJ.showStatus("Correcting for beam hardening in slice #" + (i + 1) + "/" + nowTiff.getNSlices());
			
			//extract X and Y of pRoi		
			nowRoi[i].fitSpline(mBDH.anglesChecked);
			Polygon myPoly = nowRoi[i].getPolygon();			
			int xID[] = myPoly.xpoints;
			int yID[] = myPoly.ypoints;
			
			//set image to next slice
			nowTiff.setPosition(i+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			//get matrix brightness reference values for this depth			
			double[] corrFunc = new double[standardRadius];
			double nowReference = blurIP.getPixelValue(standardRadius - 1, i);
			for (int r = 0 ; r < standardRadius ; r++) corrFunc[r] = nowReference / blurIP.getPixelValue(r, i);
			
			//get air-phase gamma reference values for this depth			
			double[] gammaFunc = new double[standardRadius];
			double nowRefGamma = gammaIP.getPixelValue(standardRadius - 1, i);
			double gammaDelta = nowReference - nowRefGamma;
			for (int r = 0 ; r < standardRadius ; r++) gammaFunc[r] = gammaDelta / (nowReference - gammaIP.getPixelValue(r, i));
					
			//do the correction
			for (int x = 0 ; x < nowIP.getWidth() ; x++) {
				for (int y = 0 ; y < nowIP.getHeight() ; y++) {			
					
					//get distance of pixel to center
					float dx = x - (float)jCO.xmid[i];
					float dy = y - (float)jCO.ymid[i];
					float nowRadius = (float)Math.sqrt((double)(dx*dx) + (double)(dy*dy));
					
					//get angle
					double alpha = Math.atan(dy/dx);
					if (dx < 0 & dy >= 0) alpha = 2 * Math.PI + alpha;
					if (dy < 0) alpha = Math.PI + alpha;
															
					//get radius at this angle
					double dx0 = xID[0] - jCO.xmid[i];
					double dy0 = yID[0] - jCO.ymid[i];
					double radiusAtThisAngle = Math.sqrt((dx0*dx0) + (dy0*dy0));;
					for (double j = 0 ; j < mBDH.anglesChecked ; j++) {
						double nowAngle  = 2 * j / mBDH.anglesChecked * Math.PI;
						if (nowAngle > alpha) {
							
							double weight = (nowAngle - alpha) / spacing;							
							
							double dx1 = xID[(int)j-1] - jCO.xmid[i];
							double dx2 = xID[(int)j] - jCO.xmid[i];
							double dy1 = yID[(int)j-1] - jCO.ymid[i];
							double dy2 = yID[(int)j] - jCO.ymid[i];
							double r1 = Math.sqrt((dx1*dx1) + (dy1*dy1));
							double r2 = Math.sqrt((dx2*dx2) + (dy2*dy2));
							
							radiusAtThisAngle = weight * r2 + (1 - weight) * r1;
							
							break;
						}
					}
					
					//check if pixel is within column and if yes apply correction
					double corrFactor = 1;
					double gammaFactor = 1;
					if (nowRadius < radiusAtThisAngle) {						
						int renormalizedRadius = (int)Math.floor(nowRadius / radiusAtThisAngle * standardRadius);
						corrFactor = corrFunc[renormalizedRadius];
						gammaFactor = gammaFunc[renormalizedRadius];
					}		
					
					//apply beam hardening correction
					float myGrey = nowIP.getPixelValue(x, y);					
					double newGrey = (int)Math.round(corrFactor * myGrey);
					
					//also apply air-phase gamma correction
					int newestGrey = (int)Math.round(nowReference - (nowReference - newGrey) * gammaFactor); 
					
					nowIP.putPixel(x, y, newestGrey);
				}	
		
			}
	
			outStack.addSlice(nowIP);
		}
		
		outTiff.setStack(outStack);
		
		return outTiff;		
		
	}
	
	public float[] correctionFunctionCalculator(int i, double[][][] bhc, int bhcEntry, float[] sR) {
		
		float[] ftot = new float[sR.length];
		int standardRadius = sR.length;
		
		double a = bhc[i][0][bhcEntry];
		double b = bhc[i][1][bhcEntry];
		double c = bhc[i][2][bhcEntry];
		double r = bhc[i][3][bhcEntry];
		double dy = bhc[i][4][bhcEntry];		
		double f1, f2, f3;
		
		for (int j = 0 ; j < standardRadius ; j++) {
			f1 = dy / Math.exp(a * sR[j] - a) + r;
			f2 = dy / Math.exp(b * sR[j] - b) + r;
			f3 = dy / Math.exp(c * sR[j] - c) + r;
			ftot[standardRadius - j - 1] = (float)((f1 + f2 + f3) / 3);  //switch vector around!!!! important!!
		}
		
		return ftot;
	}
	
	/*public ImagePlus putSteelColumnUprightInCenter(ImagePlus nowTiff, ObjectDetector.ColumnCoordinates prelimCC) {
		
		ObjectDetector jOD = new ObjectDetector();		
		Slicer sly = new Slicer();
		Median jMed = new Median();
		
		ObjectDetector.ColumnCoordinates newCC = jOD.new ColumnCoordinates();
		
		ImagePlus straightTiff = new ImagePlus();
		ImagePlus outTiff = new ImagePlus();
		
		//by-pass the cutting step when re-measuring column wall..
		//outTiff = nowTiff;
		
		if (prelimCC.tiltTotal > 0.01) {   //if tilting angle is too big then put column straight
			
			IJ.showStatus("Putting column straight and moving it to center of canvas, step 1 of 4...");
			
			ImagePlus boTiff0 = new ImagePlus();
			ImageStack boStack0 = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
			//find tilting angle relative to XY-plane
			double alpha = prelimCC.tiltInXZ;
			double beta = prelimCC.tiltInYZ;
			double dz = 1000;
			
			double dx = Math.tan(alpha) * dz;
			double dy = Math.tan(beta) * dz;
			double gamma = Math.atan(dx/dy);
			
			//correct for ambiguity in atan response			
			if (dy < 0) gamma = gamma + Math.PI; 
			
			//rotate so that the tilting is only in y-direction
			for (int i = 1 ; i < nowTiff.getNSlices() + 1 ; i++) {
				
				nowTiff.setPosition(i);
				ImageProcessor nowIP = nowTiff.getProcessor();
				nowIP.setInterpolate(true);
				double gammaInDegree = gamma * 360 / 2 / Math.PI ;
				nowIP.rotate(gammaInDegree + 90);
				
				boStack0.addSlice(nowIP);
			}
			
			boTiff0.setStack(boStack0);			
						
			//reslice
			ImagePlus vertiTiff = sly.reslice(boTiff0);
			
			IJ.showStatus("Putting column straight and moving it to center of canvas, step 2 of 4...");
			
			//put column upright
			ImagePlus boTiff = new ImagePlus();
			ImageStack boStack = new ImageStack(vertiTiff.getWidth(), vertiTiff.getHeight());			
			for (int i = 1 ; i < vertiTiff.getNSlices() + 1 ; i++) {
				
				vertiTiff.setPosition(i);
				ImageProcessor nowIP = vertiTiff.getProcessor();
				nowIP.setInterpolate(true);
				double deltaInDegree = prelimCC.tiltTotal * 360 / 2 / Math.PI ;
				nowIP.rotate(-deltaInDegree);
				
				boStack.addSlice(nowIP);
			}
		
			boTiff.setStack(boStack);	
			
			//reslice back to XY-plane view
			ImagePlus straightTiff0 = sly.reslice(boTiff);
			
			IJ.showStatus("Putting column straight and moving it to center of canvas, step 3 of 4...");
			
			//rotate back to original orientation
			ImageStack straightStack = new ImageStack(straightTiff0.getWidth(), straightTiff0.getHeight());
			for (int i = 1 ; i < straightTiff0.getNSlices() + 1 ; i++) {
				
				straightTiff0.setPosition(i);
				ImageProcessor nowIP = straightTiff0.getProcessor();
				nowIP.setInterpolate(true);
				double gammaInDegree = gamma * 360 / 2 / Math.PI ;
				nowIP.rotate(-gammaInDegree - 90);
				
				straightStack.addSlice(nowIP);
			}
			
			straightTiff.setStack(straightStack);

		}
		else {
			straightTiff = nowTiff;
		}
		
		//straightTiff.show();
		//straightTiff.draw();
		
		//re-find column outlines
		newCC = jOD.findOrientationSteelColumn(straightTiff);
		
		//move column into the center of the canvas and cut out unnecessary parts of the canvas
		IJ.showStatus("Putting column straight and moving it to center of canvas, step 4 of 4...");
		double rim = 25;
		double mRadius = jMed.evaluate(newCC.medianOuterRadius);
		double toBeLeft = mRadius + rim;
		
		//check if cut-out Roi is not larger than the image..
		if (2 * toBeLeft > straightTiff.getWidth()) toBeLeft = straightTiff.getWidth() / 2;
		if (2 * toBeLeft > straightTiff.getHeight()) toBeLeft = straightTiff.getHeight() / 2;
			
		//do the cutting
		double diameter = 2 * toBeLeft;
		ImageStack outStack = new ImageStack((int)Math.round(diameter), (int)Math.round(diameter));
		//Boolean cutSuccessful = true;
		for (int i = 1 ; i < straightTiff.getNSlices() + 1 ; i++) {
			
			IJ.showStatus("Adding straightened slice " + i + "/" + straightTiff.getNSlices());
			
			straightTiff.setPosition(i);
			ImageProcessor nowIP = straightTiff.getProcessor();
			nowIP.setInterpolate(true);
			
			double XC = straightTiff.getWidth() / 2;
			double YC = straightTiff.getHeight() / 2;
			
			double dx = XC - jMed.evaluate(newCC.xmid);
			double dy = YC - jMed.evaluate(newCC.ymid);
			
			nowIP.translate(dx, dy);
			
			Roi cutRoi = new Roi(XC - toBeLeft, YC - toBeLeft, (int)Math.round(diameter), (int)Math.round(diameter));			
			nowIP.setRoi(cutRoi);
			ImageProcessor cutIP = nowIP.crop();
			
			outStack.addSlice(cutIP);
		}
		
		outTiff.setStack(outStack);
		
		
		//outTiff.draw();
		//outTiff.show();		
		
		return outTiff;
		
	}*/
	
	public ImagePlus putColumnUprightInCenter(ImagePlus nowTiff, ObjectDetector.ColumnCoordinates prelimCC, MenuWaiter.ColumnFinderMenuReturn jCFS) {
		
		ObjectDetector jOD = new ObjectDetector();		
		Slicer sly = new Slicer();
		Median jMed = new Median();
		
		ObjectDetector.ColumnCoordinates newCC = jOD.new ColumnCoordinates();
		
		ImagePlus straightTiff = new ImagePlus();
		ImagePlus outTiff = new ImagePlus();
		
		//by-pass the cutting step when re-measuring column wall..
		//outTiff = nowTiff;
		
		if (prelimCC.tiltTotal > 0.01) {   //if tilting angle is too big then put column straight
			
			IJ.showStatus("Putting column straight and moving it to center of canvas, step 1 of 4...");
			
			ImagePlus boTiff0 = new ImagePlus();
			ImageStack boStack0 = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
			//find tilting angle relative to XY-plane
			double alpha = prelimCC.tiltInXZ;
			double beta = prelimCC.tiltInYZ;
			double dz = 1000;
			
			double dx = Math.tan(alpha) * dz;
			double dy = Math.tan(beta) * dz;
			double gamma = Math.atan(dx/dy);
			
			//correct for ambiguity in atan response			
			if (dy < 0) gamma = gamma + Math.PI; 
			
			//rotate so that the tilting is only in y-direction
			for (int i = 1 ; i < nowTiff.getNSlices() + 1 ; i++) {
				
				nowTiff.setPosition(i);
				ImageProcessor nowIP = nowTiff.getProcessor();
				nowIP.setInterpolate(true);
				double gammaInDegree = gamma * 360 / 2 / Math.PI ;
				nowIP.rotate(gammaInDegree + 90);
				
				boStack0.addSlice(nowIP);
			}
			
			boTiff0.setStack(boStack0);			
						
			//reslice
			ImagePlus vertiTiff = sly.reslice(boTiff0);
						
			//free memory ..
			boTiff0.flush();
			IJ.freeMemory();IJ.freeMemory();
			
			IJ.showStatus("Putting column straight and moving it to center of canvas, step 2 of 4...");
			
			//put column upright
			ImagePlus boTiff = new ImagePlus();
			ImageStack boStack = new ImageStack(vertiTiff.getWidth(), vertiTiff.getHeight());			
			for (int i = 1 ; i < vertiTiff.getNSlices() + 1 ; i++) {
				
				vertiTiff.setPosition(i);
				ImageProcessor nowIP = vertiTiff.getProcessor();
				nowIP.setInterpolate(true);
				double deltaInDegree = prelimCC.tiltTotal * 360 / 2 / Math.PI ;
				nowIP.rotate(-deltaInDegree);
				
				boStack.addSlice(nowIP);
			}

			//free memory ..
			vertiTiff.flush();
			IJ.freeMemory();IJ.freeMemory();
			
			//transfer images to boTiff
			boTiff.setStack(boStack);	
			
			//reslice back to XY-plane view
			ImagePlus straightTiff0 = sly.reslice(boTiff);
			
			//free memory from boTiff again..
			boTiff.flush();
			IJ.freeMemory();IJ.freeMemory();
			
			IJ.showStatus("Putting column straight and moving it to center of canvas, step 3 of 4...");
			
			//rotate back to original orientation
			ImageStack straightStack = new ImageStack(straightTiff0.getWidth(), straightTiff0.getHeight());
			for (int i = 1 ; i < straightTiff0.getNSlices() + 1 ; i++) {
				
				straightTiff0.setPosition(i);
				ImageProcessor nowIP = straightTiff0.getProcessor();
				nowIP.setInterpolate(true);
				double gammaInDegree = gamma * 360 / 2 / Math.PI ;
				nowIP.rotate(-gammaInDegree - 90);
				
				straightStack.addSlice(nowIP);
			}
			
			//free memory from straightTiff0 again..
			straightTiff0.flush();
			IJ.freeMemory();IJ.freeMemory();
			
			straightTiff.setStack(straightStack);

		}
		else {
			straightTiff = nowTiff;
		}
				
		//straightTiff.draw();
		//straightTiff.show();
		
		//re-find column outlines
		newCC = jOD.findOrientationOfPVCOrAluColumn(straightTiff, jCFS);
		
		//move column into the center of the canvas and cut out unnecessary parts of the canvas
		IJ.showStatus("Putting column straight and moving it to center of canvas, step 4 of 4...");
		double rim = 25;
		double mRadius = jMed.evaluate(newCC.outerMajorRadius);
		double toBeLeft = mRadius + rim;
		
		//check if cut-out Roi is not larger than the image..
		if (2 * toBeLeft > straightTiff.getWidth()) toBeLeft = straightTiff.getWidth() / 2;
		if (2 * toBeLeft > straightTiff.getHeight()) toBeLeft = straightTiff.getHeight() / 2;
			
		//do the cutting
		double diameter = 2 * toBeLeft;
		ImageStack outStack = new ImageStack((int)Math.round(diameter), (int)Math.round(diameter));
		
		//Boolean cutSuccessful = true;
		for (int i = 1 ; i < straightTiff.getNSlices() + 1 ; i++) {
			
			IJ.showStatus("Moving column to center of canvas " + i + "/" + straightTiff.getNSlices());
			
			straightTiff.setPosition(i);
			ImageProcessor nowIP = straightTiff.getProcessor();
			nowIP.setInterpolate(true);
			
			double XC = straightTiff.getWidth() / 2;
			double YC = straightTiff.getHeight() / 2;
			
			double dx = XC - jMed.evaluate(newCC.xmid);
			double dy = YC - jMed.evaluate(newCC.ymid);
			
			nowIP.translate(dx, dy);
			
			Roi cutRoi = new Roi(XC - toBeLeft, YC - toBeLeft, (int)Math.round(diameter), (int)Math.round(diameter));			
			nowIP.setRoi(cutRoi);
			ImageProcessor cutIP = nowIP.crop();
			
			outStack.addSlice(cutIP);
		}
		
		//free memory from straightTiff again..
		straightTiff.flush();
		IJ.freeMemory();IJ.freeMemory();
		
		outTiff.setStack(outStack);		
		
		//outTiff.draw();
		//outTiff.show();		
		
		return outTiff;
		
	}
	
	public ImagePlus cutColumnsEndsOff(ObjectDetector.ColumnContainer colCon) {
		
		int z;
		
		ImageStack outStack = new ImageStack(colCon.nowTiff.getWidth(), colCon.nowTiff.getHeight());
		ImagePlus outTiff = new ImagePlus();
	
		for (z = colCon.mendedCC.topOfColumn ; z < colCon.mendedCC.bottomOfColumn ; z++) { 
			
			IJ.showStatus("Finalizing straightened column " + (z - colCon.mendedCC.topOfColumn) + "/" + colCon.mendedCC.heightOfColumn);
			
			//set stack position to the correct depth
			colCon.nowTiff.setPosition(z);
			ImageProcessor localIP=colCon.nowTiff.getProcessor();			
			outStack.addSlice(localIP);
		}
		
		outTiff.setStack(outStack);	
		
		return outTiff;
	}

	/*public ImagePlus doSteelIlluCorrAndCutOffEnds(ImagePlus nowTiff, ObjectDetector.ColumnCoordinates preciseCC) {
		
		int z;
				
		double targetValue = preciseCC.steelGreyValue;
		double corrFactor;		
		
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImagePlus outTiff = new ImagePlus();
	
		for (z = preciseCC.topOfColumn ; z < preciseCC.bottomOfColumn + 1 ; z++) { 
			
			IJ.showStatus("Correcting illumination of slice " + (z - preciseCC.topOfColumn) + "/" + (preciseCC.bottomOfColumn - preciseCC.topOfColumn));
			
			//set stack position to the correct depth
			nowTiff.setPosition(z + 1);
			
			//calc corrfactor
			corrFactor = targetValue / preciseCC.wallGreyValues[z];
			
			//get and correct current slice
			ImageProcessor localIP = nowTiff.getProcessor();
			localIP.multiply(corrFactor);
			outStack.addSlice(localIP);
		}
		
		outTiff.setStack(outStack);		
		
		return outTiff;
	}*/
	
	///////////////////////////////////////////////
	// binarize
	///////////////////////////////////////////////
		
	public ImagePlus[] applyAChosenThresholdingMethod(ImagePlus nowTiff, String nowGaugePath, MenuWaiter.ThresholderMenuReturn mTMR, String myOutPath) {
		
		Median jMed = new Median();
		InputOutput jIO = new InputOutput();	
		HistogramStuff hist = new HistogramStuff();
		RoiHandler rH = new RoiHandler();
		ObjectDetector jOD = new ObjectDetector();
		
		ImagePlus[] outImg = {null, null, null};
		ImageStack myStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());

		ImagePlus airAndWaterImg = new ImagePlus();
						
		int[] myHist = new int[256 * 256];
		int[] newHist;
		int myThresh = 0;
		int i, j;
		
		AutoThresholder myAuto = new AutoThresholder();		
				
		//read gauge file
		ObjectDetector.ColumnCoordinates jCO = jIO.readGaugeFile(nowGaugePath);	
		PolygonRoi[] pRoi = rH.makeMeAPolygonRoiStack("inner", "exact", jCO, 0);
		PolygonRoi[] oRoi = rH.makeMeAPolygonRoiStack("inner", "wide", jCO, 0);
		
		//find illumination of this column	
		double[] illumination = jOD.findMedianWallGreyValues(nowTiff, jCO);
		double myIllumination = jMed.evaluate(illumination);	
		int[] onePercQuantile = new int[nowTiff.getNSlices()];
		
		//get the stacked histogram		
		for (i = 1 ; i < nowTiff.getNSlices() + 1 ; i++) {
			
			IJ.showStatus("Getting 16-bit histogram of slice #" + i + "/" + nowTiff.getNSlices());
			
			nowTiff.setPosition(i);		
			
			ImageProcessor myIP = nowTiff.getProcessor();
			ImageProcessor modIP = myIP.duplicate();
			
			//cut out everything outside column
			modIP.setRoi(pRoi[i - 1]);
			modIP.setColor(0);
			modIP.fillOutside(pRoi[i - 1]);	
			
			newHist=modIP.getHistogram();
			newHist[0] = 0; //set zero GV to zero
			for (j = 0 ; j < myHist.length ; j++) {
				myHist[j] = myHist[j] + newHist[j];
			}			
			
			//save 1% perc and wallOfGrey
			double[] cumHist = hist.calcCumulativeHistogram(newHist);			
			onePercQuantile[i-1] = hist.findPercentileFromCumHist(cumHist, 0.001);
	
		}		
		
		//save onePerc and wall of grey
		jIO.saveCriticalGreyValues(nowGaugePath, onePercQuantile, illumination);
			
		//prepare possible scaling nowTiff to 8-bit
		double[] cumHist = hist.calcCumulativeHistogram(myHist);
		double lBound = 0.0000001;double uBound = 0.99999;			
		int lowestVal = (int)Math.round(hist.findPercentileFromCumHist(cumHist, lBound));			//set lower grey value
		int largestVal = (int)Math.round(hist.findPercentileFromCumHist(cumHist, uBound));		
		float valSpan = (float)myIllumination - lowestVal;
		if (!mTMR.setMaxGrey2WallGrey | mTMR.useConstantThreshold == true) valSpan = largestVal - lowestVal;
			
		//find threshold	
		if (!mTMR.useConstantThreshold) {
		
			ImageStack eightBitStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
			
			for (i = 1 ; i < nowTiff.getNSlices() + 1 ; i++) {
				
				IJ.showStatus("Scaling 16-bit image to 8-bit at slice #" + i + "/" + nowTiff.getNSlices());
				
				nowTiff.setPosition(i);
				ImageProcessor myIP = nowTiff.getProcessor();
				ImageProcessor modIP = myIP.duplicate();
				modIP.max(largestVal);
				modIP.subtract(lowestVal);
				modIP.multiply(256 / valSpan);
				modIP.min(0);
				modIP.max(255);
				ImageProcessor eightIP = modIP.convertToByte(false);
				
				eightBitStack.addSlice(eightIP);
			}
			
			ImagePlus eightBitTiff = new ImagePlus();
			eightBitTiff.setStack(eightBitStack);
	
			//get the 8-bit stacked histogram
			int[] my8Hist = new int[256];
			int[] new8Hist;
			for (i = 1 ; i < eightBitTiff.getNSlices() + 1 ; i++) {
				
				IJ.showStatus("Getting 8-bit histogram of slice #" + i + "/" + (eightBitTiff.getNSlices()));
				
				eightBitTiff.setPosition(i);			
				
				ImageProcessor myIP = eightBitTiff.getProcessor();
				ImageProcessor modIP = myIP.duplicate();
				
				//cut out everything outside column
				modIP.setRoi(pRoi[i - 1]);
				modIP.setColor(0);
				modIP.fillOutside(pRoi[i - 1]);	
				
				new8Hist=modIP.getHistogram();	
				for (j = 0 ; j < my8Hist.length ; j++) {
					my8Hist[j] = my8Hist[j] + new8Hist[j];
				}			
				my8Hist[0] = 0; //set zero GV to zero
				my8Hist[255] = 0; //set last GV to zero	
			}				
		
			//find primary threshold
			myThresh = myAuto.getThreshold(mTMR.myPrimaryMethod, my8Hist);
		
			//find secondary threshold
			if(mTMR.mySecondaryMethod != null) {
				for(i = myThresh ; i < my8Hist.length ; i++) my8Hist[i] = 0;
				myThresh = myAuto.getThreshold(mTMR.mySecondaryMethod, my8Hist);
			}			
				
			//do binarization
			for (i = 1 ; i < eightBitTiff.getNSlices() + 1 ; i++) {
				
				IJ.showStatus("Binarizing slice #" + i + "/" + (eightBitTiff.getNSlices()));			
							
				eightBitTiff.setPosition(i);
				ImageProcessor myIP = eightBitTiff.getProcessor();
				
				//apply the threshold
				myIP.threshold(myThresh);
				
				//create a binary where the less dense phase is 255 and everything else 0
				myIP.invert();
				//for (int x = 0 ; x < myIP.getWidth() ; x++) {
				//	for (int y = 0 ; y < myIP.getHeight() ; y++) {						
				//		int thisImgPixel = myIP.getPixel(x, y);					
				//		myIP.putPixel(x, y, thisImgPixel);
				//	}
				//}
				
				//cut out everything outside column
				myIP.setRoi(pRoi[i - 1]);
				myIP.setColor(0);
				myIP.fillOutside(pRoi[i - 1]);	
			
				myStack.addSlice(myIP);
			}
			
			//create binary out image 
			outImg[0].setStack(myStack);
			
			//outImg[0].updateAndDraw();
			//outImg[0].show();
			
			//also save threshold		
			String thresholdSaverPath = myOutPath + "\\EightBitThresholds.txt";
			jIO.writeThreshold(thresholdSaverPath, nowGaugePath, myThresh);	
			
			//also save threshold comparison images	
			jIO.writeSnapshots4Comparison(myOutPath, nowGaugePath, nowTiff, outImg[0], 1, largestVal, lowestVal, valSpan, oRoi);
		
		}
		
		else { 	//in case that a constant threshold is used for all images..
	
			if (mTMR.airThreshold > 0) {
				
				ImageStack nowStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());		
				
				for (i = 1 ; i < nowTiff.getNSlices() + 1 ; i++) {
				
					IJ.showStatus("Binarizing slice for air-phase #" + i + "/" + nowTiff.getNSlices());
				
					nowTiff.setPosition(i);
					ImageProcessor myIP = nowTiff.getProcessor();
					ImageProcessor modIP = myIP.duplicate();
					
					modIP.threshold(mTMR.airThreshold);
					
					modIP.invert();
					
					modIP.setRoi(pRoi[i - 1]);
					modIP.setColor(0);
					modIP.fillOutside(pRoi[i - 1]);	
					
					ImageProcessor eightIP = modIP.convertToByte(false);
				
					nowStack.addSlice(eightIP);				
				}
				
				//create binary out image 
				ImagePlus airTiff = new ImagePlus();
				airTiff.setStack(nowStack);
				outImg[0] = airTiff;
				
				//airTiff.updateAndDraw();
				//airTiff.show();
				
			}
			
			if (mTMR.waterThreshold > 0) {
				
				ImageStack nowStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
				
				for (i = 1 ; i < nowTiff.getNSlices() + 1 ; i++) {
				
					IJ.showStatus("Binarizing slice for water-phase #" + i + "/" + nowTiff.getNSlices());
				
					nowTiff.setPosition(i);
					ImageProcessor myIP = nowTiff.getProcessor();
					ImageProcessor modIP = myIP.duplicate();
					
					modIP.threshold(mTMR.waterThreshold);
					
					modIP.invert();
					
					modIP.setRoi(pRoi[i - 1]);
					modIP.setColor(0);
					modIP.fillOutside(pRoi[i - 1]);	
					
					ImageProcessor eightIP = modIP.convertToByte(false);
				
					nowStack.addSlice(eightIP);				
				}
				
				//create binary out image 
				airAndWaterImg.setStack(nowStack);
				outImg[1] = airAndWaterImg;
				
				//subtract air from air and water to only receive water..
				ImageStack waterStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
				
				for (i = 1 ; i < nowTiff.getNSlices() + 1 ; i++) {
					
					IJ.showStatus("Separating water and fresh organic matter from air in slice #" + i + "/" + nowTiff.getNSlices());
				
					outImg[0].setPosition(i);
					ImageProcessor airIP = outImg[0].getProcessor();
					
					airAndWaterImg.setPosition(i);
					ImageProcessor waterIP = airAndWaterImg.getProcessor();
										
					ImageProcessor modIP = waterIP.duplicate();
					
					for (int x = 0 ; x < nowTiff.getWidth() ; x++) {
						for (int y = 0 ; y < nowTiff.getWidth() ; y++) {
							
							int airPx = airIP.getPixel(x, y);
							int waterPx = waterIP.getPixel(x, y);
							
							if (waterPx == 255 & airPx == 255) modIP.putPixel(x, y, 0);
							else modIP.putPixel(x, y, waterPx);
						}
					}					
					
					ImageProcessor eightIP = modIP.convertToByte(false);				
					
					waterStack.addSlice(eightIP);				
				}
				
				//create binary out image 
				ImagePlus waterTiff = new ImagePlus();
				waterTiff.setStack(waterStack);		
				outImg[1] = waterTiff;
				
			}
			
			if (mTMR.stoneThreshold > 0) {
				
				ImageStack nowStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
				
				for (i = 1 ; i < nowTiff.getNSlices() + 1 ; i++) {
				
					IJ.showStatus("Binarizing slice for stones #" + i + "/" + nowTiff.getNSlices());
				
					nowTiff.setPosition(i);
					ImageProcessor myIP = nowTiff.getProcessor();
					ImageProcessor modIP = myIP.duplicate();
					
					modIP.threshold(mTMR.stoneThreshold);
					
					modIP.setRoi(pRoi[i - 1]);
					modIP.setColor(0);
					modIP.fillOutside(pRoi[i - 1]);	
					
					ImageProcessor eightIP = modIP.convertToByte(false);
				
					nowStack.addSlice(eightIP);				
				}
				
				//create binary out image 
				ImagePlus stoneTiff = new ImagePlus();
				stoneTiff.setStack(nowStack);				
				outImg[2] = stoneTiff;
								
			}
			
			//also save threshold comparison images				
			jIO.writeSnapshots4ComparisonMega(myOutPath, nowGaugePath, nowTiff, outImg, 1, largestVal, lowestVal, valSpan, oRoi);
			
		}
					
		return outImg;
	}
	
	public ImagePlus binarizeInTwoSteps(ImagePlus myImg) {
		
		ImagePlus outImg = new ImagePlus();
						
		ImageStack myStack = new ImageStack(myImg.getWidth(), myImg.getHeight());
		
		ImageProcessor myIP = myImg.getProcessor().convertToByte(true);
				
		int[] myHist = new int[256];
		int[] newHist;
		int myThresh;
		int i, j;
		
		AutoThresholder myAuto = new AutoThresholder();
		String[] myMethods = AutoThresholder.getMethods();
		
		for (i = 1 ; i < myImg.getStackSize() + 1 ; i++) {
			myImg.setPosition(i);			
			myIP = myImg.getProcessor().convertToByte(true);
			newHist=myIP.getHistogram();	
			for (j = 0 ; j < myHist.length ; j++) {
				myHist[j] = myHist[j] + newHist[j];
			}			
		}
		
		myThresh = myAuto.getThreshold(myMethods[7], myHist);  //13 Renyi, 7 mean 
		
		for (i = myThresh; i < myHist.length ; i++) {
			myHist[i]=0;
		}
		
		int[] cutHist = new int[myThresh];
		for (i = 0 ; i < myThresh ; i++) {
			cutHist[i] = myHist[i];
		}
		
		myThresh = myAuto.getThreshold(myMethods[13], myHist);
				
		//do binarization
		for (i = 1 ; i < myImg.getStackSize() + 1 ; i++) {
			IJ.showStatus("Binarizing slice " + i + "/" + myImg.getNSlices());			
			myImg.setPosition(i);	
			myIP = myImg.getProcessor().convertToByte(true);
			myIP.threshold(myThresh);
			myIP.invert();
			myStack.addSlice(myIP);
		}
		
		//create binary outimage 
		outImg.setStack(myStack);
			
		return outImg;
	}
	
	public ImagePlus normalizeImages2AirAndWallGrey(ImagePlus nowTiff, String nowGaugePath, MenuWaiter.NormalizerReferences myNR, String myOutPath) {
		
		//init units
		Random rn = new Random();
		RoiHandler roi = new RoiHandler();
		HistogramStuff hist = new HistogramStuff();
		TailoredMaths maths = new TailoredMaths();		
		InputOutput jIO = new InputOutput();	

		//init vars
		ImagePlus outTiff = new ImagePlus();								
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());						
		int[] myHist = new int[256 * 256];
		int[] wallHist = new int[256 * 256];
		int[] newHist;
		int i, j;
				
		//read gauge file
		ObjectDetector.ColumnCoordinates jCO = jIO.readGaugeFile(nowGaugePath);		
		PolygonRoi[] pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 3);		
		PolygonRoi[] oRoi = roi.makeMeAPolygonRoiStack("outer", "manual", jCO, 3);		
		
		if (myNR.material.equalsIgnoreCase("Alu")) {
			pRoi = roi.makeMeAPolygonRoiStack("inner", "manual", jCO, 3);		
			oRoi = roi.makeMeAPolygonRoiStack("outerFromInner", "manual", jCO, 3);
		}
				
		//find illumination of this column	
		double[] lowerQuantile = new double[nowTiff.getNSlices()];
		double[] upperQuantile = new double[nowTiff.getNSlices()];
		double[] wall = new double[nowTiff.getNSlices()];
		double[] normLower = new double[nowTiff.getNSlices()];
	
		//get the stacked histogram		
		for (i = 1 ; i < nowTiff.getNSlices() + 1 ; i++) {
			
			IJ.showStatus("Getting 16-bit histogram of slice #" + i + "/" + nowTiff.getNSlices());
			
			nowTiff.setPosition(i);		
			
			ImageProcessor myIP = nowTiff.getProcessor();
			ImageProcessor modIP = myIP.duplicate();
			ImageProcessor wallIP = myIP.duplicate();
			
			//cut out everything outside column
			modIP.setRoi(pRoi[i - 1]);
			modIP.setColor(0);
			modIP.fillOutside(pRoi[i - 1]);	
			
			//get histogram of soil
			newHist=modIP.getHistogram();
			newHist[0] = 0; //set zero GV to zero
			for (j = 0 ; j < myHist.length ; j++) {
				myHist[j] = myHist[j] + newHist[j];
			}
			
			//find reference quantiles
			double[] cumHist = hist.calcCumulativeHistogram(newHist);			
			if (myNR.lowerReference > 0) lowerQuantile[i-1] = hist.findPercentileFromCumHist(cumHist, myNR.lowerReference);
			if (myNR.upperReference > 0) upperQuantile[i-1] = hist.findPercentileFromCumHist(cumHist, myNR.upperReference);
			
			//cut out everything but the wall
			wallIP.setRoi(pRoi[i-1]);
			wallIP.setColor(0);
			wallIP.fill(pRoi[i - 1]);
			
			wallIP.setRoi(oRoi[i - 1]);
			wallIP.setColor(0);
			wallIP.fillOutside(oRoi[i - 1]);
			
			//get histogram of wall
			newHist=wallIP.getHistogram();
			newHist[0] = 0; //set zero GV to zero
			for (j = 0 ; j < wallHist.length ; j++) {
				wallHist[j] = wallHist[j] + newHist[j];
			}	
			
			/*ImagePlus newImg = new ImagePlus();
			ImageStack newStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
			newStack.addSlice(wallIP);
			newImg.setStack(newStack);
			newImg.updateAndDraw();
			
			jOD.showMeMyRoi(newImg, wallIP, pRoi[i - 1], 1);
			jOD.showMeMyRoi(newImg, wallIP, oRoi[i - 1], 2);*/
				
			//find median wall grey			
			cumHist = hist.calcCumulativeHistogram(newHist);			
			wall[i-1] = hist.findPercentileFromCumHist(cumHist, 0.5);		
			
			//normalize lower reference
			if (myNR.lowerReference > 0) {
				if (myNR.upperReference > 0) normLower[i - 1] = lowerQuantile[i-1] / upperQuantile[i-1];
				else normLower[i - 1] = lowerQuantile[i-1] / wall[i-1];				
			}
			else {
				normLower[i - 1] = wall[i-1] / upperQuantile[i-1];					
			}
			  
		}		
		
		//write results into file
		int windowHalfSize = 100;  //window half-size for LOESS filter
		if (myNR.lowerReference > 0) myNR.originalLower = maths.LinearLOESSFilter(lowerQuantile, windowHalfSize);
		else myNR.originalLower = wall;
		if (myNR.upperReference > 0) myNR.originalUpper = maths.LinearLOESSFilter(upperQuantile, windowHalfSize);
		else myNR.originalUpper = wall;
			
		String myFileName = nowTiff.getTitle().substring(0, nowTiff.getTitle().length() - 4);
		jIO.writeIlluminationCorrectionIntoAsciiFile(myNR, myOutPath, myFileName);
		
		//apply correction image
		for (i = 0 ; i < nowTiff.getNSlices() ; i++) {
			
			IJ.showStatus("Normalizing slice #" + i + "/" + nowTiff.getNSlices());
			
			nowTiff.setPosition(i + 1);
			
			ImageProcessor myIP = nowTiff.getProcessor();
			
			double nowIntercept = myNR.lowerTarget;
			double nowSlope = (myNR.upperTarget - myNR.lowerTarget) / (myNR.originalUpper[i] - myNR.originalLower[i]);
						
			for (int x = 0 ; x < myIP.getWidth() ; x++) {
				for (int y = 0 ; y < myIP.getHeight() ; y++) {
					
					double nowX = myIP.getPixelValue(x, y) - myNR.originalLower[i];
					double newPixelValue = nowIntercept + nowSlope * nowX;
					
					//introduce a correction approach for filtering discretization artifact
					double fudger = nowSlope / (1 + Math.abs(nowX));
					double artifactSmoother = 2 * rn.nextDouble() * fudger - fudger;
					
					//put new value
					myIP.putPixelValue(x, y, newPixelValue + artifactSmoother);					
					
				}
			}
			
			outStack.addSlice(myIP);

		}
		
		outTiff.setStack(outStack);
		
		return outTiff;
		
	}
	
	
	///////////////////////////////////////////////
	// gimmicks
	//////////////////////////////////////////////
	
	public ImageProcessor invertSoilBinary(ImageProcessor myIP, PolygonRoi pRoi) {
	
		ImageProcessor modIP = myIP;
		
		modIP.invert();
		
		modIP.setRoi(pRoi);
		modIP.setColor(0);
		modIP.fillOutside(pRoi);
		modIP.resetRoi();
		
		return modIP;
				
	}
	
	public ImagePlus flipTiff(ImagePlus nowTiff) {
		
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImagePlus outTiff = new ImagePlus();
		
		int i;
		
		for (i = nowTiff.getNSlices() ; i > 0 ; i--) {
			
			nowTiff.setPosition(i);
			
			ImageProcessor myIP = nowTiff.getProcessor();
			
			outStack.addSlice(myIP);			
		}
		outTiff.setStack(outStack);
			
		return outTiff;
		
	}
	
	public ImagePlus cutOutRealThickness(ImagePlus nowTiff, ImagePlus rawThickImp) {
		
		ImagePlus outTiff = rawThickImp.duplicate();
		//ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		//rawThickImp.updateAndDraw();
		//rawThickImp.show();
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			nowTiff.setPosition(z+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			ImageProcessor outIP = outTiff.getProcessor();
			for (int x = 0 ; x < nowTiff.getWidth() ; x++) {		
				for (int y = 0 ; y < nowTiff.getHeight() ; y++) {
					int nowPix = nowIP.getPixel(x, y);
					if (nowPix == 0) outIP.putPixel(x, y, 0);  
				}
			}
		}
		
		//outTiff.updateAndDraw();
		//outTiff.show();
		
		return outTiff;
	}
	
	public ImagePlus binarizeGradientMask(ImagePlus nowTiff, int myThresh) {
		
		ImagePlus outTiff = new ImagePlus();
		ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());	
		
		//create mask for cutting out the rest..
		for (int i = 0 ; i < nowTiff.getNSlices() ; i++) {
			nowTiff.setPosition(i+1);		
		
			ImageProcessor binIP = nowTiff.getProcessor();
			binIP.threshold(myThresh);
			//binIP.dilate();			
			binIP.multiply(1/255.0);
			
			outStack.addSlice(binIP);
		}
		outTiff.setStack(outStack);
		
		return outTiff;		
	}

	public ImageProcessor fillHolesWithNeighboringValues(ImageProcessor surfaceIP, PolygonRoi[] pRoi, String topOrBottom) {
		
		JParticleCounter dPC = new JParticleCounter();
		AndAllTheRest aa = new AndAllTheRest();
		DisplayThings disp = new DisplayThings();
		
		//segment out deep roots			
		ImageProcessor justHoles = surfaceIP.duplicate();
		justHoles.threshold(1);
		justHoles.invert();
		
		justHoles.setBackgroundValue(0);
		if (topOrBottom.equalsIgnoreCase("top")) justHoles.fillOutside(pRoi[0]);
		else justHoles.fillOutside(pRoi[pRoi.length - 1]);
		ImageProcessor justHoles8Bit = justHoles.convertToByte(false);
		
		//identify the locations 
		ImagePlus holeOne = new ImagePlus("", justHoles8Bit);
		Object[] result = dPC.getParticles(holeOne, 1, 0d, Double.POSITIVE_INFINITY,-1, false);			
		int[][] particleLabels = (int[][]) result[1];
		long[] particleSizes = dPC.getParticleSizes(particleLabels);
		final int nParticles = particleSizes.length;			
		int[][] limits = dPC.getParticleLimits(holeOne, particleLabels, nParticles);  //xmin, xmax, ymin, ymax, zmin, zmax
		ImagePlus myHoles = disp.jDisplayParticleLabels(particleLabels, holeOne);
		IJ.freeMemory();IJ.freeMemory();
		
		//fill Holes
		ImageProcessor outIP = surfaceIP.duplicate();
		for (int i = 1 ; i < nParticles ; i++) {
			
			int xmin = limits[i][0] - 1;
			int xmax = limits[i][1] + 1;
			int ymin = limits[i][2] - 1;
			int ymax = limits[i][3] + 1;
							
			float[] X = new float[]{xmin, xmin, xmax, xmax, xmin};
			float[] Y = new float[]{ymin, ymax, ymax, ymin, ymin};
			
			PolygonRoi mRoi = new PolygonRoi(X, Y, Roi.POLYGON);
	
			//init list for hole-neighboring pixels
			ArrayList<Integer> neighbors = new ArrayList<Integer>();
			
			//cut out canvas around hole #i, remove all holes with other IDs and binarize the image.
			ImageProcessor hoIP = myHoles.getProcessor();
			
			ImageProcessor thIP = new ByteProcessor(xmax - xmin + 1, ymax - ymin + 1); 
			for (int x = xmin ; x < xmax + 1 ; x++) {
				for (int y = ymin ; y < ymax + 1 ; y++) {	
					int nowVal = (int)Math.round(hoIP.getPixelValue(x, y));
					if (nowVal == i) thIP.putPixel(x - xmin, y - ymin, 255); 		
					else thIP.putPixel(x, y, 0); 	
				}
			}
			
			// do the same for the surface elevation map
			ImageProcessor vlIP = surfaceIP.duplicate();
			vlIP.setRoi(mRoi);
			vlIP.crop();
			
			//make a copy of the binarized hole, dilate it and sample grey value of neighbors
			ImageProcessor cpIP = thIP.duplicate();
			cpIP.erode();
	
			for (int x = 0 ; x < cpIP.getWidth() ; x++) {
				for (int y = 0 ; y < cpIP.getHeight() ; y++) {						
					int nowTH = thIP.getPixel(x, y);
					int nowCP = cpIP.getPixel(x, y);
					int nowVL = (int)Math.round(vlIP.getPixelValue(x + xmin, y + ymin));
					if (nowTH == 0 & nowCP > 0 & nowVL > 0) {
						neighbors.add(nowVL);
					}
				}
			}
	
			//calculate median neighbor value
			double[] neighborsAsArray = new double[neighbors.size()]; 
			for (int j = 0 ; j < neighbors.size() ; j++) neighborsAsArray[j] = neighbors.get(j);
			int medianNeighbor = (int)StatUtils.percentile(neighborsAsArray, 10);
			
			//fill hole with neighborvalue
			for (int x = 0 ; x < cpIP.getWidth() ; x++) {
				for (int y = 0 ; y < cpIP.getHeight() ; y++) {
					int nowTH = thIP.getPixel(x, y);
					if (nowTH > 0) {
						outIP.putPixel(x + xmin, y + ymin, medianNeighbor);
					}
				}
			}				
		}
		
		return outIP;
	}
	
}
	


	


