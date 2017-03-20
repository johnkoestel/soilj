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

import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.AutoThresholder;
import ij.process.AutoThresholder.Method;

import org.apache.commons.math3.stat.StatUtils;

import SoilJ.tools.ObjectDetector.ApproximateColumnIllumination;

/**
 * MenuWaiter is a SoilJ class in which all pop-up menus are collected that are featured in SoilJ
 *
 * @author John Koestel
 *
 */

public class MenuWaiter implements PlugIn {


	public class BeamDeHardeningReturn {

		public boolean isSteelColumn = true;
		public int approximateRadius = 700;
		public int anglesChecked = 72;
		public double upperQ = 0.6;
		public double lowerQ = 0.8;
		public double correctionFunctionMedianFilter = 20;
		public String maskingMethod;

		double radiusOfMedianFilter = 2;
		double stepsize = 0.05;
		double cutoff = 0.2;
		int skip = 2;
		int maxEval = 200;
		double goodnessCriterion = 0.9;
		double severeGoodnessCriterion = 0.98;
		double gammaGoodnessCriterion = 0.7;
		double airPhaseClassMemberMinimum = 2;

	}

	public class ChiCalculatorReturn {

		public int wholeColumn = 0;
		public Boolean allClusters = false;
		public int clustersTakenIntoAccount = 101;
		public int minClusterSize = 5000;

	}


	public class MedianFilterAndUnsharpMaskReturn {

		public float medianFilterSizeXDir = 2;
		public float medianFilterSizeYDir = 2;
		public float medianFilterSizeZDir = 2;
		public double uMaskStandardDeviationXDir = 2;
		public double uMaskStandardDeviationYDir = 2;
		public double uMaskStandardDeviationZDir = 2;
		public float uMaskSharpeningWeight = 0.6f;

	}

	public class REVAnalyzerOptions {

		public String choiceOfRoi;
		public String choiceOfMethod;

		public int cubeX1;
		public int cubeX2;
		public int cubeY1;
		public int cubeY2;
		public int cubeZ1;
		public int cubeZ2;

		public int edgeX;
		public int edgeY;
		public int edgeZ;

		public int divNumber;
		public int stepNumber;
		public int stepLength;

		public int cylX;
		public int cylY;
		public int cylZ1;
		public int cylZ2;
		public int cylRadius;

		public boolean globVolume;
		public boolean globThickness;
		public boolean calcFractal;
		public boolean globAnisotropy;

		public boolean performParticleAnalyses;

		public boolean calcVolume;
		public boolean calcEuler;
		public boolean calcThickness;
		public boolean calcCriticalPoreDiameter;
		public boolean calcAnisotropy;
		public boolean calcPercolation;

		public boolean plotLabels;
		public boolean plotVolume;
		public boolean plotThickness;
		public boolean plotPercolation;

	}

	public class PoreSpaceAnalyzerOptions {

		public String choiceOfRoi;
		public String choiceOfZRoi;
		public String choiceOfXYRoi;

		public boolean cutZPercent; 
		public boolean cutXYPercent;		
		
		public int heightOfRoi;
		public int cutAwayFromTop;
		public int cutAwayFromBottom;
		
		public int cutAwayFromWall;
		public int cutAwayFromCenter;		
		
		public boolean includeSurfaceTopography;
		public boolean useInnerCircleFiles;
		public boolean cutCanvas = false;
		
		public double areaOfInterest;

		public int cubeX1;
		public int cubeX2;
		public int cubeY1;
		public int cubeY2;
		public int cubeZ1;
		public int cubeZ2;

		public int cylX;
		public int cylY;
		public int cylZ1;
		public int cylZ2;
		public int cylRadius;
		
		public boolean globVolume;
		public boolean globSurface;
		public boolean globThickness;
		public boolean calcCriticalPoreDiameter;
		public boolean calcChi;
		public boolean calcFractal;
		public boolean globAnisotropy;

		public boolean performParticleAnalyses;

		public boolean calcVolume;
		public boolean calcSurface;
		public boolean calcMoments;
		public boolean calcUnitVectors;
		public boolean calcEuler;
		public boolean calcThickness;
		public boolean calcCorrLength;
		public boolean calcAnisotropy;
		public boolean calcSkeleton;
		public boolean calcInclination;
		public boolean calcPercolation;

		public boolean plotLabels;
		public boolean plotVolume;
		public boolean plotSurface;
		public boolean plotThickness;
		public boolean plotPercolation;
		public boolean plotInclination;
		public boolean plotCorrelationLength;
		public boolean plotAnisotropy;

	}

	public class SurfaceFinderReturn {

		int neglectCrumbsOfLessThan;

	}

	public class ThresholderMenuReturn {

		public boolean useInnerCircle;

		public boolean filterImages = false;
		public String filterTag;
		public boolean setMaxgray2Wallgray;

		public Method myPrimaryMethod = AutoThresholder.Method.Default;
		public Method mySecondaryMethod = AutoThresholder.Method.Default;

		public int minThreshold = 0;
		public int maxThreshold = 0;

		public boolean useConstantThreshold;

		public boolean save3DImage;
		public boolean save4Evaluation;
		public boolean save4GeoDict;

	}

	public class OMFinderSettings {

		int mingrayValue;
		int maxgrayValue;
		double overlap;
		double windowSize;

	}

	public class RandomClusterGenerator {

		public final double pc =  0.0976d;

		public String shape;
		public int domainX;
		public int domainY;
		public int domainZ;
		public double[] porosityBounds;
		public double standardDeviation;
		public int numOfCopies;
		public String mode;
		public double[] porosityList;

	}

	public void run(String arg) {
				//ok, this is not needed..
	}

	public BeamDeHardeningReturn showBeamDeHardeningMenu() {

		GenericDialog gd = new GenericDialog("Beam de-hardening menu");

		BeamDeHardeningReturn mBDH = new BeamDeHardeningReturn();
		
		gd.addMessage("WARNING! This SoilJ module is still in an experimental stage and does therefore not offer a user-friendly interface.\n");
		gd.addMessage("Neither is the approach for fitting smooth function to radial intensity biases fully perfected nor is there\n");
		gd.addMessage("a tool to for a rapid checking of the quality and precision of the artifact removal, yet.\n");
		gd.addMessage("Admittedly, also a documentation of what how to tune the settings is still missing. \n");
		gd.addMessage("If you want to use this module, start with trying the default settings and proceed by trial-and-error if necessary.\n");
		gd.addMessage("If you are really desperate to use the SoilJ feature, try contacting john.koestel@slu.se and hope that he has time.\n");
		gd.addMessage("\n");

		//gd.addCheckbox("Do you have steel columns?", false);

		gd.addNumericField("What is the approximate radius of your columns in pixel? ", 550, 0, 6, "");

		//gd.addNumericField("How many radial profiles do you want to check per horizontal cross-section? ", 72, 0, 6, "");

		//gd.addNumericField("Please specify the lower percentile that you want to use for the beam-hardening correction", 60, 0, 6, "");

		//gd.addNumericField("Please specify the upper percentile that you want to use for the beam-hardening correction ", 80, 0, 6, "");

		//gd.addNumericField("Please specify the median filter radius applied to the correction function ", 20, 0, 6, "");

		gd.addNumericField("Specify the median filter radius applied before estimating the correction fucntions!", mBDH.radiusOfMedianFilter, 2, 6, "");

		gd.addNumericField("Please specify the radial discretization fraction you want to use for the beam de-hardening!", 0.05, 3, 6, "");

		gd.addNumericField("Please specify at which radial fraction that defines the centre of each horizontal crosssection!", 0.2, 3, 6, "");

		gd.addNumericField("Please specify how many radial steps you want to skip from the wall!", mBDH.skip, 0, 6, "");

		gd.addNumericField("Give a maximum of correction function evaluations during the fit!", mBDH.maxEval, 0, 6, "");

		gd.addNumericField("Give a minimum R^2 for accepting the correction function fit to be reliable!", mBDH.goodnessCriterion, 2, 6, "");

		gd.addNumericField("Give a minimum R^2 for accepting the correction function fit to be highly reliable!", mBDH.severeGoodnessCriterion, 2, 6, "");

		gd.addNumericField("Give a minimum R^2 for accepting the air-phase gamma correction function fit to be reliable!", mBDH.gammaGoodnessCriterion, 2, 6, "");

		gd.addNumericField("Give a number that the gray value sampled for the air phase has to be represneted with at least!", mBDH.airPhaseClassMemberMinimum, 0, 6, "");

		String[] choiceOfMasking = new String[3];
		choiceOfMasking[0] = "None";
		choiceOfMasking[1] = "Otsu";
		choiceOfMasking[2] = "MaxEntropy";
		gd.addChoice("Please choose a masking method for the non-matrix!", choiceOfMasking, "Otsu");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	mBDH.isSteelColumn = false; //gd.getNextBoolean();
	    	mBDH.approximateRadius = (int)Math.round(gd.getNextNumber());
	    	//mBDH.anglesChecked = (int)Math.round(gd.getNextNumber());
	    	//mBDH.lowerQ = gd.getNextNumber();
	    	//mBDH.upperQ = gd.getNextNumber();
	    	//mBDH.correctionFunctionMedianFilter = gd.getNextNumber();
	    	mBDH.radiusOfMedianFilter = gd.getNextNumber();
	    	mBDH.stepsize = gd.getNextNumber();
	    	mBDH.cutoff = gd.getNextNumber();
	    	mBDH.skip = (int)gd.getNextNumber();
	    	mBDH.maxEval = (int)gd.getNextNumber();
	    	mBDH.goodnessCriterion = gd.getNextNumber();
	    	mBDH.severeGoodnessCriterion = gd.getNextNumber();
	    	mBDH.gammaGoodnessCriterion = gd.getNextNumber();
	    	mBDH.airPhaseClassMemberMinimum = gd.getNextNumber();
	    	mBDH.maskingMethod = gd.getNextChoice();
	    }

		return mBDH;

	}

	public ChiCalculatorReturn showChiCalculatorMenu() {

		GenericDialog gd = new GenericDialog("chi calculator menu");

		ChiCalculatorReturn mCCR = new ChiCalculatorReturn();

		String[] choiceOfROI = new String[2];
		choiceOfROI[0] = "sub-sample";
		choiceOfROI[1] = "whole column";
		choiceOfROI[2] = "random field";
		gd.addChoice("Please choose a region of interest!", choiceOfROI, "sub-sample");

		gd.addCheckbox("Do you want to consider all pore-clusters in the sample?", false);

		String[] choiceOfCutoff = new String[3];
		choiceOfCutoff[0] = "11 largest pore-clusters";
		choiceOfCutoff[1] = "101 largest pore-clusters";
		choiceOfCutoff[2] = "1001 largest pore-clusters";
		gd.addChoice("If not, how many do you want to consider?", choiceOfCutoff, "101 largest pore-clusters");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	mCCR.wholeColumn = gd.getNextChoiceIndex();
	    	mCCR.allClusters = gd.getNextBoolean();
	    	int howManyIndex = gd.getNextChoiceIndex();
	    	switch (howManyIndex) {
	    		case 0: mCCR.clustersTakenIntoAccount = 11;break;
	    		case 1: mCCR.clustersTakenIntoAccount = 101;break;
	    		case 2: mCCR.clustersTakenIntoAccount = 1001;break;
	    	}
	    	if (mCCR.allClusters == true) mCCR.clustersTakenIntoAccount = 0;

	    }

		return mCCR;

	}

	public class ClipperMenuReturn {

		public boolean isSoilColumn = false;
		public boolean isCylinder = false;
		public boolean isRectangular = false;
		public boolean preserveOriginialCanvasSize = false;
		public boolean referenceIsSoilSurface = false;
		public boolean referenceIsTopMostSlice = false;
		public boolean addCanvasExccedance2Bottom = false;

		public int clipFromInnerPerimeter = 0;
		public int clipFromCanvasEdge = 0;
		public int canvasExceedsBy = 0;

		public int startAtSlice = 0;
		public int stopAtSlice = 0;
		public int heightOfROI = 0;

	}

	/*public ClipperMenuReturn showClipperDialog() {

		//construct objects
		GenericDialog gd = new GenericDialog("Clipper dialog");

		ClipperMenuReturn cMR = new ClipperMenuReturn();

		//build menu
		gd.addCheckbox("Do you want to use soil column outline coordinates saved as the Inner Circle?", true);

		String[] myChoices = new String[2];
		myChoices[0]="Cylindrical";myChoices[1]="Rectangular";
		gd.addRadioButtonGroup("What kind of shape do you want to cut out?", myChoices, 1, 2, myChoices[0]);

		gd.addCheckbox("Do you want to preserve the original canvas size?", false);

		gd.addNumericField("How many voxel do you want to cut away from the inner perimeter (canvas edge)?", 0, 0, 5, "");
		gd.addNumericField("By how many voxels should the canvas exceed the clipped image (if canvas size is not to be preserved)?", 2, 0, 3, "");
		gd.addCheckbox("Shall the blank canvas be also exceeded at the bottom of the column?", false);

		String[] myRefs = new String[2];
		myRefs[0]="The median soil surface as detected in SurfaceOfColumn";myRefs[1]="The topmost slice of the 3-D image";
		gd.addRadioButtonGroup("Choose a reference slice:", myRefs, 1, 2, myRefs[0]);

		gd.addNumericField("How many voxel below the reference depth is the upper boundary of your ROI?", 0, 0, 5, "");
		gd.addNumericField("How tall is your ROI in voxels?", 0, 0, 3, "");


		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	//get segmentation options
	    	cMR.isSoilColumn = gd.getNextBoolean();

	    	int shapeChoice = 0;
	    	String choice = gd.getNextRadioButton();
	    	for (int i = 0 ; i < myChoices.length ; i++) if (myChoices[i].equalsIgnoreCase(choice)) {
	    		shapeChoice = i;
	    		break;
	    	}
	    	switch (shapeChoice) {
				case 0 : {cMR.isCylinder = true; cMR.isRectangular = false; break;}
				case 1 : {cMR.isCylinder = false; cMR.isRectangular = true; break;}
	    	}

		    cMR.preserveOriginialCanvasSize = gd.getNextBoolean();

		    cMR.clipFromInnerPerimeter = (int)Math.round(gd.getNextNumber());
		    cMR.clipFromCanvasEdge = cMR.clipFromInnerPerimeter;
		    cMR.canvasExceedsBy = (int)Math.round(gd.getNextNumber());
		    cMR.addCanvasExccedance2Bottom = gd.getNextBoolean();

	    	int refChoice = 0;
	    	String rchoice = gd.getNextRadioButton();
	    	for (int i = 0 ; i < myRefs.length ; i++) if (myRefs[i].equalsIgnoreCase(rchoice)) {
	    		refChoice = i;
	    		break;
	    	}
	    	switch (refChoice) {
				case 0 : {cMR.referenceIsSoilSurface = true; cMR.referenceIsTopMostSlice = false; break;}
				case 1 : {cMR.referenceIsSoilSurface = false; cMR.referenceIsTopMostSlice = true; break;}
	    	}

		    cMR.startAtSlice = (int)Math.round(gd.getNextNumber());
		    cMR.heightOfROI = (int)Math.round(gd.getNextNumber());
		    cMR.stopAtSlice = cMR.startAtSlice + cMR.heightOfROI;

	    	return cMR;
	    }

	}*/

	public ThresholderMenuReturn showManualThresholdingDialog(ThresholderMenuReturn mTMR) {

		//construct objects
		GenericDialog gd = new GenericDialog("Thank you for choosing a user defined (constant) threshold for all images!");

		gd.addMessage("Make sure that you have run 'NormalizeTheseImages' when you apply this option!!!");

		gd.addNumericField("Enter the lower threshold for segmenting the image ", 0, 0, 5, "");

		gd.addNumericField("Enter the upper threshold for segmenting the image ", 0, 0, 5, "");
	
		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	mTMR.minThreshold = (int)Math.round(gd.getNextNumber());

	    	mTMR.maxThreshold = (int)Math.round(gd.getNextNumber());
	    	
	    	return mTMR;
	    }
	}

	public Boolean[] showFinalizerDialog() {

		//construct objects
		GenericDialog gd = new GenericDialog("");
		gd.addCheckbox("I want to merge two binaries", false);
		gd.addCheckbox("I want to cut out my sample using a binary mask", true);
		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	Boolean mergeOrNot = gd.getNextBoolean();
	    	Boolean cutOrNot = gd.getNextBoolean();
	    	Boolean[] finalizer = {mergeOrNot, cutOrNot};
	    	return finalizer;
	    }
	}

	public class Calc3DMenuReturn {

		public String operation;
		public String operationTag;
		public String filterTag;

		public boolean useInnerCircle;
		public boolean filterImages;

	}

	public Calc3DMenuReturn show3DCalcDialog(String A, String B) {

		//construct objects
		GenericDialog gd = new GenericDialog("3D Calculator dialog");
		Calc3DMenuReturn m3D = new Calc3DMenuReturn();

		//construct dialog window
		gd.addMessage(A);
		String[] items = {"Plus", "Minus"};

		gd.addChoice("", items, items[1]);

		gd.addMessage(B);

		gd.addMessage("");

		//add choices to menu
		gd.addCheckbox("Do you want to use the column outlines defined in 'inner circle' as a ROI?", false);

		gd.addMessage("Search-string for only segmenting some of the files in this folder. Leave blank for segmenting them all.");
		gd.addStringField("", "", 25);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	String myChoice = gd.getNextChoice();

	    	if (myChoice.equalsIgnoreCase(items[0])) {
	    		m3D.operation = "+";
	    		m3D.operationTag = "Plus";
	    	}

	    	if (myChoice.equalsIgnoreCase(items[1])) {
	    		m3D.operation = "-";
	    		m3D.operationTag = "Minus";
	    	}

	    	m3D.useInnerCircle = gd.getNextBoolean();

	    	m3D.filterTag = gd.getNextString();
	    	if (m3D.filterTag.hashCode() == 0) m3D.filterImages = false;
	    	else m3D.filterImages = true;

	      	return m3D;
	    }

	}

	public class HistogramMenuReturn {

		public boolean useInnerCircle;
		public boolean useSoilSurfaceFile;

		public boolean filterImages = false;
		public String filterTag;

		public boolean saveAsciis = true;
		public boolean saveImages = false;
		public boolean calcLumpedHistogram = false;

	}

	public HistogramMenuReturn showHistogramDialog() {

		//construct objects
		GenericDialog gd = new GenericDialog("HistoGrammar Dialog");

		HistogramMenuReturn hMR = new HistogramMenuReturn();

		//add choices to menu
		gd.addCheckbox("Do you want to use the column outlines defined in 'inner circle' as a ROI?", false);
		gd.addCheckbox("Do you want to use the soil surface topography?", false);

		/*gd.addMessage("Search-string for only segmenting some of the files in this folder. Leave blank for segmenting them all.");
		gd.addStringField("", "", 25);

		gd.addCheckbox("Save histograms as ASCII files?", true);

		gd.addCheckbox("Save histograms as images?", false);

		gd.addCheckbox("Calculate lumped histogram of all images in this folder?", true);
		gd.addMessage("(This option only makes sense if your images have normalized gray values, e.g. Hounsfield units)");
*/
		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	hMR.useInnerCircle = gd.getNextBoolean();
	    	hMR.useSoilSurfaceFile = gd.getNextBoolean();

	    	hMR.filterTag = gd.getNextString();
	    	if (hMR.filterTag.hashCode() == 0) hMR.filterImages = false;
	    	else hMR.filterImages = true;

	    	hMR.saveAsciis = true;//gd.getNextBoolean();
	    	hMR.saveImages = false;//gd.getNextBoolean();
	    	hMR.calcLumpedHistogram = true;//gd.getNextBoolean();

	      	return hMR;
	    }
	}

	public MedianFilterAndUnsharpMaskReturn showMedianAndUsharpMaskMenu() {

		GenericDialog gd = new GenericDialog("3-D median filter and 3-D unsharp mask");

		MedianFilterAndUnsharpMaskReturn mMU = new MedianFilterAndUnsharpMaskReturn();

		//set up menu
		gd.addNumericField("Median filter radius in X-direction ", mMU.medianFilterSizeXDir, 1, 6, "");
		gd.addNumericField("Median filter radius in Y-direction ", mMU.medianFilterSizeYDir, 1, 6, "");
		gd.addNumericField("Median filter radius in Z-direction ", mMU.medianFilterSizeZDir, 1, 6, "");
		gd.addNumericField("Standard deviation of the Gaussian kernel for the unsharp mask in X-direction ", mMU.uMaskStandardDeviationXDir, 1, 6, "");
		gd.addNumericField("Standard deviation of the Gaussian kernel for the unsharp mask in Y-direction ", mMU.uMaskStandardDeviationYDir, 1, 6, "");
		gd.addNumericField("Standard deviation of the Gaussian kernel for the unsharp mask in Z-direction ", mMU.uMaskStandardDeviationZDir, 1, 6, "");
		gd.addNumericField("Weight for unsharp mask application ", mMU.uMaskSharpeningWeight, 2, 6, "");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	mMU.medianFilterSizeXDir = (float)gd.getNextNumber();
	    	mMU.medianFilterSizeYDir = (float)gd.getNextNumber();
	    	mMU.medianFilterSizeZDir = (float)gd.getNextNumber();

	    	mMU.uMaskStandardDeviationXDir = gd.getNextNumber();
	    	mMU.uMaskStandardDeviationYDir = gd.getNextNumber();
	    	mMU.uMaskStandardDeviationZDir = gd.getNextNumber();
	    	mMU.uMaskSharpeningWeight = (float)gd.getNextNumber();
	    }

		return mMU;

	}

	public class RingArtifactRemoverOptions {

		public boolean useInnerCircle;
		public int incrementBetweenSlices2BeCheckedForRings;

	}

	public RingArtifactRemoverOptions returnRingArtifactRemoverOptions() {

		GenericDialog gd = new GenericDialog("3-D median filter and 3-D unsharp mask");

		RingArtifactRemoverOptions rARO = new RingArtifactRemoverOptions();

		// set up dialog box
		gd.addCheckbox("Do you want to use the InnerCircle files?", false);

		gd.addNumericField("Increment between horizontal to be checked for ring artifacts", 10, 0);


		//show dialog
		gd.showDialog();
		if (gd.wasCanceled()) return null;
		else {

			rARO.useInnerCircle = gd.getNextBoolean();

			rARO.incrementBetweenSlices2BeCheckedForRings = (int)Math.round(gd.getNextNumber());
		}

		return rARO;

	}

	public class CalibrationReferences {

		public boolean useInnerCircle;
		public boolean useSoilSurfaceFile;

		public String material;
		public String lowRef;
		public String hiRef;

		public double lowerTarget;
		public double upperTarget;

		public double lowerReference;
		public double upperReference;

		public double divisorMinimumMode;

		public boolean sampleLowerWithinSoil;
		public boolean sampleUpperWithinSoil;

		public double[] originalLower;
		public double[] originalUpper;

		public String lowerTag;
		public String upperTag;

	}

	public CalibrationReferences showCalibrationMenu() {

		MenuWaiter.CalibrationReferences mNR = new CalibrationReferences();

		//construct objects
		GenericDialog gd = new GenericDialog("Choose references for calibrating the gray-values!");

		// set up dialog box
		gd.addCheckbox("Do you want to use the InnerCircle files?", false);
	
		String[] columnMaterials = {"PVC", "aluminium"};
		gd.addRadioButtonGroup("What kind of material is your column made of?", columnMaterials, 2, 1, columnMaterials[1]);

		gd.addMessage("");

		gd.addNumericField("lower normalization target value", 5000, 0, 6, "");
		gd.addNumericField("upper normalization target value", 20000, 0, 6, "");

		gd.addMessage("");

		String[] lowRef = {"quantile", "wall"};
		gd.addRadioButtonGroup("How do you want to define you lower reference gray value?", lowRef, 2, 1, lowRef[0]);

		String[] hiRef = {"quantile", "wall"};
		gd.addRadioButtonGroup("How do you want to define you upper reference gray value?", hiRef, 2, 1, hiRef[1]);

		String[] samplingLocation = {"inside the column", "outside the column, close to the wall", "outside the column, with distance to the wall"};

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	mNR.useInnerCircle = gd.getNextBoolean();	    	

	    	mNR.material = gd.getNextRadioButton();

	    	mNR.lowerTarget = gd.getNextNumber();
	    	mNR.upperTarget = gd.getNextNumber();

	    	mNR.lowRef = gd.getNextRadioButton();
	    	mNR.hiRef = gd.getNextRadioButton();

	    }

	    if (mNR.lowRef.equalsIgnoreCase("quantile") | (mNR.hiRef.equalsIgnoreCase("quantile"))) {

	    	GenericDialog newGD = new GenericDialog("Please fill in missing input parameters!");

		    if (mNR.lowRef.equalsIgnoreCase("quantile")) {
		    	newGD.addNumericField("quantile for lower normalization reference value", 0.001, 5, 7, "");
		    	newGD.addRadioButtonGroup("where do you want to sample the lower reference value?", samplingLocation, 3, 1, samplingLocation[0]);
			}

		    if (mNR.hiRef.equalsIgnoreCase("quantile")) {
		    	newGD.addNumericField("quantile for upper normalization reference value", 0.001, 5, 7, "");
		    	newGD.addRadioButtonGroup("where do you want to sample the upper reference value?", samplingLocation, 3, 1, samplingLocation[0]);
			}

		    //show dialog
		    String sampleLowerWithinSoil = "";
		    String sampleUpperWithinSoil = "";

			myReference = "If you are using this plugin please cite the following references: \n\n";
			newGD.setInsets(40, 0, 0);newGD.addMessage(myReference);
			myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
			myReference += "Vadose Zone Journal, ????";
			newGD.setInsets(0, 0, 0);newGD.addMessage(myReference);

		    newGD.showDialog();

	  	    if (newGD.wasCanceled()) return null;
	  	    else {
	  	    	if (mNR.lowRef.equalsIgnoreCase("quantile")) {
	  	    		mNR.lowerReference = newGD.getNextNumber();
	  	    		sampleLowerWithinSoil = newGD.getNextRadioButton();

	  	    		String mylo = String.format("%1.4f", mNR.lowerReference);
	  		    	mNR.lowerTag = "Quantile" + mylo.substring(2, 5);

	  		    	if (sampleLowerWithinSoil == samplingLocation[0]) {
	  		    		mNR.sampleLowerWithinSoil = true;
	  		    		mNR.lowerTag += "Inside";
	  		    	}
	  		    	else {
	  		    		mNR.sampleLowerWithinSoil = false;
	  		    		if (sampleLowerWithinSoil == samplingLocation[1]) {	  		    	   		
	  		    	   		mNR.lowerTag += "Outside";
	  		    		}
	  		    		else {
	  		    			mNR.lowerTag += "FarOutside";
	  		    		}
	  		    	}
	  	    	}
  	    	
	  	    	if (mNR.hiRef.equalsIgnoreCase("quantile")) {
	  	    		mNR.upperReference = newGD.getNextNumber();
	  	    		sampleUpperWithinSoil = newGD.getNextRadioButton();

	  	    		String myup = String.format("%1.4f", mNR.upperReference);
	  		    	mNR.lowerTag = "Quantile" + myup.substring(2, 5);

	  		    	if (sampleUpperWithinSoil == samplingLocation[0]) {
	  		    		mNR.sampleUpperWithinSoil = true;
	  		    		mNR.upperTag += "Inside";
	  		    	}
	  		    	else {
	  		    		mNR.sampleUpperWithinSoil = false;
	  		    		if (sampleUpperWithinSoil == samplingLocation[1]) {	  		    	   		
	  		    	   		mNR.upperTag += "Outside";
	  		    		}
	  		    		else {
	  		    			mNR.upperTag += "FarOutside";
	  		    		}
	  		    	}
	  	    	}
	  	    }
  	    }

  	    if (mNR.lowRef.equalsIgnoreCase("wall")) mNR.lowerTag = "Wall";

  	    if (mNR.hiRef.equalsIgnoreCase("wall")) mNR.upperTag = "Wall";
  
	    return mNR;

	}

	public OMFinderSettings showOMFinderMenu() {

		OMFinderSettings oMF = new OMFinderSettings();

		GenericDialog gd = new GenericDialog("Fresh organic matter extractor dialog");

		gd.addNumericField("What is the minimum gray value of the organic matter? ", 9500, 0, 6, "");
		gd.addNumericField("What is the maximum gray value of the organic matter?", 15000, 0, 6, "");
		//gd.addNumericField("Window size relative to span between min and max gray values?", 0.33, 2, 4, "");
		//gd.addNumericField("Please choose an overlap (%)", 0.5, 2, 4, "");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	oMF.mingrayValue = (int)Math.round(gd.getNextNumber());
	    	oMF.maxgrayValue = (int)Math.round(gd.getNextNumber());
	    	oMF.windowSize = 0.33f;//gd.getNextNumber();
	    	oMF.overlap = 0.5f;//gd.getNextNumber();

			return oMF;
	    }

	}

	public PoreSpaceAnalyzerOptions showPoreSpaceAnalyzerMenu() {

		//construct objects
		GenericDialog gd = new GenericDialog("Please tell me what you want to get analysed!");

		PoreSpaceAnalyzerOptions mPSAO = new PoreSpaceAnalyzerOptions();
	    mPSAO.useInnerCircleFiles = false;
	    mPSAO.includeSurfaceTopography = false;

		String[] choiceOfRoi = new String[6];
		choiceOfRoi[0] = "No ROI but the entire image!";
		choiceOfRoi[1] = "Cylindrical column with outlines defined under 'Inner Circle'";
		choiceOfRoi[2] = "Cylindrical ROI";
		choiceOfRoi[3] = "Cuboid ROI";
		choiceOfRoi[4] = "No ROI but the top half of the entire image!";
		choiceOfRoi[5] = "No ROI but the bottom half of the entire image!";
		gd.addRadioButtonGroup("Please choose a region of interest!", choiceOfRoi, 4, 1, "No ROI but the entire image!");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Doube M, Kłosowski MM, Arganda-Carreras I, Cordeliéres F, Dougherty RP, Jackson J, Schmid B, Hutchinson JR, Shefelbine SJ. (2010)\n";
		myReference += "BoneJ: free and extensible bone image analysis in ImageJ. Bone 47:1076-9. doi: 10.1016/j.bone.2010.08.023\n\n";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	String myChoice = gd.getNextRadioButton();
	    	int myChoiceIndex = 0;
	    	for (int i = 0 ; i < choiceOfRoi.length ; i++) if (myChoice.equalsIgnoreCase(choiceOfRoi[i])) {
	    		myChoiceIndex = i;
	    		break;
	    	}	    	
	    	switch (myChoiceIndex) {
	    		case 0: mPSAO.choiceOfRoi = "Everything!"; break;
	    		case 1: mPSAO.choiceOfRoi = "RealSample"; break;
	    		case 2: mPSAO.choiceOfRoi = "Cylinder"; break;
	    		case 3: mPSAO.choiceOfRoi = "Cuboid"; break;
	    		case 4: mPSAO.choiceOfRoi = "TopOfEveryThing"; break;
	    		case 5: mPSAO.choiceOfRoi = "BottomOfEveryThing"; break;
	    	}
	    }

	    gd.removeAll();

	    if (mPSAO.choiceOfRoi.equalsIgnoreCase("Everything!") | mPSAO.choiceOfRoi.equalsIgnoreCase("TopOfEveryThing") | mPSAO.choiceOfRoi.equalsIgnoreCase("BottomOfEveryThing")) {

	    	gd.addNumericField("What is the cross-sectional area of the object of interest in pixels? ", 0, 0, 8, "");
	    	gd.addMessage("(This is needed in case you want SoilJ to calculate the porosity but your object of interest does not cover the entire XY-plane, e.g. in the case of a cylindrical object)");
	    	gd.addMessage("");
			gd.addNumericField("How many slices do you want to cut away from the top? ", 0, 0, 6, "");
			gd.addNumericField("How many slices do you want to cut away from the bottom? ", 0, 0, 6, "");

			myReference = "If you are using this plugin please cite the following references: \n\n";
			gd.setInsets(50, 0, 0);gd.addMessage(myReference);
			myReference = "Doube M, Kłosowski MM, Arganda-Carreras I, Cordeliéres F, Dougherty RP, Jackson J, Schmid B, Hutchinson JR, Shefelbine SJ. (2010)\n";
			myReference += "BoneJ: free and extensible bone image analysis in ImageJ. Bone 47:1076-9. doi: 10.1016/j.bone.2010.08.023\n\n";
			gd.setInsets(0, 0, 0);gd.addMessage(myReference);
			myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
			myReference += "Vadose Zone Journal, ????";
			gd.setInsets(0, 0, 0);gd.addMessage(myReference);

	    	//show dialog
			gd.showDialog();
		    if (gd.wasCanceled()) return null;
		    else {
		    	//domain definitions
		    	mPSAO.areaOfInterest = gd.getNextNumber();
		    	mPSAO.cutAwayFromTop = (int)Math.round(gd.getNextNumber());
		    	mPSAO.cutAwayFromBottom = (int)Math.round(gd.getNextNumber());
		    }
	    }

		if (mPSAO.choiceOfRoi.equalsIgnoreCase("RealSample")) {
			
			GenericDialog gdR = new GenericDialog("Please tell me more");

			
			mPSAO.useInnerCircleFiles = true;
			
			String[] choiceOfSubRoi = new String[5];
			choiceOfSubRoi[0] = "The entire soil column";
			choiceOfSubRoi[1] = "The central part of the soil column (define depth under soil surface and column height in voxel)";
			choiceOfSubRoi[2] = "The central part of the soil column (define number of voxels to be cut away from top and/or bottom)";
			choiceOfSubRoi[3] = "The central part of the soil column (define percentage to be cut away from top and/or bottom)";
			choiceOfSubRoi[4] = "The top or bottom half of the soil column";
			gdR.addRadioButtonGroup("Please choose a region of interest along the Z-axis", choiceOfSubRoi, 5, 1, "The entire soil column");

			String[] choiceOfSubRoi2 = new String[5];
			choiceOfSubRoi2[0] = "The entire soil column";
			choiceOfSubRoi2[1] = "The central part of the soil column (define number of voxels to be cut away from the wall)";
			choiceOfSubRoi2[2] = "The central part of the soil column (define percentage to be cut away from the wall)";
			choiceOfSubRoi2[3] = "A hollow cylinder cut out from the soil column (define number of voxels to be cut away from the center and/or the wall)";
			choiceOfSubRoi2[4] = "A hollow cylinder cut out from the soil column (define percentage to be cut away from the center and/or the wall)";
			gdR.addRadioButtonGroup("Please choose a region of interest in the XY-plane", choiceOfSubRoi2, 5, 1, "The entire soil column");
			
			String[] referencePoint = new String[2];
			referencePoint[0] = "The top and bottom soil surfaces as found by the SoilJ SoilSurfaceFinder";
			referencePoint[1] = "The top and bottom of the 3-D image canvas";
			gdR.addRadioButtonGroup("Please choose reference points for the top and bottom of your column", referencePoint, 2, 1, "The top and bottom soil surfaces as found by the SoilJ SoilSurfaceFinder");
		
			myReference = "If you are using this plugin please cite the following references: \n\n";
			gdR.setInsets(50, 0, 0);gdR.addMessage(myReference);
			myReference = "Doube M, Kłosowski MM, Arganda-Carreras I, Cordeliéres F, Dougherty RP, Jackson J, Schmid B, Hutchinson JR, Shefelbine SJ. (2010)\n";
			myReference += "BoneJ: free and extensible bone image analysis in ImageJ. Bone 47:1076-9. doi: 10.1016/j.bone.2010.08.023\n\n";
			gdR.setInsets(0, 0, 0);gdR.addMessage(myReference);
			myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
			myReference += "Vadose Zone Journal, ????";
			gdR.setInsets(0, 0, 0);gdR.addMessage(myReference);

			//show dialog
			gdR.showDialog();
		    if (gdR.wasCanceled()) return null;
		    else {
		    	
    			mPSAO.cutZPercent = false;
    			mPSAO.cutXYPercent = false;
    			
    			mPSAO.cutAwayFromBottom = 0;
    			mPSAO.cutAwayFromTop = 0;
    			mPSAO.cutAwayFromCenter = 0;
    			mPSAO.cutAwayFromWall = 0;
    			mPSAO.heightOfRoi = 0;
		    	
		    	//domain definitions Z
		    	String myChoice = gdR.getNextRadioButton();
		    	int myChoiceIndex = 0;
		    	for (int i = 0 ; i < choiceOfSubRoi.length ; i++) if (myChoice.equalsIgnoreCase(choiceOfSubRoi[i])) {
		    		myChoiceIndex = i;
		    		break;
		    	}
		    
		    	switch (myChoiceIndex) {
		    		case 0: {
		    			mPSAO.choiceOfZRoi = "everything";		    			
		    			break;
		    		}
		    		case 1: {
		    			
		    			mPSAO.choiceOfZRoi = "fixedHeight";		    			
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			gd3.addNumericField("Cut away from top", 0, 0, 5, "voxel");
		    			gd3.addNumericField("Column height", 0, 0, 5, "voxel");
		    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	mPSAO.cutAwayFromTop = (int) Math.round(gd3.getNextNumber());
		    		    	mPSAO.heightOfRoi = (int) Math.round(gd3.getNextNumber());
		    		    }
		    		
		    			break;
		    			
		    		}
		    		case 2: {
		    		
		    			mPSAO.choiceOfZRoi= "relativeHeightVoxels"; 
		    		    			
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			gd3.addNumericField("Cut away from top", 0, 0, 5, "voxel");
		    			gd3.addNumericField("Cut away from bottom", 0, 0, 5, "voxel");
		    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	mPSAO.cutAwayFromTop = (int) Math.round(gd3.getNextNumber());
		    		    	mPSAO.cutAwayFromBottom = (int) Math.round(gd3.getNextNumber());
		    		    }	    		    
		    			
		    			break;
		    		}
		    		case 3: {
		    			mPSAO.choiceOfZRoi = "relativeHeightPercent";
			    	
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			gd3.addNumericField("Cut away from top", 0, 0, 5, "%");
		    			gd3.addNumericField("Cut away from bottom", 0, 0, 5, "%");
		    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	mPSAO.cutZPercent = true;
		    		    	mPSAO.cutAwayFromTop = (int) Math.round(gd3.getNextNumber());
		    		    	mPSAO.cutAwayFromBottom = (int) Math.round(gd3.getNextNumber());
		    		    }
		    			break;
		    		}
		    		case 4: {
		    			mPSAO.choiceOfZRoi = "TopOrBottom";
		    			
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			
		    			String[] topBot= new String[2];
		    			topBot[0] = "Top half of column";
		    			topBot[1] = "Bottom half of column";
		    			gd3.addRadioButtonGroup("Please choose your preference", topBot, 2, 1, "Top half of column");
		    					    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	myChoice = gd.getNextRadioButton();
		    		    	myChoiceIndex = 0;
		    		    	for (int i = 0 ; i < topBot.length ; i++) if (myChoice.equalsIgnoreCase(topBot[i])) {
		    		    		myChoiceIndex = i;
		    		    		break;
		    		    	}
		    		    
		    		    	switch (myChoiceIndex) {
		    		    		case 0: {
		    		    			mPSAO.cutZPercent = true;
				    		    	mPSAO.cutAwayFromTop = 0;
				    		    	mPSAO.cutAwayFromBottom = 50;		    		    		
		    		    			break;
		    		    		}
		    		    		case 1: {
		    		    			mPSAO.cutZPercent = true;
				    		    	mPSAO.cutAwayFromTop = 50;
				    		    	mPSAO.cutAwayFromBottom = 0;		    		    		
		    		    			break;
		    		    		}
		    		    	}
		    		    }
		    			break;
		    		}
		    	}
		    	
		    	myChoice = gdR.getNextRadioButton();
		    	myChoiceIndex = 0;
		    	for (int i = 0 ; i < choiceOfSubRoi2.length ; i++) if (myChoice.equalsIgnoreCase(choiceOfSubRoi2[i])) {
		    		myChoiceIndex = i;
		    		break;
		    	}
		    
		    	switch (myChoiceIndex) {
		    		case 0: {
		    			mPSAO.choiceOfXYRoi = "everything"; break;
		    		}
		    		case 1: {
		    			mPSAO.choiceOfXYRoi = "centralVoxels"; 
		    			
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			gd3.addNumericField("Cut away from wall", 0, 0, 5, "voxel");
		    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	mPSAO.cutAwayFromWall = (int) Math.round(gd3.getNextNumber());		    		    
		    		    }
		    			
		    			break;
		    		}
		    		case 2: {		    			
		    			
		    			mPSAO.choiceOfXYRoi= "centralPercent";	
		    			mPSAO.cutXYPercent = true;
		    			
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			gd3.addNumericField("Cut away from wall", 0, 0, 5, "%");
		    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	mPSAO.cutAwayFromWall = (int) Math.round(gd3.getNextNumber());		   
		    		    }
		    			
		    			break;
		    		}
		    		case 3: {
		    			mPSAO.choiceOfXYRoi = "donutVoxels"; 
		    			
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			gd3.addNumericField("Cut away from wall", 0, 0, 5, "voxel");
		    			gd3.addNumericField("Cut away from center", 0, 0, 5, "voxel");	
		    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	mPSAO.cutAwayFromWall = (int) Math.round(gd3.getNextNumber());
		    		    	mPSAO.cutAwayFromCenter = (int) Math.round(gd3.getNextNumber());
		    		    }
		    			
		    			break;
		    		}
		    		case 4: {
		    			mPSAO.choiceOfXYRoi = "donutPercent"; 
		    			
		    			mPSAO.cutXYPercent = true;
		    			
		    			GenericDialog gd3 = new GenericDialog("Please tell me more!");
		    			gd3.addNumericField("Cut away from wall", 0, 0, 5, "%");
		    			gd3.addNumericField("Cut away from center", 0, 0, 5, "%");			    		
		    			
		    			gd3.showDialog();
		    		    if (gd3.wasCanceled()) return null;
		    		    else {
		    		    	mPSAO.cutAwayFromWall = (int) Math.round(gd3.getNextNumber());
		    		    	mPSAO.cutAwayFromCenter = (int) Math.round(gd3.getNextNumber());
		    		    }
		    			
		    			break;		    		
		    		}
		    	}
		    	
		    	myChoice = gdR.getNextRadioButton();
		    	myChoiceIndex = 0;
		    	for (int i = 0 ; i < referencePoint.length ; i++) if (myChoice.equalsIgnoreCase(referencePoint[i])) {
		    		myChoiceIndex = i;
		    		break;
		    	}
		    
		    	switch (myChoiceIndex) {
		    		case 0: mPSAO.includeSurfaceTopography = true;break;
		    		case 1: mPSAO.includeSurfaceTopography = false;break;		    	    		
		    	}		
		    }
		    
		}

	    if (mPSAO.choiceOfRoi.equalsIgnoreCase("Cylinder")) {
			gd.addNumericField("X-coordinate of cylinder center ", 0, 0, 6, "");
			gd.addNumericField("Y-coordinate of cylinder center ", 0, 0, 6, "");
			gd.addNumericField("Z-coordinate of cylinder top ", 0, 0, 6, "");
			gd.addNumericField("Z-coordinate of cylinder bottom ", 0, 0, 6, "");
			gd.addNumericField("Radius of cylinder ", 0, 0, 6, "");

			myReference = "If you are using this plugin please cite the following references: \n\n";
			gd.setInsets(50, 0, 0);gd.addMessage(myReference);
			myReference = "Doube M, Kłosowski MM, Arganda-Carreras I, Cordeliéres F, Dougherty RP, Jackson J, Schmid B, Hutchinson JR, Shefelbine SJ. (2010)\n";
			myReference += "BoneJ: free and extensible bone image analysis in ImageJ. Bone 47:1076-9. doi: 10.1016/j.bone.2010.08.023\n\n";
			gd.setInsets(0, 0, 0);gd.addMessage(myReference);
			myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
			myReference += "Vadose Zone Journal, ????";
			gd.setInsets(0, 0, 0);gd.addMessage(myReference);

			//show dialog
			gd.showDialog();
		    if (gd.wasCanceled()) return null;
		    else {
		    	//domain definitions
		    	mPSAO.cylX = (int)Math.round(gd.getNextNumber());
		    	mPSAO.cylY = (int)Math.round(gd.getNextNumber());
		    	mPSAO.cylZ1 = (int)Math.round(gd.getNextNumber());
		    	mPSAO.cylZ2 = (int)Math.round(gd.getNextNumber());
		    	mPSAO.cylRadius	 = (int)Math.round(gd.getNextNumber());
		    }
		}

		if (mPSAO.choiceOfRoi.equalsIgnoreCase("Cuboid")) {
			gd.addNumericField("Left X-coordinate of cuboid ", 0, 0, 6, "");
			gd.addNumericField("Right X-coordinate of cuboid (0 corresponds to the maximum x)", 0, 0, 6, "");
			gd.addNumericField("Topmost Y-coordinate of cuboid ", 0, 0, 6, "");
			gd.addNumericField("Bottommost Y-coordinate of cuboid (0 corresponds to the maximum y)", 0, 0, 6, "");
			gd.addNumericField("Z-coordinate of cuboid top ", 0, 0, 6, "");
			gd.addNumericField("Z-coordinate of cuboid bottom ", 0, 0, 6, "");

			//show dialog
			gd.showDialog();
		    if (gd.wasCanceled()) return null;
		    else {
		    	//domain definitions
		    	mPSAO.cubeX1 = (int)Math.round(gd.getNextNumber());
		    	mPSAO.cubeX2 = (int)Math.round(gd.getNextNumber());
		    	mPSAO.cubeY1 = (int)Math.round(gd.getNextNumber());
		    	mPSAO.cubeY2 = (int)Math.round(gd.getNextNumber());
		    	mPSAO.cubeZ1 = (int)Math.round(gd.getNextNumber());
		    	mPSAO.cubeZ2 = (int)Math.round(gd.getNextNumber());
		    }
		}

		gd.removeAll();

		GenericDialog gd2 = new GenericDialog("Please tell me what you want to get analysed!");

		String[] whatIntegralMeasuresShallIInvestigate = new String[]{"Macropore Volume","Average Thickness","Critical Pore Diameter","Fractal Dimension","Anisotropy"};
		boolean[] myIntChoices = new boolean[]{true,true,true,true,true};
		gd2.setInsets(20, 100, 0);gd2.addMessage("\nWhich global morphological measures should I calculate?");
		gd2.setInsets(0, 100, 0);gd2.addCheckboxGroup(3, 3, whatIntegralMeasuresShallIInvestigate, myIntChoices);

		String[] whatShallIInvestigate = new String[]{"Volume", "Unit Vectors","Euler Number","Thickness",
				"Anisotropy","Percolation","Crit. Pore Diameter"};
		boolean[] myChoices = new boolean[]{true,false,true,true,true,true,true};
		gd2.setInsets(20, 100, 0);gd2.addMessage("\nWhich morphological measures should I calculate for each cluster?");
		gd2.setInsets(0, 100, 0);gd2.addCheckboxGroup(4, 3, whatShallIInvestigate, myChoices);

		String[] whichImagesShallIPlot = new String[]{"Cluster Label","Volume","Thickness","Percolating Clusters"};
		boolean[] myPlotChoices = new boolean[]{true,false,true,false};
		gd2.setInsets(20, 100, 0);gd2.addMessage("\nWhich images shall I save?");
		gd2.setInsets(0, 100, 0);gd2.addCheckboxGroup(3, 3, whichImagesShallIPlot, myPlotChoices);

		myReference = "If you are using this plugin please cite the following references: \n\n";
		gd2.setInsets(50, 0, 0);gd2.addMessage(myReference);
		myReference = "Doube M, Kłosowski MM, Arganda-Carreras I, Cordeliéres F, Dougherty RP, Jackson J, Schmid B, Hutchinson JR, Shefelbine SJ. (2010)\n";
		myReference += "BoneJ: free and extensible bone image analysis in ImageJ. Bone 47:1076-9. doi: 10.1016/j.bone.2010.08.023\n\n";
		gd2.setInsets(0, 0, 0);gd2.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd2.setInsets(0, 0, 0);gd2.addMessage(myReference);

		//show dialog
		gd2.showDialog();
	    if (gd2.wasCanceled()) return null;	    
	    else {

	    	mPSAO.calcCriticalPoreDiameter = false;
	    	
	    	//global measures
	    	mPSAO.globVolume = gd2.getNextBoolean();
			mPSAO.globThickness = gd2.getNextBoolean();
			mPSAO.calcCriticalPoreDiameter = gd2.getNextBoolean();
			mPSAO.calcFractal = gd2.getNextBoolean();
			mPSAO.globAnisotropy = gd2.getNextBoolean();

			//local measures
	    	mPSAO.calcVolume = gd2.getNextBoolean();
			mPSAO.calcUnitVectors = gd2.getNextBoolean();
			mPSAO.calcEuler = gd2.getNextBoolean();
			mPSAO.calcThickness = gd2.getNextBoolean();
			mPSAO.calcAnisotropy = gd2.getNextBoolean();
			mPSAO.calcPercolation = gd2.getNextBoolean();
			mPSAO.calcCriticalPoreDiameter = gd2.getNextBoolean();

			//plotting options
			mPSAO.plotLabels = gd2.getNextBoolean();
			mPSAO.plotVolume = gd2.getNextBoolean();
			mPSAO.plotThickness = gd2.getNextBoolean();
			mPSAO.plotPercolation = gd2.getNextBoolean();
			//mPSAO.plotAnisotropy = gd2.getNextBoolean();

			mPSAO.performParticleAnalyses = false;
			if (mPSAO.calcVolume == true ||
					mPSAO.calcUnitVectors == true ||
					mPSAO.calcEuler == true ||
					mPSAO.calcThickness == true ||
					mPSAO.calcAnisotropy == true ||
					mPSAO.calcPercolation == true ||
					mPSAO.calcCriticalPoreDiameter == true ||
					mPSAO.calcCriticalPoreDiameter == true ||
				
					//plotting options
					mPSAO.plotLabels == true ||
					mPSAO.plotVolume == true ||
					mPSAO.plotThickness == true ||
					mPSAO.plotPercolation == true ||
					mPSAO.plotAnisotropy == true) {

				mPSAO.performParticleAnalyses = true;

			}

	    }

	    return mPSAO;

	}

	public REVAnalyzerOptions showREVAnalyzerMenu() {

		///////////////////////////////////////
		// type of ROI and Method
		///////////////////////////////////////

		//construct objects
		GenericDialog gd = new GenericDialog("Please tell me what you want to get analysed!");

		REVAnalyzerOptions mRA = new REVAnalyzerOptions();

		String[] choiceOfRoi = new String[2];
		choiceOfRoi[0] = "Cuboid ROI";
		choiceOfRoi[1] = "Cylindrical ROI (defunct)";
		gd.addRadioButtonGroup("Please choose a region of interest!", choiceOfRoi, 1, 2, "Cuboid ROI");

		String[] choiceOfMethod = new String[3];
		choiceOfMethod[0] = "Sub-ROI by division";
		choiceOfMethod[1] = "Sub-ROI by shrinkage";
		choiceOfMethod[2] = "Sub-ROI by division and shrinkage";
		gd.addRadioButtonGroup("Please choose a region of interest!", choiceOfMethod, 3, 1, "Sub-ROI by division");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	//get ROI typ
	    	String myChoice = gd.getNextRadioButton();
	    	int myChoiceIndex = 0;
	    	for (int i = 0 ; i < choiceOfRoi.length ; i++) if (myChoice.equalsIgnoreCase(choiceOfRoi[i])) {
	    		myChoiceIndex = i;
	    		break;
	    	}
	    	switch (myChoiceIndex) {
	    		case 0: mRA.choiceOfRoi = "Cuboid"; break;
	    		//case 1: mRA.choiceOfRoi = "Cylinder"; break;
	    	}

	    	//get Method
	    	myChoice = gd.getNextRadioButton();
	    	myChoiceIndex = 0;
	    	for (int i = 0 ; i < choiceOfMethod.length ; i++) if (myChoice.equalsIgnoreCase(choiceOfMethod[i])) {
	    		myChoiceIndex = i;
	    		break;
	    	}
	    	switch (myChoiceIndex) {
	    		case 0: mRA.choiceOfMethod = "Sub-ROI by division"; break;
	    		case 1: mRA.choiceOfMethod = "Sub-ROI by shrinkage"; break;
	    		case 2: mRA.choiceOfMethod = "Sub-ROI by division and shrinkage"; break;
	    	}
	    }

	    gd.removeAll();

		///////////////////////////////////////
		// ROI coordinates
		///////////////////////////////////////

	    if (mRA.choiceOfRoi.equalsIgnoreCase("Cuboid")) {

	    	GenericDialog cubeGD = new GenericDialog("Please input cube coordinates that result in even edge lengths!");

	    	cubeGD.addNumericField("Left X-coordinate of cuboid ", 50, 0, 6, "");
	    	cubeGD.addNumericField("Right X-coordinate of cuboid (0 corresponds to the maximum x)", 650, 0, 6, "");
	    	cubeGD.addNumericField("Topmost Y-coordinate of cuboid ", 50, 0, 6, "");
	    	cubeGD.addNumericField("Bottommost Y-coordinate of cuboid (0 corresponds to the maximum y)", 650, 0, 6, "");
	    	cubeGD.addNumericField("Z-coordinate of cuboid top ", 50, 0, 6, "");
	    	cubeGD.addNumericField("Z-coordinate of cuboid bottom ", 650, 0, 6, "");

			//show dialog
	    	cubeGD.showDialog();
		    if (cubeGD.wasCanceled()) return null;
		    else {
		    	//domain definitions
		    	mRA.cubeX1 = (int)Math.round(cubeGD.getNextNumber());
		    	mRA.cubeX2 = (int)Math.round(cubeGD.getNextNumber());
		    	mRA.cubeY1 = (int)Math.round(cubeGD.getNextNumber());
		    	mRA.cubeY2 = (int)Math.round(cubeGD.getNextNumber());
		    	mRA.cubeZ1 = (int)Math.round(cubeGD.getNextNumber());
		    	mRA.cubeZ2 = (int)Math.round(cubeGD.getNextNumber());

		    	mRA.edgeX = mRA.cubeX2 - mRA.cubeX1;
		    	mRA.edgeY = mRA.cubeY2 - mRA.cubeY1;
		    	mRA.edgeZ = mRA.cubeZ2 - mRA.cubeZ1;
		    }
		    cubeGD.removeAll();
		}

	    if (mRA.choiceOfRoi.equalsIgnoreCase("Cylinder")) {
			gd.addNumericField("X-coordinate of cylinder center ", 0, 0, 6, "");
			gd.addNumericField("Y-coordinate of cylinder center ", 0, 0, 6, "");
			gd.addNumericField("Z-coordinate of cylinder top ", 0, 0, 6, "");
			gd.addNumericField("Z-coordinate of cylinder bottom ", 0, 0, 6, "");
			gd.addNumericField("Radius of cylinder ", 0, 0, 6, "");

			//show dialog
			gd.showDialog();
		    if (gd.wasCanceled()) return null;
		    else {
		    	//domain definitions
		    	mRA.cylX = (int)Math.round(gd.getNextNumber());
		    	mRA.cylY = (int)Math.round(gd.getNextNumber());
		    	mRA.cylZ1 = (int)Math.round(gd.getNextNumber());
		    	mRA.cylZ2 = (int)Math.round(gd.getNextNumber());
		    	mRA.cylRadius = (int)Math.round(gd.getNextNumber());
		    }

		    //check for evenness
	    	mRA.edgeX = 2 * mRA.cylRadius;
	    	mRA.edgeY = 2 * mRA.cylRadius;
	    	mRA.edgeZ = mRA.cubeZ2 - mRA.cubeZ1;
	    	if (Math.floorMod(mRA.edgeZ, 2) > 0) {
	    		mRA.edgeZ--;
	    		mRA.cubeZ2--;
	    	}
		}

		///////////////////////////////////////
		// the number of steps to be investigated for the ROI
		///////////////////////////////////////

		GenericDialog stepDG = new GenericDialog("Please tell me what you want to get analysed!");

		if (mRA.choiceOfMethod.equalsIgnoreCase("Sub-ROI by division")) {
			stepDG.addNumericField("Please enter the number of division steps you want to investigate ", 2, 0);

			//show dialog
			stepDG.showDialog();
		    if (stepDG.wasCanceled()) return null;
		    else {
		    	mRA.divNumber = (int)Math.round(stepDG.getNextNumber());
		    }

		    //check if ROI is divisible by 2^stepnumber
		    String myMesg = "";

		    int maxdivisor = (int)Math.round(Math.pow(8, mRA.divNumber));
		    int restX = mRA.edgeX % maxdivisor; //Math.floorMod(mRA.edgeX,maxdivisor);
		    int restY = mRA.edgeY % maxdivisor; //Math.floorMod(mRA.edgeY,maxdivisor);
		    int restZ = mRA.edgeY % maxdivisor; //Math.floorMod(mRA.edgeZ,maxdivisor);

		    if (restX > 0) {
		    	myMesg += "X-edge is not divisible by " + maxdivisor + ". Changing rightmost X to " + (mRA.cubeX2 - restX) + ".\n\n";
		    	mRA.cubeX2 -= restX; mRA.edgeX -= restX;
		    }

		    if (restY > 0) {
		    	myMesg += "Y-edge is not divisible by " + maxdivisor + ". Changing bottommost Y to " + (mRA.cubeY2 - restY) + ".\n\n";
		    	mRA.cubeY2 -= restY; mRA.edgeY -= restY;
		    }

		    if (restZ > 0) {
		    	myMesg += "Z-edge is not divisible by " + maxdivisor + ". Changing bottommost Z to " + (mRA.cubeZ2 - restZ) + ".\n\n";
		    	mRA.cubeZ2 -= restZ; mRA.edgeZ -= restZ;
		    }

		    //show message recapitulating the choices
		    myMesg += "You are about to run a ROI analysis using " + mRA.divNumber + " ROI divisons.\n\n" +
		    		"This corresponds to 1 ROI with edge lengths of " + mRA.edgeX + " vx, " +
		    		mRA.edgeY  + " vx and " + mRA.edgeZ + " vx in X, Y, and Z directions.\n\n";

		    int roiNumber = 1;

		    for (int i = 1 ; i <= mRA.divNumber ; i++) {
		    	myMesg += "This corresponds to "+ String.format("%3.0f",Math.pow(8, i)) + " ROIs with edge lengths of " + mRA.edgeX/(i*2) +
		    			" vx, " + mRA.edgeY/(i*2)  + " vx and " + mRA.edgeZ/(i*2) +
		    			" vx in X, Y, and Z directions.\n\n";
		    	roiNumber += Math.pow(8, i);
		    }

		    myMesg += "You will analyze " + roiNumber + " different ROIs.";

		    if (IJ.showMessageWithCancel("Are you sure?", myMesg));
		    else return null;

		}

		if (mRA.choiceOfMethod.equalsIgnoreCase("Sub-ROI by shrinkage")) {
			stepDG.addNumericField("Please enter the number shrinkage steps you want to investigate ", 10, 0);
			stepDG.addNumericField("Please enter the number size of the shrinkage steps ", 10, 0);

			//show dialog
			stepDG.showDialog();
		    if (stepDG.wasCanceled()) return null;
		    else {
		    	mRA.stepNumber = (int)Math.round(stepDG.getNextNumber());
		    	mRA.stepLength = (int)Math.round(stepDG.getNextNumber());
		    }

		    //check whether step length is even
		    String myMesg = "";
		    boolean hasBeenDecreased = false;
		    if (mRA.stepLength % 2 > 0) {
		    	mRA.stepLength++;
		    	hasBeenDecreased = true;
		    }

		    //check whether starting volume is enough to conduct all steps..
		    double[] edges = {mRA.edgeX, mRA.edgeY, mRA.edgeZ};
		    if (mRA.stepNumber * mRA.stepLength >= StatUtils.min(edges)) {
		    	mRA.stepLength = (int)Math.round(Math.floor(StatUtils.min(edges) / mRA.stepNumber)) - 1;
		    	if (mRA.stepLength % 2 > 0) mRA.stepLength--;		//make sure that step-length is even.
		    	myMesg += "Step-number times step-length exceeds minimum edge length!\n";
		    	myMesg += "Step-length has therefore been reduced to " + mRA.stepLength + "!\n";
		    	myMesg += " \n";
		    }
		    else if (hasBeenDecreased) {
		    	myMesg += "Step-length must be even!\n";
		    	myMesg += "Step-length has therefore been increased to " + mRA.stepLength + "!\n";
		    	myMesg += " \n";
		    }

		    //show message recapitulating the choices
		    myMesg += "You are about to run a ROI analysis using " + mRA.stepNumber + " ROI shrinkages\n" +
		    		"with an increment of " + mRA.stepLength + " vx.\n\n" +

		    		"This corresponds to " + (mRA.stepNumber +1) + " ROI with edge lengths of ...\n";

		    for (int i = 0 ; i <= mRA.stepNumber ; i++) {
		    	myMesg += (mRA.edgeX - i*mRA.stepLength) + " vx, " + (mRA.edgeY - i*mRA.stepLength)  + " vx and " + (mRA.edgeZ - i*mRA.stepLength)+ " vx in X, Y, and Z directions.\n";
		    }

		    myMesg += "\nYou will analyze " + (mRA.stepNumber + 1) + " different ROIs.";

		    if (IJ.showMessageWithCancel("Are you sure?", myMesg));
		    else return null;

		}

		if (mRA.choiceOfMethod.equalsIgnoreCase("Sub-ROI by division and shrinkage")) {
			stepDG.addNumericField("Please enter the number division steps you want to investigate ", 1, 0);
			stepDG.addNumericField("Please enter the number shrinkage steps you want to investigate ", 15, 0);
			stepDG.addNumericField("Please enter the number size of the shrinkage steps ", 20, 0);

			//show dialog
			stepDG.showDialog();
		    if (stepDG.wasCanceled()) return null;
		    else {
		    	mRA.divNumber = (int)Math.round(stepDG.getNextNumber());
		    	mRA.stepNumber = (int)Math.round(stepDG.getNextNumber());
		    	mRA.stepLength = (int)Math.round(stepDG.getNextNumber());
		    }

		    //check if ROI is divisible by 2^stepnumber
		    String myMesg = "";

		    int maxdivisor = (int)Math.round(Math.pow(8, mRA.divNumber));
		    int restX = mRA.edgeX % maxdivisor; //Math.floorMod(mRA.edgeX,maxdivisor);
		    int restY = mRA.edgeY % maxdivisor; //Math.floorMod(mRA.edgeY,maxdivisor);
		    int restZ = mRA.edgeY % maxdivisor; //Math.floorMod(mRA.edgeZ,maxdivisor);

		    if (restX > 0) {
		    	myMesg += "X-edge is not divisible by " + maxdivisor + ". Changing rightmost X to " + (mRA.cubeX2 - restX) + ".\n\n";
		    	mRA.cubeX2 -= restX; mRA.edgeX -= restX;
		    }

		    if (restX > 0) {
		    	myMesg += "Y-edge is not divisible by " + maxdivisor + ". Changing bottommost Y to " + (mRA.cubeY2 - restY) + ".\n\n";
		    	mRA.cubeY2 -= restY; mRA.edgeY -= restY;
		    }

		    if (restX > 0) {
		    	myMesg += "Z-edge is not divisible by " + maxdivisor + ". Changing bottommost Z to " + (mRA.cubeZ2 - restZ) + ".\n\n";
		    	mRA.cubeZ2 -= restZ; mRA.edgeZ -= restZ;
		    }

		    //check whether step length is even
		    boolean hasBeenDecreased = false;
		    if (mRA.stepLength % 2 > 0) {
		    	mRA.stepLength++;
		    	hasBeenDecreased = true;
		    }

		    //check whether starting volume is enough to conduct all steps..
		    double[] edges = {mRA.edgeX, mRA.edgeY, mRA.edgeZ};
		    if (mRA.stepNumber * mRA.stepLength >= StatUtils.min(edges) / Math.pow(2, mRA.divNumber) - 1) {
		    	double smallestEdge = StatUtils.min(edges) / Math.pow(2, mRA.divNumber);
		    	mRA.stepLength = (int)Math.floor(smallestEdge / mRA.stepNumber) - 1;
		    	if (mRA.stepLength % 2 > 0) mRA.stepLength--;		//make sure that step-length is even.
		    	myMesg += "Step-number times step-length exceeds minimum edge length!\n";
		    	myMesg += "Step-length has therefore been reduced to " + mRA.stepLength + "!\n";
		    	myMesg += " \n";
		    }
		    else if (hasBeenDecreased) {
		    	myMesg += "Step-length must be even!\n";
		    	myMesg += "Step-length has therefore been increased to " + mRA.stepLength + "!\n";
		    	myMesg += " \n";
		    }

		    //check if step length is still 2 or larger
		    if (mRA.stepLength < 2) {
		    	IJ.error("This is silly! Please choose more reasonable parameters!\nBailing now...!");
		    	return null;
		    }

		    //recapitulate
		    myMesg += "You are about to run a ROI analysis using " + mRA.divNumber + " ROI divisons followed by " + mRA.stepNumber + " ROI shrinkages\n\n";
		    myMesg += "with an increment of " + mRA.stepLength + " vx.\n\n";
		    myMesg += 		"This corresponds to 1 ROI with edge lengths of " + mRA.edgeX + " vx, " +
		    		mRA.edgeY  + " vx and " + mRA.edgeZ + " vx in X, Y, and Z directions.\n\n";

		    int roiNumber = 1;

		    for (int i = 1 ; i <= mRA.divNumber ; i++) {
		    	myMesg += "This corresponds to "+ String.format("%3.0f",Math.pow(8, i)) + " ROIs with edge lengths of " + mRA.edgeX/(i*2) +
		    			" vx, " + mRA.edgeY/(i*2)  + " vx and " + mRA.edgeZ/(i*2) +
		    			" vx in X, Y, and Z directions.\n\n";
		    	roiNumber += Math.pow(8, i);
		    }

		    for (int i = 1 ; i <= mRA.stepNumber ; i++) {
		    	myMesg += (int)Math.round(Math.pow(8, mRA.divNumber)) + " ROIs with "+ (mRA.edgeX/(mRA.divNumber*2) - i*mRA.stepLength) +
		    			" vx, " + (mRA.edgeY/(mRA.divNumber*2) - i*mRA.stepLength)  + " vx and " +
		    			(mRA.edgeZ/(mRA.divNumber*2) - i*mRA.stepLength)+ " vx in X, Y, and Z directions.\n";
		    }

		    myMesg += "\nYou will analyze " + (roiNumber + mRA.stepNumber*(int)Math.round(Math.pow(8, mRA.divNumber))) + " different ROIs.";

		    if (IJ.showMessageWithCancel("Are you sure?", myMesg));
		    else return null;

		}

		stepDG.removeAll();


		///////////////////////////////////////
		// the morphological measures needed
		///////////////////////////////////////

		GenericDialog gd2 = new GenericDialog("Please tell me what you want to get analysed!");

		String[] whatIntegralMeasuresShallIInvestigate = new String[]{"Macropore Volume","Average Thickness","Fractal Dimension","Anisotropy"};
		boolean[] myIntChoices = new boolean[]{true,true,true,true};
		gd2.setInsets(20, 200, 0);gd2.addMessage("\nWhich global morphological measures should I calculate?");
		gd2.setInsets(0, 200, 0);gd2.addCheckboxGroup(3, 3, whatIntegralMeasuresShallIInvestigate, myIntChoices);

		String[] whatShallIInvestigate = new String[]{"Volume", "Euler Number","Thickness",
				"Anisotropy","Percolation","Crit. Pore Diameter"};
		boolean[] myChoices = new boolean[]{true,true,true,true,true,true};
		gd2.setInsets(20, 200, 0);gd2.addMessage("\nWhich morphological measures should I calculate for each cluster?");
		gd2.setInsets(0, 200, 0);gd2.addCheckboxGroup(3, 3, whatShallIInvestigate, myChoices);

		String[] whichImagesShallIPlot = new String[]{"Cluster Label","Volume","Thickness","Percolating Clusters","Anisotropy"};
		boolean[] myPlotChoices = new boolean[]{true,false,true,false,false};
		gd2.setInsets(20, 200, 0);gd2.addMessage("\nWhich images shall I save?");
		gd2.setInsets(0, 200, 0);gd2.addCheckboxGroup(3, 3, whichImagesShallIPlot, myPlotChoices);

		//show dialog
		gd2.showDialog();
	    if (gd2.wasCanceled()) return null;
	    else {

	    	//global measures
	    	mRA.globVolume = gd2.getNextBoolean();
			mRA.globThickness = gd2.getNextBoolean();
			mRA.calcFractal = gd2.getNextBoolean();
			mRA.globAnisotropy = gd2.getNextBoolean();

			//local measures
	    	mRA.calcVolume = gd2.getNextBoolean();
	    	mRA.calcEuler = gd2.getNextBoolean();
			mRA.calcThickness = gd2.getNextBoolean();
			mRA.calcAnisotropy = gd2.getNextBoolean();
			mRA.calcPercolation = gd2.getNextBoolean();
			mRA.calcCriticalPoreDiameter = gd2.getNextBoolean();

			//plotting options
			mRA.plotLabels = gd2.getNextBoolean();
			mRA.plotVolume = gd2.getNextBoolean();
			mRA.plotThickness = gd2.getNextBoolean();
			mRA.plotPercolation = gd2.getNextBoolean();

			mRA.performParticleAnalyses = false;
			if (mRA.calcVolume == true ||
					mRA.calcEuler == true ||
					mRA.calcThickness == true ||
					mRA.calcAnisotropy == true ||
					mRA.calcPercolation == true ||
					mRA.calcCriticalPoreDiameter == true ||

					//plotting options
					mRA.plotLabels == true ||
					mRA.plotVolume == true ||
					mRA.plotThickness == true ||
					mRA.plotPercolation == true) {

				mRA.performParticleAnalyses = true;

			}

	    }

	    //IJ.error("Plot thickness is " + mRA.plotThickness);

	    return mRA;

	}

	public class AirFilledPoresFinderMenu {

		public double resolutionInMicroMeter;
		public double matrixPotentialAtColumnBottomInCM;

	}

	public AirFilledPoresFinderMenu showAPFinderMenu() {

		//construct objects
		GenericDialog gd = new GenericDialog("PVC Column Finder Dialog");

		AirFilledPoresFinderMenu aPF = new AirFilledPoresFinderMenu();


		gd.addNumericField("Please enter the image resolution in micrometer", 0.65, 3, 6, "");
		gd.addNumericField("Please enter the tension at the bottom boundary in centimeter", 5, 3, 6, "");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	aPF.resolutionInMicroMeter = gd.getNextNumber();
	    	aPF.matrixPotentialAtColumnBottomInCM = gd.getNextNumber();

	    }

		return aPF;

	}

	public class ColumnFinderMenuReturn {

		public boolean isAlreadyNormalized;
		public boolean isSteel;
		public boolean isPVC;
		public boolean isAlu;
		public boolean try2FindColumnTopAndBottom;
		public boolean hasBevel;

		//debug settings
		public boolean debug;
		public boolean showRadialProfiles;
		public boolean showFit;

		//???
		public int segmentLength;

		//a priori fixed column properties..
		public int wallThickness;
		public int topOfColumn;
		public int bottomOfColumn;

		//fitting parameters needed for finding outer wall
		public double airWallContrast;

		//fitting parameter for finding inner wall
		public double wallSoilStdContrastThreshold;

		//parameters needed to decide whether the column is at this depth
		public double CVWallBrightnessThresh;
		public double maxCVOfWallThickness;
		public double stdThreshold;
		public double ratioBetweenInnerAndOuterRadius;
		public double percentToleratedDifference;

		//if column has already been normalized
		public double absoluteWallgrayValue;

		//optimaization parameters
		public int medianFilter;
		public double r2Thresh;
		public int maxFittingAttempts;		//attempts allowed to fiddle with "airWallContrast" and "wallSoilStdContrastThreshold" to find the edge
		public int maxNumberOfOutliers4PerimeterFits;

	}

	public ColumnFinderMenuReturn showColumnStraightenerMenu() {

		//construct objects
		GenericDialog gd = new GenericDialog("PVC Column Finder Dialog");

		ColumnFinderMenuReturn mCFS = new ColumnFinderMenuReturn();

		//for now (Feb 2017), do not show anything..
		
		//add choices to menu
	/*	String[] columnMaterials = {"steel", "PVC", "aluminium"};
		gd.addRadioButtonGroup("What kind of material is your column made of?", columnMaterials, 3, 1, columnMaterials[2]);

		gd.addNumericField("Minimum contrast between air and wall relative to image contrast (0..1)", 0.25, 3, 6, "");
		gd.addNumericField("Footprint of median filter for radial wall position search", 5, 0, 3, "");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	String myChoice = gd.getNextRadioButton();
	    	int myNumericChoice = 0;
	    	for (int i = 0 ; i < columnMaterials.length ; i++) if (myChoice.equalsIgnoreCase(columnMaterials[i])) myNumericChoice = i;
	    	switch (myNumericChoice) {
	    		case 0:	{
	    			mCFS.isSteel = true;
	    			mCFS.isPVC = false;
	    			mCFS.isAlu = false;
	    			break;
	    		}
	    		case 1: {
	    			mCFS.isSteel = false;
	    			mCFS.isPVC = true;
	    			mCFS.isAlu = false;
	    			break;
	    		}
	    		case 2: {
	    			mCFS.isSteel = false;
	    			mCFS.isPVC = false;
	    			mCFS.isAlu = true;
	    			break;
	    		}
	    	}

	    	mCFS.airWallContrast = (float)gd.getNextNumber();
	    	mCFS.medianFilter = (int)Math.round(gd.getNextNumber());
	    }*/

		//for now, assign the properties manually 
		mCFS.isSteel = false;
		mCFS.isPVC = false;
		mCFS.isAlu = true;
		
		mCFS.airWallContrast = 0.25f;
    	mCFS.medianFilter = 5;
		
	    return mCFS;

	}

	public ColumnFinderMenuReturn showColumnFinderDialog() {

		//construct objects
		GenericDialog gd = new GenericDialog("PVC Column Finder Dialog");

		ColumnFinderMenuReturn mCFS = new ColumnFinderMenuReturn();

		//add choices to menu
		gd.addCheckbox("This column has already standardized gray-values", false);
		gd.addCheckbox("These columns are bevelled!", true);
		gd.addCheckbox("Do you want the program to find top and bottom of your column?", true);

		gd.addMessage("");

		String[] columnMaterials = {"aluminium or PVC column and a BAD contrast between soil matrix and wall",
				"aluminium or PVC column and a GOOD contrast between soil matrix and wall"};
		gd.addRadioButtonGroup("What kind of material is your column made of?", columnMaterials, 2, 1, columnMaterials[1]);

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	mCFS.isAlreadyNormalized = gd.getNextBoolean();
			mCFS.hasBevel = gd.getNextBoolean();
	    	mCFS.try2FindColumnTopAndBottom = gd.getNextBoolean();

	    	String myChoice = gd.getNextRadioButton();
	    	int myNumericChoice = 0;
	    	for (int i = 0 ; i < columnMaterials.length ; i++) if (myChoice.equalsIgnoreCase(columnMaterials[i])) myNumericChoice = i;
	    	switch (myNumericChoice) {
	    		case -1:	{  //until the steel option is re-introduced
	    			mCFS.isSteel = true;
	    			mCFS.isPVC = false;
	    			mCFS.isAlu = false;
	    			break;
	    		}
	    		case 0: {
	    			mCFS.isSteel = false;
	    			mCFS.isPVC = true;
	    			mCFS.isAlu = false;
	    			break;
	    		}
	    		case 1: {
	    			mCFS.isSteel = false;
	    			mCFS.isPVC = false;
	    			mCFS.isAlu = true;
	    			break;
	    		}
	    	}

	    }

	    gd.removeAll();

		/*gd.addMessage("");
		gd.addMessage("Advanced Options");
		gd.addMessage("");*/
	    
	    GenericDialog gd2 = new GenericDialog("PVC Column Finder Dialog");

		if (mCFS.isPVC == true) gd2.addNumericField("How thick is the column wall in pixels? ", 0, 0, 6, "");
		if (!mCFS.try2FindColumnTopAndBottom) {
			gd2.addNumericField("In which layer is the top of the column? (if you do not want to search for it)", 1, 0, 4, "");
			gd2.addNumericField("In which layer is the bottom of the column? (enter '1' for the bottommost layer)", 1, 0, 4, "");
		}

		//gd2.addMessage("Criteria for finding the outer wall perimeter");
		if (mCFS.isAlreadyNormalized) gd2.addNumericField("Gray-value of column wall ", 20000, 0, 6, "");
	/*	else {
			gd2.addNumericField("Minimum contrast between air and wall relative to image contrast (0..1)", 0.3, 3, 6, "");
			gd2.addNumericField("Absolute minimum contrast between air and wall", 100, 0, 6, "");
		}*/

		/*if (mCFS.isAlu == true) {
			gd2.addMessage("Criteria for finding the inner wall perimeter");
			gd2.addNumericField("StDev. contrast threshold between wall and soil (50 - 2000) ", 500, 0, 6, "");
		}*/

		//criteria to check whether there is a wall
		gd2.addMessage("Criteria to check whether there is a wall");
		gd2.addNumericField("Ratio between inner and outer radius must be larger than ...", 0.75, 3, 6, "");
		gd2.addNumericField("Maximum allowed coefficient of variation of wall-thickness", 0.14, 4, 6, "");
		gd2.addNumericField("Maximal CV of wall gray values", 0.14, 4, 6, "");
		gd2.addNumericField("Neighboring column outline coordinates may be different by max", 0.5, 4, 6, "%");

	/*	gd2.addMessage("Filtering options");
		gd2.addNumericField("Footprint of median filter for radial wall position search", 5, 0, 3, "");
*/
		gd2.addMessage("Criteria for goodnes of fit");
		gd2.addNumericField("Maximal number of edge finding trials", 10, 0, 6, "");
		gd2.addNumericField("Maximal number of outliers during ellipse-fit", 10, 0, 6, "");
		gd2.addNumericField("Minumum R2 required for successful ellipse-fit", 0.99999, 5, 6, "");
		
		gd2.addCheckbox("Visualize column wall finds", false);

		myReference = "If you are using this plugin please cite the following references: \n\n";
		gd2.setInsets(40, 0, 0);gd2.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd2.setInsets(0, 0, 0);gd2.addMessage(myReference);
		

		//show dialog
		gd2.showDialog();

		if (gd2.wasCanceled()) return null;
		else {
			if (mCFS.isPVC) mCFS.wallThickness = (int)Math.round(gd2.getNextNumber());
			if (!mCFS.try2FindColumnTopAndBottom) {
				mCFS.topOfColumn = (int)Math.round(gd2.getNextNumber());
				mCFS.bottomOfColumn = (int)Math.round(gd2.getNextNumber());
			}

			//get parameters to get outer edge
			if (mCFS.isAlreadyNormalized) mCFS.absoluteWallgrayValue = (float)gd2.getNextNumber();
			else {
				mCFS.airWallContrast = 0.3;//(float)gd2.getNextNumber();
				mCFS.stdThreshold = 1;//(float)gd2.getNextNumber();
			}

			// get parameters to find inner edge
			if (mCFS.isAlu) mCFS.wallSoilStdContrastThreshold = 500;//gd2.getNextNumber();

			//get parameters to check whether the wall is there
			mCFS.ratioBetweenInnerAndOuterRadius = gd2.getNextNumber();
			mCFS.maxCVOfWallThickness = gd2.getNextNumber();
			mCFS.CVWallBrightnessThresh = gd2.getNextNumber();
			mCFS.percentToleratedDifference = gd2.getNextNumber();

			//get optimization parameters;
			mCFS.medianFilter = 5;	//(int)Math.round(gd2.getNextNumber());
			mCFS.maxFittingAttempts = (int)Math.round(gd2.getNextNumber());
			mCFS.maxNumberOfOutliers4PerimeterFits = (int)Math.round(gd2.getNextNumber());
			mCFS.r2Thresh = gd2.getNextNumber();

			mCFS.debug = gd2.getNextBoolean();
			
			return mCFS;
		}
	}

	public ColumnFinderMenuReturn showTuneThreshold4ColumnDetectionMenu(ApproximateColumnIllumination aCI) {

		//construct objects
		GenericDialog gd = new GenericDialog("PVC Column Finder Dialog");

		ColumnFinderMenuReturn mCFS = new ColumnFinderMenuReturn();

		gd.addMessage("The likely median gray value outside the column is " + aCI.outside);
		gd.addMessage("The likely median gray value of the column wall is " + aCI.wall);
		gd.addMessage("The likely median gray value inside the column is " + aCI.inside);
		gd.addMessage("The global contrast for this image is " + aCI.globalContrast);
		gd.addMessage("");

		double suggestedAirWallThreshold = 0.7 * ((double)(aCI.wall - aCI.outside) / (double)aCI.globalContrast);
		double suggestedWallSoilThreshold = 0.7 * ((double)(aCI.wall - aCI.inside) / (double)aCI.globalContrast);

		gd.addNumericField("Minimum contrast between air and wall relative to average gray level (0..1)",suggestedAirWallThreshold, 3, 6, "");
		gd.addNumericField("Minimum contrast between wall and soil relative to average gray level (0..1)", suggestedWallSoilThreshold, 3, 6, "");

		//show dialog
		gd.showDialog();

		if (gd.wasCanceled()) return null;
		else {

			mCFS.airWallContrast = (float)gd.getNextNumber();

			return mCFS;
		}

	}

	public RandomClusterGenerator showRandomClusterGeneratorMenu() {

		GenericDialog gd = new GenericDialog("random cluster generator menu");

		RandomClusterGenerator mRCG = new RandomClusterGenerator();

		gd.addNumericField("What is the width (x) of your domain (vx)? ", 200, 0, 4, "");
		gd.addNumericField("What is the length (y) of your domain (vx)? ", 200, 0, 4, "");
		gd.addNumericField("What is the height (z) of your domain (vx)? ", 200, 0, 4, "");

		String[] choiceOfPoroClasses = new String[3];
		choiceOfPoroClasses[0] = "Range between two bounds";
		choiceOfPoroClasses[1] = "Gaussian around pc";
		choiceOfPoroClasses[2] = "List of porosities in a file";
		gd.addChoice("Please choose porosity classes!", choiceOfPoroClasses, "Gaussian around pc");

		String[] choiceOfShape = new String[2];
		choiceOfShape[0] = "Cubic";
		choiceOfShape[1] = "Cylindric";
		gd.addChoice("Please choose a shape of your domain!", choiceOfShape, "cylinder");

		gd.addNumericField("In case you chose 'Range between two bounds', please give me a lower bound for the random porosities! ", 0.0f, 3, 5, "");
		gd.addNumericField("In case you chose 'Range between two bounds', please give me a upper bound for the random porosities! ", 0.2f, 3, 5, "");

		gd.addNumericField("In case you chose 'Gaussian around pc', please give me a relative standard deviation! ", 0.1f, 3, 5, "");

		gd.addNumericField("How many copies of this random field do you want to create? ", 1, 0, 4, "");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	mRCG.domainX = (int)Math.round(gd.getNextNumber());
	    	mRCG.domainY = (int)Math.round(gd.getNextNumber());
	    	mRCG.domainZ = (int)Math.round(gd.getNextNumber());
	    	double[] porosityRange = new double[]{gd.getNextNumber(), gd.getNextNumber()};
	    	double standardDeviation = gd.getNextNumber();
	    	int poroClassIndex = gd.getNextChoiceIndex();
	    	switch (poroClassIndex) {
				case 0: {
					mRCG.mode = "range";
					mRCG.porosityBounds = porosityRange;
					break;
					}
				case 1: {
					mRCG.mode = "Gaussian";
					mRCG.standardDeviation = standardDeviation;
					break;
				}
				case 2: {
					mRCG.mode = "predefinedList";
					break;
				}
	    	}
	    	int typeChoiceIndex = gd.getNextChoiceIndex();
		    switch (typeChoiceIndex) {
		    	case 0: {
					mRCG.shape = "Cubic";
					break;
					}
				case 1: {
					mRCG.shape = "Cylindric";
					break;
				}
	    	}
	    	mRCG.numOfCopies = (int)Math.round(gd.getNextNumber());
	    }

		return mRCG;

	}

	public SurfaceFinderReturn showSurfaceFinderMenu() {

		//construct objects
		GenericDialog gd = new GenericDialog("Soil surface finder");

		//init SurfaceFinderReturn
		SurfaceFinderReturn mSFR = new SurfaceFinderReturn();
	
		//assign pores that are open to the surface to the soil bulk volume?
		gd.addMessage("Do you want filter out small soil crumbs loosely on the soil surface?");

		gd.addNumericField("Diameter of soil crumbs to filter out (vx) ", 8, 0, 3, "");

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show the dialog and harvest the information
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	mSFR.neglectCrumbsOfLessThan = (int)Math.round(gd.getNextNumber());

	    	return mSFR;
	    }

	}

	public ThresholderMenuReturn showThresholdingDialog() {

		//construct objects
		GenericDialog gd = new GenericDialog("Thresholding dialog");

		ThresholderMenuReturn mTMR = new ThresholderMenuReturn();

		mTMR.useConstantThreshold = false; //set constant threshold option to a default false

		int filterChoices = 0;

		String[] myChoices = new String[10];
		myChoices[0] = "Otsu";
		myChoices[1] = "Renyi Entropy";
		myChoices[2] = "Maximum Entropy";
		myChoices[3] = "Mean";
		myChoices[4] = "Moments";
		myChoices[5] = "Huang";
		myChoices[6] = "IsoData";
		myChoices[7] = "DefaultIsoData";
		myChoices[8] = "IJ_IsoData";
		myChoices[9] = "UserDefinedThreshold";

		String[] mySecondaryChoices = new String[myChoices.length];
		mySecondaryChoices[0]="no two-step segmentation";
		for (int i = 1 ; i < myChoices.length ; i++) mySecondaryChoices[i] = myChoices[i - 1];

		//add threshold method choices
		gd.addMessage("Segmentation Options");

		gd.addCheckbox("Do you want to use the column outlines defined in 'inner circle' as a ROI?", false);

		//gd.addMessage("Search-string for only segmenting some of the files in this folder. Leave blank for segmenting them all.");
		//gd.addStringField("", "", 25);

		gd.addCheckbox("Do you want to cut away all gray values brighter than the wall?", true);

		gd.addMessage("Please choose a primary thresholding algorithm!");
		gd.addRadioButtonGroup("", myChoices, 5, 2, myChoices[9]);

		gd.addMessage("In case you want to perform a two-step segmentation approach, choose a secondary thresholding algorithm!");
		gd.addRadioButtonGroup("", mySecondaryChoices, 5, 2, mySecondaryChoices[9]);

		//saving options
		gd.addMessage("");
		gd.addMessage("Saving options");

		gd.addCheckbox("Save the segmented, 3-D binary images", false);

		gd.addCheckbox("Save segmentation results as overlays in some sample slices", true);

		gd.addCheckbox("Export TIFF stacks for GeoDict", false);

		String myReference = "If you are using this plugin please cite the following references: \n\n";
		gd.setInsets(40, 0, 0);gd.addMessage(myReference);
		myReference = "Koestel, J. (2017) SoilJ: An ImageJ plugin for semi-automatized image-processing of 3-D X-ray images of soil columns\n";
		myReference += "Vadose Zone Journal, ????";
		gd.setInsets(0, 0, 0);gd.addMessage(myReference);

		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {

	    	//get segmentation options
	    	mTMR.useInnerCircle = gd.getNextBoolean();

	    /*	mTMR.filterTag = gd.getNextString();
	    	if(mTMR.filterTag.isEmpty()) mTMR.filterImages = false;
	    	else mTMR.filterImages = true;*/

	    	mTMR.setMaxgray2Wallgray = gd.getNextBoolean();

	    	String choice = gd.getNextRadioButton();
	    	for (int i = 0 ; i < myChoices.length ; i++) if (myChoices[i].equalsIgnoreCase(choice)) {
	    		filterChoices = i;
	    		break;
	    	}
	    	switch (filterChoices) {
				case 0 : {mTMR.myPrimaryMethod = AutoThresholder.Method.Otsu; break;}
				case 1 : {mTMR.myPrimaryMethod = AutoThresholder.Method.RenyiEntropy; break;}
				case 2 : {mTMR.myPrimaryMethod = AutoThresholder.Method.MaxEntropy; break;}
				case 3 : {mTMR.myPrimaryMethod = AutoThresholder.Method.Mean; break;}
				case 4 : {mTMR.myPrimaryMethod = AutoThresholder.Method.Moments; break;}
				case 5 : {mTMR.myPrimaryMethod = AutoThresholder.Method.Huang; break;}
				case 6 : {mTMR.myPrimaryMethod = AutoThresholder.Method.IsoData; break;}
				case 7 : {mTMR.myPrimaryMethod = AutoThresholder.Method.Default; break;}
				case 8 : {mTMR.myPrimaryMethod = AutoThresholder.Method.IJ_IsoData; break;}
				case 9 : {mTMR.myPrimaryMethod = null;
							mTMR.useConstantThreshold = true;
							break;}
	    	}

	    	choice = gd.getNextRadioButton();
	    	for (int i = 0 ; i < myChoices.length ; i++) if (mySecondaryChoices[i].equalsIgnoreCase(choice)) {
	    		filterChoices = i;
	    		break;
	    	}
	    	switch (filterChoices) {
	    		case 0 : {mTMR.mySecondaryMethod = null; break;}
	    		case 1 : {mTMR.mySecondaryMethod = AutoThresholder.Method.Otsu; break;}
				case 2 : {mTMR.mySecondaryMethod = AutoThresholder.Method.RenyiEntropy; break;}
				case 3 : {mTMR.mySecondaryMethod = AutoThresholder.Method.MaxEntropy; break;}
				case 4 : {mTMR.mySecondaryMethod = AutoThresholder.Method.Mean; break;}
				case 5 : {mTMR.mySecondaryMethod = AutoThresholder.Method.Moments; break;}
				case 6 : {mTMR.mySecondaryMethod = AutoThresholder.Method.Huang; break;}
				case 7 : {mTMR.mySecondaryMethod = AutoThresholder.Method.IsoData; break;}
				case 8 : {mTMR.mySecondaryMethod = AutoThresholder.Method.Default; break;}
				case 9 : {mTMR.mySecondaryMethod = AutoThresholder.Method.IJ_IsoData; break;}
	    	}

	    	//get saving options
	    	mTMR.save3DImage = gd.getNextBoolean();
	    	mTMR.save4Evaluation = gd.getNextBoolean();
	    	mTMR.save4GeoDict = gd.getNextBoolean();

	      	return mTMR;
	    }
	}

}



