package SoilJ.tools;

import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.AutoThresholder;
import ij.process.AutoThresholder.Method;
import SoilJ.tools.ObjectDetector.ApproximateColumnIllumination;

public class MenuWaiter implements PlugIn {

	
	public class BeamDeHardeningReturn {
	
		public Boolean isSteelColumn = true; 
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

	public class PoreSpaceAnalyserOptions {
		
		public String choiceOfRoi;
		
		public int heightOfRoi;
		public int edgeLengthOfRoi;
		public int cutAwayFromTop;
		public int cutAwayFromWall;
		
		public boolean hasSurfaceFiles;
		
		public boolean globVolume;
		public boolean	globSurface;
		public boolean globThickness;
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
		
		public boolean bareSoil;
		public boolean croppedSoil;
		
		public boolean binary;
		public boolean greyscale;
		
		public int airThresh;
		public int threshold;
		
		int neglectCrumbsOfLessThan;
		
		public double radius4PoreFilter;
		public double radiusIncreaser;
		
	}

	public class ThresholderMenuReturn {
		
		public boolean filterImages = false;
		public String filterTag;
		public boolean setMaxGrey2WallGrey;
		
		public Method myPrimaryMethod = AutoThresholder.Method.Default;
		public Method mySecondaryMethod = AutoThresholder.Method.Default;
		
		public int airThreshold = 0;
		public int waterThreshold = 0;
		public int stoneThreshold = 0;
		
		public boolean useConstantThreshold;
		
		public boolean save3DImage;
		public boolean save4Evaluation;
		public boolean save4GeoDict;	
		
		public boolean saveWaterPhase;
		public boolean saveStonePhase;
		
	}

	public class OMFinderSettings {
		
		int minGreyValue;
		int maxGreyValue;
		int overlap;
		int windowSize;
		
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
		
		gd.addCheckbox("Do you have steel columns?", false);
		
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
		
		gd.addNumericField("Give a number that the grey value sampled for the air phase has to be represneted with at least!", mBDH.airPhaseClassMemberMinimum, 0, 6, "");
		
		String[] choiceOfMasking = new String[3];
		choiceOfMasking[0] = "None";
		choiceOfMasking[1] = "Otsu";
		choiceOfMasking[2] = "MaxEntropy";
		gd.addChoice("Please choose a masking method for the non-matrix!", choiceOfMasking, "Otsu");	
		
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else { 
	    	mBDH.isSteelColumn = gd.getNextBoolean();	    
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
	
	public ClipperMenuReturn showClipperDialog() {
		
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
		
	}

	public ThresholderMenuReturn showConstantThresholdingDialog(ThresholderMenuReturn mTMR) {
		
		//construct objects
		GenericDialog gd = new GenericDialog("Thank you for choosing a user defined (constant) threshold for all images!");
				
		gd.addMessage("Make sure that you have run 'NormalizeTheseImages' when you apply this option!!!");
		
		gd.addNumericField("Enter the constant threshold for segmenting the air-phase ", 0, 0, 5, "");
		
		gd.addNumericField("Enter the constant threshold for segmenting the water-phase and fresh organic matter (optional) ", 0, 0, 5, "");
		
		gd.addNumericField("Enter the constant threshold for segmenting the mineral-phase, i.e. sand and gravel (optional) ", 0, 0, 5, "");
		
		gd.addMessage("Which 3-D binary images do you want to save?");
		
		gd.addCheckbox("Air-phase", false);
		
		gd.addCheckbox("Water-phase", false);
		
		gd.addCheckbox("Stone-phase", false);
		
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else { 
	    	
	    	mTMR.airThreshold = (int)Math.round(gd.getNextNumber());
	    	
	    	mTMR.waterThreshold = (int)Math.round(gd.getNextNumber());
	    	
	    	mTMR.stoneThreshold = (int)Math.round(gd.getNextNumber());
	    	
	    	mTMR.save3DImage = gd.getNextBoolean();
	    	
	    	mTMR.saveWaterPhase = gd.getNextBoolean();
	    	
	    	mTMR.saveStonePhase = gd.getNextBoolean();
	    	
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

	public class HistogramMenuReturn {
		
		public Boolean filterImages = false;
		public String filterTag;
		
		public Boolean saveAsciis = true;
		public Boolean saveImages = false;
		public Boolean calcLumpedHistogram = false;
	
	}
	
	public HistogramMenuReturn showHistogramDialog() {		
	
		//construct objects
		GenericDialog gd = new GenericDialog("Thresholding dialog");
		
		HistogramMenuReturn hMR = new HistogramMenuReturn();
		
		//add choices to menu
		gd.addMessage("Search-string for only segmenting some of the files in this folder. Leave blank for segmenting them all.");	
		gd.addStringField("", "", 25);
	
		gd.addCheckbox("Save histograms as ASCII files?", true);
		
		gd.addCheckbox("Save histograms as images?", false);
		
		gd.addCheckbox("Calculate lumped histogram of all images in this folder?", false);
		gd.addMessage("(This option only makes sense if your images have normalized grey values, e.g. Hounsfield units)");	
		
		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else {
	    	
	    	hMR.filterTag = gd.getNextString();
	    	if (hMR.filterTag.hashCode() == 0) hMR.filterImages = false;
	    	else hMR.filterImages = true;
	    	
	    	hMR.saveAsciis = gd.getNextBoolean();
	    	hMR.saveImages = gd.getNextBoolean();
	    	hMR.calcLumpedHistogram = gd.getNextBoolean();
	    	
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
	
	public class NormalizerReferences {
		
		public String material;
		
		public double lowerReference;
		public double upperReference;
		public double lowerTarget;
		public double upperTarget;
		
		public double[] originalLower;
		public double[] originalUpper;		
		
		public String lowerTag;
		public String upperTag;
		
	}

	public NormalizerReferences showNormalizerReferencesMenu() {
		
		MenuWaiter.NormalizerReferences mNR = new NormalizerReferences();
		
		//construct objects
		GenericDialog gd = new GenericDialog("Choose parameters for statistical region merging!");
		
		// set up dialog box		
		String[] columnMaterials = {"steel", "PVC", "aluminium"};
		gd.addRadioButtonGroup("What kind of material is your column made of?", columnMaterials, 3, 1, columnMaterials[1]);
		
		gd.addNumericField("Quantile for lower normalization reference value (put 0 for grey-values of wall)", 0.001, 3, 6, "");
		gd.addNumericField("Quantile for upper normalization reference value (put 0 for grey-values of wall)", 0, 3, 6, "");
		
		gd.addMessage("");
		
		gd.addNumericField("Quantile for lower normalization target value", 5000, 0, 6, "");
		gd.addNumericField("Quantile for upper normalization target value", 20000, 0, 6, "");
		
		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else { 
	    	
	    	mNR.material = gd.getNextRadioButton();
	    	
	    	mNR.lowerReference = gd.getNextNumber();
	    	mNR.upperReference = gd.getNextNumber();
	    	mNR.lowerTarget = gd.getNextNumber();
	    	mNR.upperTarget = gd.getNextNumber();
	    	
	    	String mylo = String.format("%1.3f", mNR.lowerReference);	    	
	    	if (mNR.lowerReference == 0) mNR.lowerTag = "Wall";
	    	else mNR.lowerTag = "Quantile" + mylo.substring(2, 5);
	    	
	    	String myup = String.format("%1.3f", mNR.upperReference);
	    	if (mNR.upperReference == 0) mNR.upperTag = "Wall";
	    	else mNR.upperTag = "Quantile" + myup.substring(2, 5);
	    	
	    	if (mNR.upperReference == 0 ) {
	    		return mNR;
	    	}
	    	else {	    			
	    		if (mNR.lowerReference >= mNR.upperReference | mNR.lowerTarget >= mNR.upperTarget) {	    			
	    			return null;
	    		}
	    		else {
	    			return mNR;
	    		}
	    	}	    	      	
	    }	
	}

	public OMFinderSettings showOMFinderMenu() {
		
		OMFinderSettings oMF = new OMFinderSettings();
		
		GenericDialog gd = new GenericDialog("Fresh organic matter extractor dialog");
		
		gd.addNumericField("What is the minimum grey value of the organic matter? ", 20, 0, 6, "");
		gd.addNumericField("What is the maximum grey value of the organic matter?", 65, 0, 4, "");
		gd.addNumericField("Please tell me the window sizes you want to apply for the extraction ", 15, 0, 4, "");
		gd.addNumericField("Please choose an overlap", 5, 0, 4, "");		
		
		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else { 
	    	
	    	oMF.minGreyValue = (int)Math.round(gd.getNextNumber());
	    	oMF.maxGreyValue = (int)Math.round(gd.getNextNumber());
	    	oMF.windowSize = (int)Math.round(gd.getNextNumber()); 
	    	oMF.overlap = (int)Math.round(gd.getNextNumber());			
	    				
			return oMF;
	    }		
		
	}

	public PoreSpaceAnalyserOptions showPoreSpaceAnalyzerMenu() {
		
		//construct objects
		GenericDialog gd = new GenericDialog("Please tell me what you want to get analysed!");		
		
		PoreSpaceAnalyserOptions mPSAO = new PoreSpaceAnalyserOptions();
		
		String[] choiceOfRoi = new String[4];
		choiceOfRoi[0] = "No ROI but the entire image!";
		choiceOfRoi[1] = "Cylindrical column with outlines defined under 'Inner Circle'";
		choiceOfRoi[2] = "Cylindrical ROI";		
		choiceOfRoi[3] = "Cuboid ROI";
		gd.addRadioButtonGroup("Please choose a region of interest!", choiceOfRoi, 4, 1, "No ROI but the entire image!");		
		
		gd.addMessage("In case you have chosen 'Cylindrical column...' ");
		gd.addNumericField("How many voxels shall I cut away from the wall? ", 0, 0, 6, "");
		gd.addNumericField("How many voxels shall I cut away from the top? ", 0, 0, 6, "");
		gd.addNumericField("How tall shall the subsample be in voxels? ", 0, 0, 6, "");	//Offer sub samples were 615
		
		gd.addMessage("In case you have chosen one of the other ROI options (cylindrical or cuboid) ...");
		gd.addNumericField("What should the diameter or edge length of your region of interest be in voxels? ", 0, 0, 6, "");
		gd.addNumericField("How tall shall the subsample be in voxels? ", 0, 0, 6, "");	//Offer sub samples were 615
	
		String[] whatIntegralMeasuresShallIInvestigate = new String[]{"Macropore Volume","Macropore Surfaces","Average Thickness",
				"Chi","Fractal Dimension","Anisotropy"};		
		boolean[] myIntChoices = new boolean[]{false,false,false,false,true,false};
		gd.setInsets(20, 200, 0);gd.addMessage("\nWhich global morphological measures should I calculate?");		
		gd.setInsets(0, 200, 0);gd.addCheckboxGroup(3, 3, whatIntegralMeasuresShallIInvestigate, myIntChoices);	
		
		String[] whatShallIInvestigate = new String[]{"Volume", "Surface Area","Moments of Inertia","Unit Vectors","Euler Number","Thickness",
				"Correlation Length","Anisotropy","Skeleton","Inclination","Percolation"};		
		boolean[] myChoices = new boolean[]{false,false,false,false,false,false,false,false,false,false,false};
		gd.setInsets(20, 200, 0);gd.addMessage("\nWhich morphological measures should I calculate for each cluster?");		
		gd.setInsets(0, 200, 0);gd.addCheckboxGroup(4, 3, whatShallIInvestigate, myChoices);	
		
		String[] whichImagesShallIPlot = new String[]{"Cluster Label","Volume","Surface","Thickness","Percolating Clusters","Inclination","Correlation Length","Anisotropy"};		
		boolean[] myPlotChoices = new boolean[]{true,false,false,false,false,false,false,false,false};		
		gd.setInsets(20, 200, 0);gd.addMessage("\nWhich images shall I save?");		
		gd.setInsets(0, 200, 0);gd.addCheckboxGroup(3, 3, whichImagesShallIPlot, myPlotChoices);	
		
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
	    	mPSAO.hasSurfaceFiles = false;
	    	switch (myChoiceIndex) {
	    		case 0: mPSAO.choiceOfRoi = "Everything!"; break;
	    		case 1: mPSAO.choiceOfRoi = "RealSample"; mPSAO.hasSurfaceFiles = true; break;
	    		case 2: mPSAO.choiceOfRoi = "Cylinder"; break;
	    		case 3: mPSAO.choiceOfRoi = "Cuboid"; break;	    		
	    	}	    	
	    	
	    	//domain definitions
	    	mPSAO.cutAwayFromWall = (int)Math.round(gd.getNextNumber());
	    	mPSAO.cutAwayFromTop = (int)Math.round(gd.getNextNumber());	    	
	    	mPSAO.heightOfRoi = (int)Math.round(gd.getNextNumber());
	    	
	    	if (myChoiceIndex > 1) {
	    		mPSAO.edgeLengthOfRoi = (int)Math.round(gd.getNextNumber());	    	
	    		mPSAO.heightOfRoi = (int)Math.round(gd.getNextNumber());
	    	}
	    		    	
	    	//global measures
	    	mPSAO.globVolume = gd.getNextBoolean();
	    	mPSAO.globSurface = gd.getNextBoolean();
			mPSAO.globThickness = gd.getNextBoolean();
			mPSAO.calcChi = gd.getNextBoolean();
			mPSAO.calcFractal = gd.getNextBoolean();
			mPSAO.globAnisotropy = gd.getNextBoolean();			
			
			//local measures
	    	mPSAO.calcVolume = gd.getNextBoolean();
	    	mPSAO.calcSurface = gd.getNextBoolean();
			mPSAO.calcMoments = gd.getNextBoolean();
			mPSAO.calcUnitVectors = gd.getNextBoolean();			
			mPSAO.calcEuler = gd.getNextBoolean();
			mPSAO.calcThickness = gd.getNextBoolean();
			mPSAO.calcCorrLength = gd.getNextBoolean();
			mPSAO.calcAnisotropy = gd.getNextBoolean();
			mPSAO.calcSkeleton = gd.getNextBoolean();
			mPSAO.calcInclination = gd.getNextBoolean();
			mPSAO.calcPercolation = gd.getNextBoolean();
			
			//plotting options
			mPSAO.plotLabels = gd.getNextBoolean();
			mPSAO.plotVolume = gd.getNextBoolean();
			mPSAO.plotSurface = gd.getNextBoolean();	
			mPSAO.plotThickness = gd.getNextBoolean();
			mPSAO.plotPercolation = gd.getNextBoolean();
			mPSAO.plotInclination = gd.getNextBoolean();
			mPSAO.plotCorrelationLength = gd.getNextBoolean();
			mPSAO.plotAnisotropy = gd.getNextBoolean();
			
			mPSAO.performParticleAnalyses = false;
			if (mPSAO.calcVolume == true || 
					mPSAO.calcSurface == true || 
					mPSAO.calcMoments == true || 
					mPSAO.calcUnitVectors == true || 			
					mPSAO.calcEuler == true || 
					mPSAO.calcThickness == true || 
					mPSAO.calcCorrLength == true || 
					mPSAO.calcAnisotropy == true || 
					mPSAO.calcSkeleton == true || 
					mPSAO.calcInclination == true || 
					mPSAO.calcPercolation == true || 
			
					//plotting options
					mPSAO.plotLabels == true || 
					mPSAO.plotVolume == true || 
					mPSAO.plotSurface == true || 	
					mPSAO.plotThickness == true || 
					mPSAO.plotPercolation == true || 
					mPSAO.plotInclination == true || 
					mPSAO.plotCorrelationLength == true || 
					mPSAO.plotAnisotropy == true) {
				
				mPSAO.performParticleAnalyses = true;
				
			}
	      	
	    }			
		
	    return mPSAO;
	    
	}

	public class ColumnFinderMenuReturn {
		
		public boolean isSteel;
		public boolean isPVC;
		public boolean isAlu;
		public boolean try2FindColumnTopAndBottom;
		public boolean debug;  
		
		public int wallThickness;
		public int topOfColumn;
		public int bottomOfColumn;
		
		public float airWallContrast;
		public float minAbsoluteGradientbetweenAirAndWall;
		public float wallSoilContrast;
		public double widthOfWallRelativeToOuterRadius;
		public double minWallGreyValue;
		public double coeffVarOfWallBrightness;
		public int medianFilter;
		
		public double minR24PerimeterFit;
		public double maxCVOfWallThickness;
		public int maxNumberOfOutliers4PerimeterFits;
	
	}
	
	public ColumnFinderMenuReturn showColumnStraightenerMenu() {
		
		//construct objects
		GenericDialog gd = new GenericDialog("PVC Column Finder Dialog");
						
		ColumnFinderMenuReturn mCFS = new ColumnFinderMenuReturn();
		
		//add choices to menu
		String[] columnMaterials = {"steel", "PVC", "aluminium"};
		gd.addRadioButtonGroup("What kind of material is your column made of?", columnMaterials, 3, 1, columnMaterials[2]);
				
		gd.addNumericField("Minimum contrast between air and wall relative to image contrast (0..1)", 0.25, 3, 6, "");
		gd.addNumericField("Footprint of median filter for radial wall position search", 5, 0, 3, "");
		
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
	    }	
		
	    return mCFS;
	    
	}
	
	public ColumnFinderMenuReturn showColumnFinderDialog() {
		
		//construct objects
		GenericDialog gd = new GenericDialog("PVC Column Finder Dialog");
						
		ColumnFinderMenuReturn mCFS = new ColumnFinderMenuReturn();
		
		//add choices to menu
		String[] columnMaterials = {"steel", "PVC", "aluminium"};
		gd.addRadioButtonGroup("What kind of material is your column made of?", columnMaterials, 3, 1, columnMaterials[2]);
		
		gd.addCheckbox("Do you want the program to find top and bottom of your column?", true);
		
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
	    	
	    	mCFS.try2FindColumnTopAndBottom = gd.getNextBoolean();
		
	    }	
	    
	    gd.removeAll();
		
		gd.addMessage("");
		gd.addMessage("Advanced Options");
		gd.addMessage("");
		
		if (mCFS.isPVC == true) gd.addNumericField("How thick is the PVC column wall in pixels? ", 0, 0, 6, "");
		if (!mCFS.try2FindColumnTopAndBottom) {
			gd.addNumericField("In which layer is the top of the column? (if you do not want to search for it)", 1, 0, 4, "");
			gd.addNumericField("In which layer is the bottom of the column? (enter '1' for the bottommost layer)", 1, 0, 4, "");
		}
		
		gd.addMessage("");		
		gd.addNumericField("Minimum contrast between air and wall relative to image contrast (0..1)", 0.04, 3, 6, "");
		gd.addNumericField("Absolute minimum contrast between air and wall", 1000, 0, 6, "");		
		gd.addNumericField("Minimum contrast between wall and soil relative to image contrast (0..1)", 0.015, 3, 6, "");
		gd.addNumericField("Ratio between inner and outer radius (0..1)", 0.87, 3, 6, "");
		gd.addNumericField("Minimum grey value of column wall relative to image contrast  (0..1)", 0.35, 3, 6, "");
		gd.addNumericField("Approximate coefficient of variations of wall grey values", 0.035, 3, 6, "");
		gd.addNumericField("Footprint of median filter for radial wall position search", 5, 0, 3, "");
		gd.addMessage("");
		gd.addNumericField("Minumum R2 required for successful ellipse-fit", 0.9999, 5, 6, "");		
		gd.addNumericField("Maximal number of outliers during ellipse-fit", 4, 0, 6, "");
		gd.addNumericField("Maximum allowed coefficient of variation of wall-thickness", 0.1, 4, 6, "");
		
		//show dialog		
		gd.showDialog();
		
		if (gd.wasCanceled()) return null;
		else {			
			if (mCFS.isPVC) mCFS.wallThickness = (int)Math.round(gd.getNextNumber());
			if (!mCFS.try2FindColumnTopAndBottom) {
				mCFS.topOfColumn = (int)Math.round(gd.getNextNumber());
				mCFS.bottomOfColumn = (int)Math.round(gd.getNextNumber());
			}					
						
			mCFS.airWallContrast = (float)gd.getNextNumber();
			mCFS.minAbsoluteGradientbetweenAirAndWall = (float)gd.getNextNumber(); 
			mCFS.wallSoilContrast = (float)gd.getNextNumber();
			mCFS.widthOfWallRelativeToOuterRadius = gd.getNextNumber();
			mCFS.minWallGreyValue = gd.getNextNumber();
			mCFS.coeffVarOfWallBrightness = gd.getNextNumber();
			mCFS.medianFilter = (int)Math.round(gd.getNextNumber());
			
			mCFS.minR24PerimeterFit = gd.getNextNumber();
			mCFS.maxNumberOfOutliers4PerimeterFits = (int)Math.round(gd.getNextNumber());
			mCFS.maxCVOfWallThickness = gd.getNextNumber();
			
			return mCFS;
		}
	}
		
	public ColumnFinderMenuReturn showTuneThreshold4ColumnDetectionMenu(ApproximateColumnIllumination aCI) {
		
		//construct objects
		GenericDialog gd = new GenericDialog("PVC Column Finder Dialog");
								
		ColumnFinderMenuReturn mCFS = new ColumnFinderMenuReturn();
		
		gd.addMessage("The likely median grey value outside the column is " + aCI.outside);	
		gd.addMessage("The likely median grey value of the column wall is " + aCI.wall);	
		gd.addMessage("The likely median grey value inside the column is " + aCI.inside);
		gd.addMessage("The global contrast for this image is " + aCI.globalContrast);
		gd.addMessage("");
		
		double suggestedAirWallThreshold = 0.7 * ((double)(aCI.wall - aCI.outside) / (double)aCI.globalContrast);
		double suggestedWallSoilThreshold = 0.7 * ((double)(aCI.wall - aCI.inside) / (double)aCI.globalContrast);
		
		gd.addNumericField("Minimum contrast between air and wall relative to average grey level (0..1)",suggestedAirWallThreshold, 3, 6, "");
		gd.addNumericField("Minimum contrast between wall and soil relative to average grey level (0..1)", suggestedWallSoilThreshold, 3, 6, "");
		
		//show dialog		
		gd.showDialog();
		
		if (gd.wasCanceled()) return null;
		else {			
		
			mCFS.airWallContrast = (float)gd.getNextNumber();
			mCFS.wallSoilContrast = (float)gd.getNextNumber();
			
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

	public double[] showSchlueterDialog() {		
	
		//construct objects
		GenericDialog gd = new GenericDialog("Choose parameters for statistical region merging!");
		
		double skip = 0;
		double cutoffGreyValue = 0.5;
		double lcutoffGreyValue = 0.005;
		double alpha = 1.5;
		
		double[] filterChoices = new double[4];
		
		// set up dialog box
		gd.addNumericField("# of horizontal slices to be skipped at both ends: ", skip, 0, 6, "");
		gd.addNumericField("upper cutoff grey value: ", cutoffGreyValue, 2, 6, "");
		gd.addNumericField("lower cutoff grey value: ", lcutoffGreyValue, 3, 6, "");		
		gd.addNumericField("Schlueter alpha: ", alpha, 2, 6, "");
	
		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else { 
	    	skip = gd.getNextNumber();	    	
	    	cutoffGreyValue = gd.getNextNumber();
	    	lcutoffGreyValue = gd.getNextNumber();
	    	alpha = gd.getNextNumber();
	    	
	    	filterChoices[0] = skip;
	      	filterChoices[1] = cutoffGreyValue;
	      	filterChoices[2] = lcutoffGreyValue;
	        filterChoices[3] = alpha;
	      	
	      	return filterChoices;
	    }			
	}

	public double[] showSchlueterSpecialDialog() {		
	
		//construct objects
		GenericDialog gd = new GenericDialog("Choose parameters for statistical region merging!");
		
		double skip = 0;
		double cutoffGreyValue = 0.5;
		double lcutoffGreyValue = 0.005;
		
		double[] filterChoices = new double[3];
		
		// set up dialog box
		gd.addNumericField("# of horizontal slices to be skipped at both ends: ", skip, 0, 6, "");
		gd.addNumericField("upper cutoff grey value: ", cutoffGreyValue, 2, 6, "");
		gd.addNumericField("lower cutoff grey value: ", lcutoffGreyValue, 3, 6, "");		
	
		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else { 
	    	skip = gd.getNextNumber();	    	
	    	cutoffGreyValue = gd.getNextNumber();
	    	lcutoffGreyValue = gd.getNextNumber();
	    	
	    	filterChoices[0] = skip;
	      	filterChoices[1] = cutoffGreyValue;
	      	filterChoices[2] = lcutoffGreyValue;
	      	
	      	return filterChoices;
	    }			
	}

	public double[] showStatisticalRegionMergingDialog() {		
	
		//construct objects
		GenericDialog gd = new GenericDialog("Choose parameters for statistical region merging!");
		
		double Q = 50;
		double vR = 3;
		double vH = 10; 
		
		double[] filterChoices = new double[3];
		
		// set up dialog box
		gd.addNumericField("Q: ", Q, 1, 6, "");
		gd.addNumericField("# of vertical regions: ", vR, 1, 6, "");
		gd.addNumericField("# of horizontal regions: ", vH, 1, 6, "");
	
		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else { 
	    	Q = gd.getNextNumber();
	    	vR = gd.getNextNumber();
	    	vH = gd.getNextNumber();
	    	
	    	filterChoices[0] = Q;
	      	filterChoices[1] = vR;
	      	filterChoices[2] = vH;
	      	
	      	return filterChoices;
	    }			
	}

	public double[] showSubDivisionDialog() {		
	
		//construct objects
		GenericDialog gd = new GenericDialog("Choose parameters for statistical region merging!");
		
		double skip = 25;
		double sT = 100;
		double vH = 1; 
		double cutoffGreyValue = 0.5;
		double lcutoffGreyValue = 0.005;
		double stepSize = 0;
		double alpha = 1.5;
		
		double[] filterChoices = new double[7];
		
		// set up dialog box
		gd.addNumericField("# of horizontal slices to be skipped at both ends: ", skip, 0, 6, "");
		gd.addNumericField("smoother thickness: ", sT, 0, 6, "");
		gd.addNumericField("# of horizontal regions: ", vH, 0, 6, "");
		gd.addNumericField("upper cutoff grey value: ", cutoffGreyValue, 2, 6, "");
		gd.addNumericField("lower cutoff grey value: ", lcutoffGreyValue, 3, 6, "");
		gd.addNumericField("step size: ", stepSize, 0, 6, "");
		gd.addNumericField("Schlueter alpha: ", alpha, 2, 6, "");
	
		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else { 
	    	skip = gd.getNextNumber();
	    	sT = gd.getNextNumber();
	    	vH = gd.getNextNumber();
	    	cutoffGreyValue = gd.getNextNumber();
	    	lcutoffGreyValue = gd.getNextNumber();
	    	stepSize = gd.getNextNumber();
	    	alpha = gd.getNextNumber();
	    	
	    	filterChoices[0] = skip;
	      	filterChoices[1] = sT;
	      	filterChoices[2] = vH;
	      	filterChoices[3] = cutoffGreyValue;
	      	filterChoices[4] = lcutoffGreyValue;
	      	filterChoices[5] = stepSize;
	        filterChoices[6] = alpha;
	      	
	      	return filterChoices;
	    }			
	}

	public SurfaceFinderReturn showSurfaceFinderMenu() {
		
		//construct objects
		GenericDialog gd = new GenericDialog("Please tell me some things about your soil columns..");		
		
		//init SurfaceFinderReturn
		SurfaceFinderReturn mSFR = new SurfaceFinderReturn();
		
		//dry or wet column
		String[] mySoils = new String[2];
		mySoils[0] = "I have bare soil without roots.";
		mySoils[1] = "I have cropped soil with many roots.";
		
		gd.addRadioButtonGroup("What kind of soil do you have?", mySoils, 2, 1, mySoils[0]);
		
		String[] myImgType = new String[2];
		myImgType[0] = "I will provide a binary image for surface detection.";
		myImgType[1] = "I will provide a grey-scale image for surface detection.";
		
		gd.addRadioButtonGroup("What kind of image do you have?", myImgType, 2, 1, myImgType[0]);
		
		gd.addMessage("In case you want to use a grey-scale image. IMPORTANT: this choice requires images with normalized grey scales!");
		
		gd.addNumericField("Please provide a segmentation threshold ", 0, 0, 6, "");
		
		//assign pores that are open to the surface to the soil bulk volume? 		
		gd.addMessage("Do you want filter out small soil crumbs loosely on the soil surface?");
		
		gd.addNumericField("Diameter of soil crumbs to filter out (vx) ", 8, 0, 3, "");
		
		//assign pores that are open to the surface to the soil bulk volume? 		
		gd.addMessage("Do you want to assign pores that open to the surface to the soil bulk volume?");
		
		gd.addNumericField("Diameter of pores you want to assign to the bulk soil (vx) ", 0, 0, 3, "");
		
		//show the dialog and harvest the information
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else { 	    	
	    	
	    	String stateOfSoil = gd.getNextRadioButton();	    	
	    	if (stateOfSoil == mySoils[0]) {
	    		mSFR.bareSoil = true;
	    	}
	    	else {
	    		mSFR.croppedSoil = true;
	    	}
	    	
	    	String typeOfImage = gd.getNextRadioButton();	    	
	    	if (typeOfImage == myImgType[0]) {
	    		mSFR.binary = true;	    		
	    	}
	    	else {	    		
	    		mSFR.greyscale = true;
	    	}
	    	
	    	mSFR.threshold = (int)Math.round(gd.getNextNumber());
	    	
	    	mSFR.neglectCrumbsOfLessThan = (int)Math.round(gd.getNextNumber());
	    	
	    	mSFR.radius4PoreFilter = gd.getNextNumber();
	
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
				
		gd.addMessage("Search-string for only segmenting some of the files in this folder. Leave blank for segmenting them all.");	
		gd.addStringField("", "", 25);
				
		gd.addCheckbox("Do you want to cut away all grey values brighter than the wall?", true);
			
		gd.addMessage("Please choose a primary thresholding algorithm!");
		gd.addRadioButtonGroup("", myChoices, 5, 2, myChoices[9]);
		
		gd.addMessage("In case you want to perform a two-step segmentation approach, choose a secondary thresholding algorithm!");
		gd.addRadioButtonGroup("", mySecondaryChoices, 5, 2, mySecondaryChoices[9]);
				
		//saving options
		gd.addMessage("");
		gd.addMessage("Saving options");
		
		gd.addCheckbox("Save the segmented, 3-D binary images", false);
		
		gd.addCheckbox("Save segmented cross-sections at 20%, 40%, 60% and 80% for evaluating the segmentation results", true);
		
		gd.addCheckbox("Export TIFF stacks for GeoDict", false);
			
		//show dialog
		gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    else { 
	
	    	//get segmentation options	    	
	    	mTMR.filterTag = gd.getNextString();	    	
	    	if(mTMR.filterTag.isEmpty()) mTMR.filterImages = false;
	    	else mTMR.filterImages = true;
	    	
	    	mTMR.setMaxGrey2WallGrey = gd.getNextBoolean();	  
	    	
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
