package soilj;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import soilj.tools.ImageManipulator;
import soilj.tools.InputOutput;
import soilj.tools.MenuWaiter;

import java.io.File;

public class SegmentThis_ extends ImagePlus implements PlugIn  {

	// loads all 3D TIFFs in the specified folder and puts them straight
	// and aligns them as such that they stand 'perfectly' straight

	public void run(String arg) {

		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ImageManipulator jIM = new ImageManipulator();
		MenuWaiter menu = new MenuWaiter();
						
		// init variables
		MenuWaiter.ThresholderMenuReturn mTMR;
		int i;
		
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();		
		ImagePlus[] outTiff = {null, null, null};
		
		//read base folder and number of 3D images
		String myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
		String[] myTiffs0 = jIO.listTiffsInFolder(new File(myBaseFolder));
		
		//ask for threshold choice
		mTMR = menu.showThresholdingDialog();
	
		//if not all tiffs shall be thresholded
		String[] myTiffs;
		if (mTMR.filterImages) {myTiffs = jIO.filterFolderList(myTiffs0, mTMR.filterTag);} else myTiffs = myTiffs0;
		
		//read gauges in the gauge folder
		String myGaugeFolder = jIO.findTheInnerCircle(myBaseFolder);		
		String[] myGauges = jIO.listGaugeFiles(myGaugeFolder, mTMR.filterTag);
		
		String myOutFolder;
		String myOutPath = null;		
		if (mTMR.useConstantThreshold == true) {
			
			mTMR = menu.showConstantThresholdingDialog(mTMR);
			String thresholds = "";
			if (mTMR.airThreshold > 0) thresholds = "" + mTMR.airThreshold;
			if (mTMR.waterThreshold > 0 & thresholds.equalsIgnoreCase("")) thresholds = "" + mTMR.waterThreshold;
			if (mTMR.waterThreshold > 0 & !thresholds.equalsIgnoreCase("")) thresholds += "_" + mTMR.waterThreshold;
			if (mTMR.stoneThreshold > 0 & thresholds.equalsIgnoreCase("")) thresholds = "" + mTMR.stoneThreshold;
			if (mTMR.stoneThreshold > 0 & !thresholds.equalsIgnoreCase("")) thresholds += "_" + mTMR.stoneThreshold;
			
			myOutFolder = "ConstantThreshold" + thresholds;
			
		}
		else {myOutFolder = mTMR.myPrimaryMethod.toString();		
			if (mTMR.filterImages) myOutFolder = mTMR.filterTag + "_" + myOutFolder;
			if (mTMR.mySecondaryMethod != null) myOutFolder = mTMR.myPrimaryMethod.toString() + "_" + mTMR.mySecondaryMethod.toString();
			if (mTMR.setMaxGrey2WallGrey) myOutFolder = myOutFolder + "_WallIsMaxGrey";
		}		
		myOutPath = myBaseFolder + myOutFolder;						
		new File(myOutPath).mkdir();	
		
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();
			
			//load image
			nowTiff = jIO.openTiff3D(myBaseFolder + "\\" + myTiffs[i]);	
			
			//select the correct gauge and surface files			
			int[] myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(myTiffs[i], myGauges, null);		
			
			//apply segmentation			
			outTiff = jIM.applyAChosenThresholdingMethod(nowTiff, myGauges[myGandS[0]], mTMR, myOutPath);			 
		
			//save result
			if (mTMR.save3DImage == true) jIO.tiffSaver(myOutPath, myTiffs[i], outTiff[0]);
			if (mTMR.saveWaterPhase == true & mTMR.waterThreshold > 0) jIO.tiffSaver(myOutPath, myTiffs[i] + "_Water", outTiff[1]);
			if (mTMR.saveStonePhase == true & mTMR.stoneThreshold > 0) jIO.tiffSaver(myOutPath, myTiffs[i] + "_Stone", outTiff[2]);
			if (mTMR.save4GeoDict == true) jIO.saveAsTiffStack4GeoDict(myOutPath, myTiffs[i], outTiff[0]);
			
			if (i == 0) jIO.writeStringIntoAsciiFile(myOutPath + "\\Path2Gauge.txt", myGaugeFolder + "\n" + mTMR.filterTag); 	//also save the gauge path..
 			
		}
		
	}
	
}