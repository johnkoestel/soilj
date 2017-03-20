package SoilJ_;

import ij.IJ;

/**
 *SoilJ is a collection of ImageJ plugins for the semi-automatized processing of 3-D X-ray images of soil columns
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

import ij.ImagePlus;
import ij.plugin.PlugIn;
import SoilJ.tools.ImageManipulator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;

import java.io.File;

/** 
 * SegmentThis is a SoilJ plugin offering several choices for image segmentation.
 * 
 * @author John Koestel
 *
 */

public class ImageSegmentation_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = ImageSegmentation_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
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
	    String[] myTiffs0 = null;String myBaseFolder = null;
	    while (myTiffs0 == null) {
	    	myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
			if (myBaseFolder == null) return;
			myTiffs0 = jIO.listTiffsInFolder(new File(myBaseFolder));
			if (myTiffs0 == null) {
				IJ.error("Please choose a folder with TIFF images or cancel. Thank you.");	
			}
		}
		
		//ask for threshold choice
		mTMR = menu.showThresholdingDialog();
		if (mTMR == null) return;
	
		//if not all tiffs shall be thresholded
		String[] myTiffs;
		if (mTMR.filterImages) {
			myTiffs = jIO.filterFolderList(myTiffs0, mTMR.filterTag);
		} 
		else myTiffs = myTiffs0;
		
		//collect all file informations
		InputOutput.MyFolderCollection mFC = jIO.new MyFolderCollection();
				
		//read gauges in the gauge folder
		if (mTMR.useInnerCircle) {
			mFC = jIO.getAllMyNeededFolders(myBaseFolder, myTiffs, "", true, false);   //"" put because no possibility for a filterTag implemented yet			
			//myGaugeFolder = jIO.findTheInnerCircle(myBaseFolder);		
			//myGauges = jIO.listInnerCircleFiles(myGaugeFolder, mTMR.filterTag);			
		}
		else mFC = jIO.getAllMyNeededFolders(myBaseFolder, myTiffs, "", false, false);   //"" put because no possibility for a filterTag implemented yet
				
		String myOutFolder;
		String myOutPath = null;		
		if (mTMR.useConstantThreshold == true) {
			
			mTMR = menu.showManualThresholdingDialog(mTMR);
			String thresholds = "";
			if (mTMR.minThreshold > 0) thresholds = "" + mTMR.minThreshold;
			if (mTMR.maxThreshold > 0 & thresholds.equalsIgnoreCase("")) thresholds = "" + mTMR.maxThreshold;
			
			myOutFolder = "ConstantThreshold" + thresholds;
			
		}
		else {myOutFolder = mTMR.myPrimaryMethod.toString();		
			if (mTMR.filterImages) myOutFolder = mTMR.filterTag + "_" + myOutFolder;
			if (mTMR.mySecondaryMethod != null) myOutFolder = mTMR.myPrimaryMethod.toString() + "_" + mTMR.mySecondaryMethod.toString();
			if (mTMR.setMaxgray2Wallgray) myOutFolder = myOutFolder + "_WallIsMaxgray";
		}		
		myOutPath = myBaseFolder + myOutFolder;						
		new File(myOutPath).mkdir();
		
		mFC.myOutFolder = myOutPath;
				
		//create folder for each horizontal cross-section for checking
		int[] eDC = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100};  //eDC stands for evaluation depth choices 
			
		//loop over 3D images
		int errorCounts = 0;
		for (i = 0 ; i < myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();
			
			//assign tiff file
			mFC.fileName = myTiffs[i];
			mFC = jIO.addCurrentFileInfo(mFC);
			
			//check if everything is in order
			if (mFC.somethingIsWrong) {
				IJ.error(mFC.eMsg);
				return;
			}
			
			//define evaluation depths
			int[] myZ = new int[eDC.length];
			for (int j = 0 ; j < eDC.length ; j++) {
				double checkZ = (double)eDC[j] * (double)mFC.nOfSlices / 100d;
				myZ[j] = (int)Math.round(checkZ);
			}
			if (myZ[0] == 0) myZ[0] = 1;
			
			try {
				
				//load image
				if (mTMR.save3DImage | mTMR.save4GeoDict) nowTiff = jIO.openTiff3D(mFC.nowTiffPath);
				else nowTiff = jIO.openTiff3DSomeSlices(mFC, myZ);
				
				//apply segmentation			
				outTiff = jIM.applyAChosenThresholdingMethod(nowTiff, mFC, mTMR, myZ);			 
			
				//save result
				if (mTMR.save3DImage == true) jIO.tiffSaver(myOutPath, myTiffs[i], outTiff[0]);
				if (mTMR.save4GeoDict == true) jIO.saveAsTiffStack4GeoDict(myOutPath, myTiffs[i], outTiff[0]);
				
				if (i == 0) jIO.writeStringIntoAsciiFile(myOutPath + "\\Path2Gauge.txt", mFC.myInnerCircleFolder + "\n" + mTMR.filterTag); 	//also save the gauge path..
				
			}
			catch(Exception e) {
					
				String eMsg = "Something went wrong when trying to segment column '" + mFC.fileName + "'!\n";
					
				errorCounts++;
					
				IJ.log(eMsg);	
			}	 			
		}
		
	}
	
}