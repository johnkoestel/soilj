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
import SoilJ.tools.InputOutput.MyFileCollection;

import java.io.File;

/** 
 * NormalizeTheseImages is a SoilJ plugin that standardizes the greyscale of an image with respect to 
 * grey value quantiles of horizontal slices and/or the grey value of the columns wall.
 * 
 * @author John Koestel
 *
 */

public class CalibrateGrayValues_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {
		
		String pathSep = "/";
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = CalibrateGrayValues_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ImageManipulator jIM = new ImageManipulator();	
		MenuWaiter menu = new MenuWaiter();
		
		MenuWaiter.CalibrationReferences myNR;
		
		// init variables
		int i;
		
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();		
		ImagePlus outTiff = new ImagePlus();
		
		//ask for threshold choice
		myNR = menu.showCalibrationMenu();
		if (myNR == null) return;

		//read file or files
		InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your image data");
		
		//add inner circle folder location if necessary
		if (myNR.useInnerCircle) {
			mFC = jIO.addInnerCircleFolder(mFC);
		}
		
		String myOutFolder;
		String myOutPath = null;		
		myOutFolder = "Standard_" + myNR.lowerTag + "And" + myNR.upperTag;		
		myOutFolder.replace(',', '.');
		myOutPath = mFC.myBaseFolder + pathSep + myOutFolder;
		new File(myOutPath).mkdir();
		mFC.myOutFolder = myOutPath;
				
		//loop over 3D images
		int errorCounts = 0;
		for (i = 0 ; i < mFC.myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();
			
			//assign tiff file
			mFC.fileName = mFC.myTiffs[i];
			mFC = jIO.addCurrentFileInfo(mFC);
			
			//check if everything is in order
			if (mFC.somethingIsWrong) {
				IJ.error(mFC.eMsg);
				return;
			}
			
			//try {
				//load image
				nowTiff = jIO.openTiff3D(mFC.myBaseFolder + pathSep + mFC.myTiffs[i]);			
				
				//select the correct gauge and surface files			
				String nowGauge = "";
				if (myNR.useInnerCircle) {
					int[] myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC);
					nowGauge = mFC.myInnerCircleFiles[myGandS[0]];
				}				
				
				//apply segmentation			
				outTiff = jIM.calibrateGrayValues(nowTiff, nowGauge, myNR, myOutPath);			
			
				//save result
				jIO.tiffSaver(myOutPath, mFC.myTiffs[i], outTiff);
			//}
			//catch(Exception e) {
				
			//	String eMsg = "Something went wrong when trying to standardize the gray values of column '" + mFC.fileName + "'!\n";
				
			//	errorCounts++;
				
			//	IJ.log(eMsg);	
			//}		
		}
		
		if (errorCounts > 0) {
			
			String eMsg = "\n\nThis crash was probably do to a bug in one of the SoilJ modules\n\n";
			eMsg += "Please feel free to try contacting John (john.koestel@slu.se).\n\n";
			eMsg += "Good luck and thank you for using SoilJ!";			
			
			IJ.log(eMsg);
		}		
		
		//finally also save the gauge path..
		if (myNR.useInnerCircle) jIO.writeStringIntoAsciiFile(myOutPath + pathSep + "Path2Gauge.txt", mFC.myInnerCircleFolder);
		
	}
	
}