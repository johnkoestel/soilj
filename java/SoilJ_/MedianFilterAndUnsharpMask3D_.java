package SoilJ_;

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

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import SoilJ.tools.ImageManipulator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;

import java.io.File;

/** 
 * MedianFilterAndUnsharpMask3D is a SoilJ plugin that applies a 3-D median filter and a 
 * 3-D unsharp mask. 
 * 
 * @author John Koestel
 *
 */

public class MedianFilterAndUnsharpMask3D_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = MedianFilterAndUnsharpMask3D_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ImageManipulator jIM = new ImageManipulator();
		MenuWaiter menu = new MenuWaiter();		
				
		// init variables
		int i;
		
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();
		ImagePlus outTiff = new ImagePlus();
		
		//open dialog box to query the filter settings
		MenuWaiter.MedianFilterAndUnsharpMaskReturn mMUS = menu.showMedianAndUsharpMaskMenu();
		if (mMUS == null) return;
		
		//load list of images to process
	    String[] myTiffs = null;String myBaseFolder = null;
	    while (myTiffs == null) {
	    	myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
			if (myBaseFolder == null) return;
			myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
			if (myTiffs == null) {
				IJ.error("Please choose a folder with TIFF images or cancel. Thank you.");	
			}
		}
		String myTiffName;
		String myOutFolder = "UnsharpMedian";
						
		//create output folder
		String mySubBaseFolder = jIO.getTheFolderAbove(myBaseFolder);
		String myOutPath = mySubBaseFolder + "\\" + myOutFolder;
		new File(myOutPath).mkdir();		
		
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  	//myTiffs.length

			//assign current TIFF
			myTiffName = myTiffs[i];
			
			//load TIFF and apply filters
			nowTiff = jIO.openTiff3D(myBaseFolder + "\\" + myTiffName);						
			outTiff = jIM.applyMedianFilterAndUnsharpMask(nowTiff, mMUS);	
			
			//save result
			jIO.tiffSaver(myOutPath, myTiffName, outTiff);
			
			//try to free up some memory
			IJ.freeMemory();IJ.freeMemory();			
		}
	}
}