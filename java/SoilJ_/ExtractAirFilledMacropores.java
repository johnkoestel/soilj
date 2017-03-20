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
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.ObjectDetector;

import java.io.File;

/** 
 * ExtractAirFilledMacropores is a SoilJ plugin for extracting the air-filled pore space 
 * from a pore-diameter image (thickness image) for soil columns with a surface exposed to te atmosphere.
 * The plugin is based in the capillary equation. The soil is assumed to be in equilibrium with a specific 
 * matrix potential defined at the lower boundary of the column.
 * The plugin does not take hysteretic effects into account.
 * The binary images of the air filled pore space are saved in a newly created folder.
 * 
 * @author John Koestel
 *
 */

public class ExtractAirFilledMacropores extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		//set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = ExtractAirFilledMacropores.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ObjectDetector jOD = new ObjectDetector();
		MenuWaiter menu = new MenuWaiter();
		
		// init variables
		int i;
		MenuWaiter.AirFilledPoresFinderMenu aPF = menu.new AirFilledPoresFinderMenu();
		aPF = menu.showAPFinderMenu();
		if (aPF == null) return;
		
		//construct image related objects
		ImagePlus nowTiff;  //input
		ImagePlus outTiff;  //output
		
		//read base folder and number of 3D images
	    String[] myTiffs = null;String myBaseFolder = null;
	    while (myTiffs == null) {
	    	myBaseFolder = jIO.chooseAFolder("Please choose the folder with thickness images.");
			if (myBaseFolder == null) return;
			myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
			if (myTiffs == null) {
				IJ.error("Please choose a folder with TIFF images or cancel. Thank you.");	
			}
		}
		String myTiffName;
		String myOutFolder = "AirFilledPoresAt_" + (int)Math.round(aPF.matrixPotentialAtColumnBottomInCM) + "_cmTension";
						
		//create output folder
		String mySubBaseFolder = jIO.getTheFolderAbove(myBaseFolder);
		String myOutPath = mySubBaseFolder + "\\" + myOutFolder;
		new File(myOutPath).mkdir();
				
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  

			//assign current TIFF
			myTiffName = myTiffs[i];
			nowTiff = jIO.openTiff3D(myBaseFolder + "\\" + myTiffName);
			
			//do the cutting
			outTiff = jOD.extractAirFilledPores(nowTiff, aPF);
			
			//save the results
			jIO.tiffSaver(myOutPath, myTiffName, outTiff);			
			
		}
	}
}