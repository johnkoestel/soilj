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
import SoilJ.tools.MorphologyAnalyzer;
import SoilJ.tools.RoiHandler;

import java.io.File;

/** 
 * ExtractPercolatingPorosity is a SoilJ plugin that extracts the percolating phase from a binary image.
 * 
 * @author John Koestel
 *
 */

public class ExtractPercolatingPorosity extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "/";
		
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = PoreSpaceAnalyzer_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();
		//MenuWaiter menu = new MenuWaiter();
		//RoiHandler roi = new RoiHandler();
		MorphologyAnalyzer morph = new MorphologyAnalyzer();
		
		// init variables
		int i;
		
		InputOutput.MyFileCollection mFC = jIO.new MyFileCollection();
		
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();
		
		//read base folder and number of 3D images
	    String[] myTiffs = null;String myBaseFolder = null;
	    while (myTiffs == null) {
	    	myBaseFolder = jIO.chooseAFolder("Please choose the folder with your cluster label images");
			if (myBaseFolder == null) return;
			myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
			if (myTiffs == null) {
				IJ.error("Please choose a folder with TIFF images or cancel. Thank you.");	
			}
		}

     	mFC.myBaseFolder = myBaseFolder;
		mFC.myTiffs = myTiffs;
		
		String myPreOutPath = myBaseFolder;	
		mFC.myPreOutFolder = myPreOutPath; 	//also add it to the folder collection
		
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();
			
			//load image
			nowTiff = jIO.openTiff3D(myBaseFolder + pathSep + myTiffs[i]);	
												
			//apply segmentation	
			morph.extractPoresizeDistro(i, mFC, nowTiff);			
			
		}
		
	}
}