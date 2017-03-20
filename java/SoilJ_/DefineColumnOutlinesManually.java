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
import ij.Prefs;
import ij.plugin.PlugIn;
import SoilJ.tools.ImageManipulator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.ObjectDetector;

import java.io.File;

/** 
 * DefineColumnOutlinesManually is a SoilJ plugin that is supposed to one fine day contain a manual column outline detection.
 * 
 * @author John Koestel
 *
 */

public class DefineColumnOutlinesManually extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = DefineColumnOutlinesManually.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//fix memory settings		
		Prefs.keepUndoBuffers=false;		
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ObjectDetector jOD = new ObjectDetector();
		ImageManipulator jIM = new ImageManipulator();
		MenuWaiter menu = new MenuWaiter();
						
		// init variables
		int i;
		
		//construct image related objects
		ImagePlus upTiff;  //perfectly upright columns
					
		//get details on PVC column finding
		MenuWaiter.ColumnFinderMenuReturn jCFS = menu.new ColumnFinderMenuReturn(); 
		jCFS = menu.showColumnStraightenerMenu();
		
		//read base folder and number of 3D images
		String myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
		String[] myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
		String myTiffName;
		String myOutFolder = "StraightAndCentered";
						
		//create output folder
		String mySubBaseFolder = jIO.getTheFolderAbove(myBaseFolder);
		String myOutPath = mySubBaseFolder + "\\" + myOutFolder;
		new File(myOutPath).mkdir();
		
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  

			//assign current TIFF
			myTiffName = myTiffs[i];
			ImagePlus nowTiff = jIO.openTiff3D(myBaseFolder + "\\" + myTiffName);
			
			IJ.freeMemory();IJ.freeMemory();
			
			//try {
				//find the inner edge of the PVC column
				boolean look4PreciseCoords = false;
				ObjectDetector.ColCoords3D prelimCC = jOD.findOrientationOfPVCOrAluColumn(nowTiff, jCFS, look4PreciseCoords);
				
				IJ.freeMemory();IJ.freeMemory();	
						
				//check if column is already upright					
				upTiff = jIM.putColumnUprightInCenter(nowTiff, prelimCC, jCFS);
				
				nowTiff.flush();
				IJ.freeMemory();IJ.freeMemory();	
			
				//save the results
				jIO.tiffSaver(myOutPath, myTiffName, upTiff);				
				
			//}
			//catch(Exception e){}
			
			upTiff.flush();
			IJ.freeMemory();IJ.freeMemory();IJ.freeMemory();IJ.freeMemory();					
						
		}
	}
}