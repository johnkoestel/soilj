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
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.io.FileInfo;
import ij.io.Opener;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import SoilJ.tools.DisplayThings;
import SoilJ.tools.ImageManipulator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.MorphologyAnalyzer;
import SoilJ.tools.ObjectDetector;
import SoilJ.tools.InputOutput.MyFolderCollection;
import SoilJ.tools.MorphologyAnalyzer.ProfileStatistics;
import SoilJ.tools.ObjectDetector.ColCoords3D;

import java.io.File;

import org.apache.commons.math3.stat.StatUtils;

/** 
 * FindColumnOutLines is a SoilJ plugin for the automatized soil column outline detection
 * 
 * @author John Koestel
 *
 */

public class FindColumnOutlines_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ObjectDetector jOD = new ObjectDetector();
		MenuWaiter menu = new MenuWaiter();
		DisplayThings disp = new DisplayThings();
				
		// init variables
		int i;
					
		//get details on PVC column finding
		MenuWaiter.ColumnFinderMenuReturn jCFS = menu.new ColumnFinderMenuReturn(); 
		jCFS = menu.showColumnFinderDialog();
		if (jCFS == null) return;
		
		//read base folder and number of 3D images
	    String[] myTiffs = null;String myBaseFolder = null;
	    while (myTiffs == null) {
	    	myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
			if (myBaseFolder == null) return;
			myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
			if (myTiffs == null) {
				IJ.error("Please choose a folder with TIFF images or cancel. Thank you.");	
			}
		}
		
		//collect all file informations
		InputOutput.MyFolderCollection mFC = jIO.new MyFolderCollection();
		mFC = jIO.getAllMyNeededFolders(myBaseFolder, myTiffs, "", false, false);   //"" put because no possibility for a filterTag implemented yet
				
		//create output folders and store information of them
		String myOutFolderName = "WallCoordinateIdentified";
		
		//create output folder
		String mySubBaseFolder = jIO.getTheFolderAbove(myBaseFolder);
		String myOutFolder = mySubBaseFolder + "\\" + myOutFolderName;
		new File(myOutFolder ).mkdir();
		mFC.myOutFolder = myOutFolder;
		
		//create innerCircle Folder
		String myInnerCircleFolder = myOutFolder + "\\InnerCircle";
		new File(myInnerCircleFolder).mkdir();
		mFC.myInnerCircleFolder = myInnerCircleFolder;
		
		//loop over 3D images
		int errorCounts = 0;
		for (i = 0 ; i < myTiffs.length ; i++) {  				
	
			//assign tiff file
			mFC.fileName = myTiffs[i];
			mFC = jIO.addCurrentFileInfo(mFC);
			
			//check if everything is in order
			if (mFC.somethingIsWrong) {
				IJ.error(mFC.eMsg);
				return;
			}
			
			//try {			
				//load a stack of sample images			
				InputOutput.SampleTiffWrapper sTW = jIO.assembleRepresentativeSample(mFC);
				
				//re-determine the column outlines			
				boolean look4PreciseCoords = true;
				ObjectDetector.ColCoords3D prelimCC = jOD.findOrientationOfPVCOrAluColumn(sTW.samTiff, jCFS, look4PreciseCoords);
				
				//find some outlines seriously
				ObjectDetector.ColCoords3D samCoords = jOD.findColumnWalls3D(sTW.samTiff, prelimCC, jCFS, sTW.samSlices);
					
				//find upper end of column				
				ObjectDetector.ColCoords3D topCoords = null;
				if (jCFS.try2FindColumnTopAndBottom) topCoords = jOD.findColumnsTop(mFC, samCoords, jCFS);
				else topCoords = jOD.findClosestXYSlice2Top(mFC, samCoords, jCFS);
				
				//find lower end of column and create complete set of slice coordinates		
				ObjectDetector.ColCoords3D roughCoords = null;
				if (jCFS.try2FindColumnTopAndBottom) roughCoords = jOD.findColumnsBottom(mFC, topCoords, jCFS);
				else roughCoords = jOD.findClosestXYSlice2Bottom(mFC, topCoords, jCFS);				
				
				//find bevel
				ObjectDetector.ColCoords3D bevCoords = jOD.new ColCoords3D();
				if (jCFS.hasBevel) bevCoords = jOD.findColumnsBevel(mFC, roughCoords, jCFS);
				else bevCoords = roughCoords;
				
				//interpolate column outlines for non-samples slices
				ObjectDetector.ColCoords3D colCoords = jOD.imputeMissingLayers(bevCoords, jCFS);
									
				//save inner circle file
				String myInnerCirclePath = mFC.myInnerCircleFolder + "\\Gauge" + mFC.colName + ".txt";
				jIO.writeInnerCircleVer1(myInnerCirclePath, colCoords);	
				
				//load complete image of soil column	
				int[] loadSlices = new int[mFC.nOfSlices];
				for (int j = 0 ; j < mFC.nOfSlices ; j++) loadSlices[j] = j;
				ImagePlus myTiff0 = jIO.openTiff3DSomeSlices(mFC, loadSlices);
				
				//create and save Z-projections depicting the outlines of the detected wall
				disp.displayColumnOutlinesByZ(myTiff0, colCoords, mFC);
				
				//cut and saveTiff
				ImagePlus myTiff = new ImagePlus();
				ImageStack outStack = new ImageStack(mFC.nowWidth, mFC.nowHeight);
				for (int j = colCoords.topOfColumn ; j < colCoords.bottomOfColumn ; j++) {
					if (colCoords.topOfColumn == 0) myTiff0.setPosition(j + 1);
					else myTiff0.setPosition(j);
					ImageProcessor nowIP = myTiff0.getProcessor();
					outStack.addSlice(nowIP);
				}
				myTiff.setStack(outStack);
				jIO.tiffSaver(mFC.myOutFolder, mFC.fileName, myTiff);
				
				//free memory
				myTiff.flush();myTiff0.flush();
				IJ.freeMemory();IJ.freeMemory();
			//}
				
/*			catch (Exception e) {
				
				String eMsg = "Something went wrong when trying to find the wall of column '" + mFC.fileName + "'!\n";
				eMsg += e.getCause() + e.getMessage() + "\n\n";
				
				errorCounts++;
				
				IJ.log(eMsg);	
				e.printStackTrace();
				IJ.log("\n\n");
				
			}	*/		
		}
		
		if (errorCounts > 0) {
			
			String eMsg = "\n\nThis plugin requires images on which the entire outer perimeter of the column wall\n";		
			eMsg += "is visible. If the wall is only partially visible, the plugin will probably not work.\n\n";		
			eMsg += "In the case that the contrast between wall and canvas gray-values is very small, the plugin is\n";
			eMsg += "doomed to fail.\n\n";
			eMsg += "In the case that the image is very noisy, use a filter to reduce the image noise\n";
			eMsg += "Suitable filters are a non-local means filter or, alternatively, a median filter \n";			
			eMsg += "followed by an unsharp mask.\n\n";
			if (jCFS.isAlu) {
				eMsg += "In the case that the contrast between column wall and soil is not very strong, \n";
				eMsg += "consider hacking in the wall thickness manually\n.";
				if(jCFS.hasBevel) {
					eMsg += "The location of the bevel should be found nevertheless.\n\n";
				}
				else eMsg += "\n";		
			}
			eMsg += "Thank you for using SoilJ!";			
			
			IJ.log(eMsg);
		}		
		
		IJ.freeMemory();IJ.freeMemory();
	}
		
}

