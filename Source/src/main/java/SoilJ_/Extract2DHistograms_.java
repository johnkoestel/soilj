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
import SoilJ.tools.HistogramStuff;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.ObjectDetector;
import SoilJ.tools.RoiHandler;
import SoilJ.tools.InputOutput.MyFileCollection;

import java.io.File;

import org.apache.commons.math3.stat.StatUtils;

/** 
 * ExtractFreshOM is a SoilJ plugin aiming at segmenting the fresh soil organic matter and/or water phase in the column
 * and saves the respective binary images in a newly created folder.  
 *  
 * @author John Koestel
 */

public class Extract2DHistograms_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "/";
		
		//set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = Extract2DHistograms_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ObjectDetector jOD = new ObjectDetector();
		MenuWaiter menu = new MenuWaiter();
		RoiHandler roi = new RoiHandler();
		HistogramStuff hist = new HistogramStuff();
		
		// init variables
		int i;
		MenuWaiter.Extract2DHistogramOptions e2DH = menu.new Extract2DHistogramOptions();
		e2DH = menu.show2DHistogramExtractionMenu();
		if (e2DH == null) return;
		MenuWaiter.ROISelectionOptions mRSO = e2DH.mRSO;
	  	    
	    //select file or files
	  	InputOutput.MyFileCollection mFC = jIO.fileSelector("Please choose a file or folder with your image data");
	  	
	  	//create also select gradient files 
	  	String myGradFolder;
	    if (e2DH.calcGradientImage) {
	    	myGradFolder = mFC.myBaseFolder + pathSep + "Gradients";
	    	new File(myGradFolder).mkdir();		
	    }
	    if (!e2DH.calcGradientImage) {
	    	if (jIO.testIfFolderIsPresent(new File(mFC.myBaseFolder), "Gradients")) {
	    		myGradFolder = mFC.myBaseFolder + pathSep + "Gradients";
	    	}
	    	else {	    	
	    		myGradFolder = null;
				while (myGradFolder == null) {
					myGradFolder = jIO.chooseAFolder("Please choose the folder with your gradient image data");
					if (myGradFolder == null) return;
					
					//remember gradient folder
					mFC.myGradientFolder = myGradFolder;
				}		
	    	}
	    }		    
	    
		//add inner circle folder location if necessary
		if (mRSO.useInnerCircleFiles) {
			mFC = jIO.addInnerCircleFolder(mFC);
		}
		
		//add surface topography folder location if necessary
		if (mRSO.includeSurfaceTopography) {
			mFC = jIO.addSurfaceTopographyFolder(mFC);
		}
	    
	    //create output paths
  		String myPreOutFolder = "";  		
  		if (mRSO.choiceOfRoi.equals("RealSample")) {	
  			if (mRSO.heightOfRoi > 0) myPreOutFolder = "Top" + mRSO.cutAwayFromTop + "Height" + mRSO.heightOfRoi + "Wall" + mRSO.cutAwayFromWall;
	  			boolean isCut = false;
	  			if (mRSO.cutAwayFromBottom > 0 | mRSO.cutAwayFromTop > 0 | mRSO.cutAwayFromWall > 0 | mRSO.cutAwayFromCenter > 0) isCut = true; 
	  			if (isCut & mRSO.heightOfRoi > 0) myPreOutFolder = "Tv" + mRSO.cutAwayFromTop + "Height" + mRSO.heightOfRoi + "Wall" + mRSO.cutAwayFromWall;
	  			if (isCut & mRSO.cutZPercent & mRSO.cutXYPercent & mRSO.heightOfRoi == 0) myPreOutFolder = 
	  					"Tp" + mRSO.cutAwayFromTop + "Bp" + mRSO.cutAwayFromBottom +	"Wp" + mRSO.cutAwayFromWall + "Cp" + mRSO.cutAwayFromCenter;
	  			if (isCut & !mRSO.cutZPercent & mRSO.cutXYPercent & mRSO.heightOfRoi == 0) myPreOutFolder = 
	  					"Tv" + mRSO.cutAwayFromTop + "Bv" + mRSO.cutAwayFromBottom +	"Wp" + mRSO.cutAwayFromWall + "Cp" + mRSO.cutAwayFromCenter;
	  			if (isCut & mRSO.cutZPercent & !mRSO.cutXYPercent & mRSO.heightOfRoi == 0) myPreOutFolder = 
	  					"Tp" + mRSO.cutAwayFromTop + "Bp" + mRSO.cutAwayFromBottom +	"Wv" + mRSO.cutAwayFromWall + "Cv" + mRSO.cutAwayFromCenter;
	  			if (isCut & !mRSO.cutZPercent & !mRSO.cutXYPercent & mRSO.heightOfRoi == 0) myPreOutFolder = 
	  					"Tv" + mRSO.cutAwayFromTop + "Bv" + mRSO.cutAwayFromBottom +	"Wv" + mRSO.cutAwayFromWall + "Cv" + mRSO.cutAwayFromCenter;
	  			if (!isCut) myPreOutFolder = "InnerCircleColumn";	
	  			  			
	  			if (mRSO.includeSurfaceTopography) myPreOutFolder = "S_"+ myPreOutFolder;
	  			else myPreOutFolder = "C_"+ myPreOutFolder;
	  			mFC = jIO.getAllMyNeededFolders(mFC.myBaseFolder, mFC.myTiffs, "", mRSO.useInnerCircleFiles, mRSO.includeSurfaceTopography);   //"" put because no possibility for a filterTag implemented yet
	  			
	  		} 
  		else {
	  			
  			if (mRSO.choiceOfRoi.equals("Cuboid")) myPreOutFolder = "Cuboid_XL" + mRSO.cubeX1 + "XR" + mRSO.cubeX2 + "YL" + mRSO.cubeY1 + "YR" + mRSO.cubeY2 + "ZT" + mRSO.cubeZ1 + "ZB" + mRSO.cubeZ2;
  			if (mRSO.choiceOfRoi.equals("Cylinder")) myPreOutFolder = "Cyl_X" + mRSO.cylX + "Y" + mRSO.cylY+ "ZT" + mRSO.cylZ1 + "ZB" + mRSO.cylZ2 + "R" + mRSO.cylRadius;
  			if (mRSO.choiceOfRoi.equals("Everything!")) myPreOutFolder = "EntireImage";
  			if (mRSO.choiceOfRoi.equals("TopOfEveryThing")) myPreOutFolder = "TopOfEntireImage";
  			if (mRSO.choiceOfRoi.equals("BottomOfEveryThing")) myPreOutFolder = "BottomOfEntireImage";
 
  		}
  		
  		String myOutPath = mFC.myBaseFolder + pathSep + "Histo_" + myPreOutFolder;
  		new File(myOutPath).mkdir();		
  		mFC.myOutFolder = myOutPath; 	//also add it to the folder collection
  		
  		//save pathSep
  		mFC.pathSep = pathSep;
  		
  		//also create a folder for the cut surfaces  
		if (mRSO.includeSurfaceTopography) {  		
			String myCutSurfacePath = mFC.myPreOutFolder + pathSep + "CutSurfaceFiles";
			new File(myCutSurfacePath).mkdir();
			mFC.myCutSurfaceFolder = myCutSurfacePath;
		}
  
		
		//saveOutPut 2D histogram
		double[][] out2DHist = new double[256][256];   //b16 == 256 * 256
				

		
		//loop over 3D images			
		for (i = 0 ; i < mFC.myTiffs.length ; i++) {  
			//assign current TIFF
			mFC.fileName = mFC.myTiffs[i];
			mFC = jIO.addCurrentFileInfo(mFC);
			mFC.nowTiffPath = mFC.myBaseFolder + pathSep + mFC.fileName;
			ImagePlus nowTiff = jIO.loadCarvedTiff(mFC, mRSO);
			
			//load or create gradient Image
			ImagePlus gradTiff = new ImagePlus();	
			if (e2DH.calcGradientImage) {
				gradTiff = jOD.createGradientImage(mFC, jIO.openTiff3D(mFC));
			}
			try {
				mFC.nowTiffPath = mFC.myGradientFolder + pathSep + mFC.fileName;
				gradTiff = jIO.loadCarvedTiff(mFC, mRSO);					
			}
			catch(Exception e) {
				IJ.error("I did not find a gradient image to " + mFC.fileName + ".\nPlease re-run this plugin and check the option to calculate a gradient image.");
				return;
			}
			
			//cut out ROIs
			RoiHandler.ColumnRoi nowRoi = roi.prepareDesiredRoi(mFC, nowTiff, mRSO);				
			RoiHandler.ColumnRoi gradRoi = roi.prepareDesiredRoi(mFC, gradTiff, mRSO);		
			
			//try to free up some memory			
			IJ.freeMemory();IJ.freeMemory();	
			
			//do the cutting
			int[][] hist2D = hist.extract2DHistogram(nowRoi.nowTiff, gradRoi.nowTiff);
			double sumEntries = 0;
			for (int j = 0 ; j < 256 ; j++) sumEntries += hist2D[j][0]; 
			
			//add to overall histogram
			double[][] indieHist = new double[256][256];
			for (int x = 0 ; x < 256 ; x++) {
				for (int y = 0 ; y < 256 ; y++) {
					out2DHist[x][y] += (double)hist2D[x][y] / (2 * sumEntries);
					indieHist[x][y] = (double)hist2D[x][y] / (2 * sumEntries);
				}
			}		
			//save the results
			mFC.nowTiffPath = mFC.myOutFolder + pathSep + mFC.fileName;
			jIO.save2DHistogramAsTiff(mFC, indieHist);
		}		
			
		//try to clear memory
		IJ.freeMemory();IJ.freeMemory();	
		
	}
}