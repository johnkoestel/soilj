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
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.MorphologyAnalyzer;
import SoilJ.tools.ObjectDetector;
import SoilJ.tools.RoiHandler;

import java.io.File;

/** 
 * PoreSpaceAnalyzer is a SoilJ plugin that calculates several morphologic properties from 
 * binary images. To a large part it makes use of plugins collected in the BoneJ library:
 * Doube M, Kłosowski MM, Arganda-Carreras I, Cordeliéres F, Dougherty RP, Jackson J, Schmid B, Hutchinson JR, Shefelbine SJ. (2010) BoneJ: free and extensible bone image analysis in ImageJ. Bone 47:1076-9. doi: 10.1016/j.bone.2010.08.023 
 * www.bonej.org
 * 
 * The results of the analyses are written into image and ASCII files. 
 * 
 * @author John Koestel
 *
 */

public class PoreSpaceAnalyzer_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = PoreSpaceAnalyzer_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();
		ObjectDetector jOD = new ObjectDetector();
		MenuWaiter menu = new MenuWaiter();
		RoiHandler roi = new RoiHandler();
		MorphologyAnalyzer morph = new MorphologyAnalyzer();
		
		// init variables
		int i;
		
		//tell me what I should do!
		MenuWaiter.PoreSpaceAnalyzerOptions mPSA = menu.showPoreSpaceAnalyzerMenu();
		if (mPSA == null) return;
		
		InputOutput.MyFolderCollection mFC = jIO.new MyFolderCollection();
		
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();
		
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
					
		//create output paths
		String myPreOutFolder = "";		
		if (mPSA.choiceOfRoi.equals("RealSample")) {			
			
			if (mPSA.heightOfRoi > 0) myPreOutFolder = "IC_Top" + mPSA.cutAwayFromTop + "Height" + mPSA.heightOfRoi + "Wall" + mPSA.cutAwayFromWall;
			boolean isCut = false;
			if (mPSA.cutAwayFromBottom > 0 | mPSA.cutAwayFromTop > 0 | mPSA.cutAwayFromWall > 0 | mPSA.cutAwayFromCenter > 0) isCut = true; 
			if (isCut & mPSA.heightOfRoi > 0) myPreOutFolder = "IC_Tv" + mPSA.cutAwayFromTop + "Height" + mPSA.heightOfRoi + "Wall" + mPSA.cutAwayFromWall;
			if (isCut & mPSA.cutZPercent & mPSA.cutXYPercent) myPreOutFolder = 
					"IC_Tp" + mPSA.cutAwayFromTop + "Bp" + mPSA.cutAwayFromBottom +	"Wp" + mPSA.cutAwayFromWall + "Cp" + mPSA.cutAwayFromCenter;
			if (isCut & !mPSA.cutZPercent & mPSA.cutXYPercent) myPreOutFolder = 
					"IC_Tv" + mPSA.cutAwayFromTop + "Bv" + mPSA.cutAwayFromBottom +	"Wp" + mPSA.cutAwayFromWall + "Cp" + mPSA.cutAwayFromCenter;
			if (isCut & mPSA.cutZPercent & !mPSA.cutXYPercent) myPreOutFolder = 
					"IC_Tp" + mPSA.cutAwayFromTop + "Bp" + mPSA.cutAwayFromBottom +	"Wv" + mPSA.cutAwayFromWall + "Cv" + mPSA.cutAwayFromCenter;
			if (isCut & !mPSA.cutZPercent & !mPSA.cutXYPercent) myPreOutFolder = 
					"IC_Tv" + mPSA.cutAwayFromTop + "Bv" + mPSA.cutAwayFromBottom +	"Wv" + mPSA.cutAwayFromWall + "Cv" + mPSA.cutAwayFromCenter;
			if (!isCut) myPreOutFolder = "InnerCircleColumn";	
			
			boolean useInnerCircle = mPSA.choiceOfRoi.equalsIgnoreCase("RealSample");
			mFC = jIO.getAllMyNeededFolders(myBaseFolder, myTiffs, "", mPSA.useInnerCircleFiles, mPSA.includeSurfaceTopography);   //"" put because no possibility for a filterTag implemented yet
			
		} else {
			
			if (mPSA.choiceOfRoi.equals("Cuboid")) myPreOutFolder = "Cuboid_XL" + mPSA.cubeX1 + "XR" + mPSA.cubeX2 + "YL" + mPSA.cubeY1 + "YR" + mPSA.cubeY2 + "ZT" + mPSA.cubeZ1 + "ZB" + mPSA.cubeZ2;
			if (mPSA.choiceOfRoi.equals("Cylinder")) myPreOutFolder = "Cyl_X" + mPSA.cylX + "Y" + mPSA.cylY+ "ZT" + mPSA.cylZ1 + "ZB" + mPSA.cylZ2 + "R" + mPSA.cylRadius;
			if (mPSA.choiceOfRoi.equals("Everything!")) myPreOutFolder = "EntireImage";
			if (mPSA.choiceOfRoi.equals("TopOfEveryThing")) myPreOutFolder = "TopOfEntireImage";
			if (mPSA.choiceOfRoi.equals("BottomOfEveryThing")) myPreOutFolder = "BottomOfEntireImage";
			
			//add basic folders to folder collection
			mFC.myBaseFolder = myBaseFolder;
			mFC.myTiffs = myTiffs;
			
		}
		String myPreOutPath = myBaseFolder + myPreOutFolder;
		new File(myPreOutPath).mkdir();		
		mFC.myPreOutFolder = myPreOutPath; 	//also add it to the folder collection
		
		//also create a folder fo the cut surfaces
		if (mPSA.includeSurfaceTopography) {
			String myCutSurfacePath = mFC.myOutFolder + "//CutSurfaceFiles";
			new File(myCutSurfacePath).mkdir();
			mFC.myCutSurfaceFolder = myCutSurfacePath;
		}
		
		if (mPSA.includeSurfaceTopography) {
			mFC = jIO.getAllMyNeededFolders(mFC.myBaseFolder, mFC.myTiffs, "", mPSA.useInnerCircleFiles, mPSA.includeSurfaceTopography);
		}			
		
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory			
			IJ.freeMemory();IJ.freeMemory();		
			
			//assemble image for analyses
			mFC.fileName = myTiffs[i];
			mFC = jIO.addCurrentFileInfo8Bit(mFC);
			
			//check if soil surface should be taken into account
			int topSurface = 0;
			int botSurface = mFC.nOfSlices;
			if (mPSA.includeSurfaceTopography) {			
				MorphologyAnalyzer.SurfaceStatistics mST = jOD.extractSurfaceStatistics(mFC);
				topSurface = mST.medianElevation;
				botSurface = mST.medianIntrusion;
			}
			
			//check if there is something to cut away
			int startSlice = topSurface + mPSA.cutAwayFromTop;
			int stopSlice = botSurface - mPSA.cutAwayFromBottom;
			if (mPSA.cutZPercent) {
				startSlice = (int)Math.round((float)mPSA.cutAwayFromTop / 100f * (float)mFC.nOfSlices);
				stopSlice = mFC.nOfSlices - (int)Math.round((float)mPSA.cutAwayFromBottom / 100f * (float)mFC.nOfSlices);
			}
			if (mPSA.heightOfRoi > 0) stopSlice = startSlice + mPSA.heightOfRoi;
					
			//load file
			int[] colSlices = new int[stopSlice - startSlice];
			for (int j = 0 ; j < colSlices.length ; j++) colSlices[j] = startSlice + j;
			nowTiff = jIO.openTiff3DSomeSlices(mFC, colSlices);
			
			//cut image			
			ImagePlus[] toAnalyseTiffs = roi.prepareDesiredRoi(mFC, nowTiff, mPSA);					
						
			//try to free up some memory			
			IJ.freeMemory();IJ.freeMemory();		
			
			//apply analyzes
			morph.tailoredPoreSpaceAnalyses(i, mFC, toAnalyseTiffs, mPSA);			
			
		}
		
	}
}