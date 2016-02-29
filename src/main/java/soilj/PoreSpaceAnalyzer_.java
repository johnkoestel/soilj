package soilj;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import soilj.tools.InputOutput;
import soilj.tools.MenuWaiter;
import soilj.tools.MorphologyAnalyzer;
import soilj.tools.RoiHandler;

import java.io.File;

public class PoreSpaceAnalyzer_ extends ImagePlus implements PlugIn  {

	// loads all 3D TIFFs in the specified folder and puts them straight
	// and aligns them as such that they stand 'perfectly' straight

	public void run(String arg) {

		//construct biggish objects
		InputOutput jIO = new InputOutput();
		MenuWaiter menu = new MenuWaiter();
		RoiHandler roi = new RoiHandler();
		MorphologyAnalyzer morph = new MorphologyAnalyzer();
		
		// init variables
		int i;
		
		//shall I cut away something?
		MenuWaiter.PoreSpaceAnalyserOptions mPSA = menu.showPoreSpaceAnalyzerMenu();
		InputOutput.MyFolderCollection mFC = jIO.new MyFolderCollection();
		
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();
		
		//read base folder and number of 3D images
		String myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
		String[] myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
					
		//create output paths
		String myPreOutFolder = "";		
		if (mPSA.choiceOfRoi.equals("RealSample")) {			
			
			if (mPSA.cutAwayFromTop > 0 | mPSA.cutAwayFromWall > 0 | mPSA.heightOfRoi > 0) myPreOutFolder = "SubSam";
			else myPreOutFolder = "WholeCol";	
			
			mFC = jIO.getAllMyNeededFolders(myBaseFolder, myTiffs, "", mPSA.hasSurfaceFiles);   //"" put because no possibility for a filterTag implemented yet
			
		} else {
			
			if (mPSA.choiceOfRoi.equals("Cuboid")) myPreOutFolder = "Cuboid";
			if (mPSA.choiceOfRoi.equals("Cylinder")) myPreOutFolder = "Cylindrical";
			if (mPSA.choiceOfRoi.equals("Everything!")) myPreOutFolder = "EntireImage";
			
			//add basic folders to folder collection
			mFC.myBaseFolder = myBaseFolder;
			mFC.myTiffs = myTiffs;
			
		}
		String myPreOutPath = myBaseFolder + myPreOutFolder;
		new File(myPreOutPath).mkdir();		
		mFC.myPreOutFolder = myPreOutPath; 	//also add it to the folder collection
		
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();
			
			//load image
			nowTiff = jIO.openTiff3D(myBaseFolder + "\\" + myTiffs[i]);	
			
			ImagePlus[] toAnalyseTiffs = roi.prepareDesiredRoi(i, mFC, nowTiff, mPSA);
												
			//apply segmentation	
			morph.tailoredPoreSpaceAnalyses(i, mFC, toAnalyseTiffs, mPSA);			
			
		}
		
	}
}