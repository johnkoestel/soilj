package soilj;

import ij.ImagePlus;
import ij.plugin.PlugIn;
import soilj.tools.ImageManipulator;
import soilj.tools.InputOutput;
import soilj.tools.MenuWaiter;

import java.io.File;

public class X1_SoilClipper_ extends ImagePlus implements PlugIn  {

	// loads all 3D TIFFs in the specified folder and puts them straight
	// and aligns them as such that they stand 'perfectly' straight

	public void run(String arg) {

		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ImageManipulator jIM = new ImageManipulator();
		MenuWaiter menu = new MenuWaiter();		 
						
		// init variables
		MenuWaiter.ClipperMenuReturn mSCM;
		int i;
		
		//construct image related objects	
		InputOutput.MyFolderCollection mFC = jIO.new MyFolderCollection();
		ImagePlus outTiff = new ImagePlus();
		
		//read base folder and number of 3D images
		String myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
		String[] myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
		ImagePlus nowTiff = new ImagePlus();
		
		//ask for threshold choice
		mSCM = menu.showClipperDialog();

		//init output folder
		String choice1; if (mSCM.isCylinder == true) choice1 = "CylHor"; else choice1 = "RectHor"; 
		String choice2; if (mSCM.referenceIsSoilSurface == true) choice2 = "_RefSoilSurf"; else choice2 = "_RefIsTop";
		String myOutFolder = "ROI_" + choice1 + mSCM.clipFromInnerPerimeter + choice2 + mSCM.startAtSlice + "_Height" + mSCM.heightOfROI;
		String myOutPath = null;				
		myOutPath = myBaseFolder + myOutFolder;		
				
		new File(myOutPath).mkdir();
		
		boolean hasSurfaceFiles = false; if (mSCM.referenceIsSoilSurface == true) hasSurfaceFiles = true;
		if (mSCM.isCylinder == true) mFC = jIO.getAllMyNeededFolders(myBaseFolder, myTiffs, "", hasSurfaceFiles);   //"" put because no possibility for a filterTag implemented yet
				
		//also make a new Inner Circle folder for the updated inner circle files		
		String myNewGaugeFolder = "InnerCircle";
		String myNewGaugePath = myOutFolder + "\\" + myNewGaugeFolder;
		new File(myNewGaugePath).mkdir();
		
		//add basic folders to folder collection
		mFC.myBaseFolder = myBaseFolder;
		mFC.myTiffs = myTiffs;		
		mFC.myOutFolder = myOutFolder;
		
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();
			
			//load image
			nowTiff = jIO.openTiff3D(myBaseFolder + "\\" + myTiffs[i]);	
					
			//clip	 
			outTiff = jIM.clipImage(i, nowTiff, mFC, mSCM); 
			
			//save result
			jIO.tiffSaver(myOutPath, myTiffs[i], outTiff);
			jIO.writeCutEllipsoidMaskFile(i, myNewGaugePath, mFC, mSCM);			
			
		}
		
	}
	
}