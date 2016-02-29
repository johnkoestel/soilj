package SoilJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import SoilJ.tools.ImageManipulator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;

import java.io.File;

public class T0_BeamDeHardening_ extends ImagePlus implements PlugIn  {

	// loads all 3D TIFFs in the specified folder and puts them straight
	// and aligns them as such that they stand 'perfectly' straight

	public void run(String arg) {

		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ImageManipulator jIM = new ImageManipulator();
		MenuWaiter menu = new MenuWaiter();
						
		// init variables
		int i;
		
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();
		ImagePlus outTiff = new ImagePlus();
		
		//ask for threshold choice
	    MenuWaiter.BeamDeHardeningReturn mBDH = menu.showBeamDeHardeningMenu();
			    
		//read base folder and number of 3D images
		String myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
		String[] myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
			
		//read gauges in the gauge folder
		InputOutput.MyFolderCollection mFC = jIO.getAllMyNeededFolders(myBaseFolder, myTiffs, "", false);

		//save and create output folder		
		String myOutFolder;
		String myOutPath = null;
		myOutFolder = "BeamDeHardened";			
		myOutPath = myBaseFolder + myOutFolder;
		new File(myOutPath).mkdir();		
				
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();		
			
			//load image
			nowTiff = jIO.openTiff3D(myBaseFolder + "\\" + myTiffs[i]);	
			
			//select the correct gauge and surface files			
			int[] myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(myTiffs[i], mFC.myGaugePaths, null);
			
			outTiff = jIM.beamDeHardenThis(nowTiff, mFC.myGaugePaths[myGandS[0]], mBDH);
			
			//save result
			jIO.tiffSaver(myOutPath, myTiffs[i], outTiff);
		}
		
	}
	
}