package soilj;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import soilj.tools.InputOutput;
import soilj.tools.MenuWaiter;
import soilj.tools.MorphologyAnalyzer;
import soilj.tools.ObjectDetector;

import java.io.File;

public class SoilSurfaceFinder_ extends ImagePlus implements PlugIn  {

	// loads all 3D TIFFs in the specified folder and puts them straight
	// and aligns them as such that they stand 'perfectly' straight

	public void run(String arg) {

		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ObjectDetector jOD = new ObjectDetector();
		MenuWaiter menu = new MenuWaiter();		
					
		// init variables
		int i;
	
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();	
		ImagePlus outTiff = new ImagePlus();
		
		//show menu for SoilSurfaceFinder
		MenuWaiter.SurfaceFinderReturn mSFR = menu.new SurfaceFinderReturn();
		mSFR = menu.showSurfaceFinderMenu();
		
		//read base folder and number of 3D images
		String myBaseFolder;
		if (mSFR.binary == true) myBaseFolder = jIO.chooseAFolder("Please choose the folder with your binarized images!");
		else myBaseFolder = jIO.chooseAFolder("Please choose the folder with your grey-scale images!");
		String[] myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
		
		//read gauges in the gauge folder
		String myGaugeFolder = jIO.findTheInnerCircle(myBaseFolder);		
		String[] myGauges = jIO.listGaugeFiles(myGaugeFolder, "");
		
		//create output path for top surface
		String myOutFolder = "SurfaceOfColumn";
		String myOutPath = myBaseFolder + "\\" + myOutFolder;
		new File(myOutPath).mkdir();

		//init surface statistics output
		MorphologyAnalyzer.SurfaceStatistics[] mSS = new MorphologyAnalyzer.SurfaceStatistics[myTiffs.length]; 
		
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  //myTiffs.length		

			//try to free up some memory
			System.gc();
			
			//load image
			nowTiff = jIO.openTiff3D(myBaseFolder + "\\" + myTiffs[i]);	
					
			//apply segmentation			
			outTiff = jOD.findSoilSurface(nowTiff, myGauges[i], mSFR);			
						
			//analyze surfaces
			mSS[i] = jOD.extractSurfaceStatistics(myOutPath, myTiffs[i], outTiff, myGauges[i], nowTiff.getNSlices());
			
			//save result
			jIO.tiffSaver(myOutPath, myTiffs[i], outTiff);			
		}
		
		//save statistics
		jIO.saveSurfaceStatistics(myOutPath, myTiffs, mSS);
		
	}
}