package SoilJ;

import ij.ImagePlus;
import ij.plugin.PlugIn;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MorphologyAnalyzer;

import java.io.File;

public class X2_SoilBulkVolumeAnalyzer_ extends ImagePlus implements PlugIn  {

	// loads all 3D TIFFs in the specified folder and puts them straight
	// and aligns them as such that they stand 'perfectly' straight

	public void run(String arg) {

		//construct biggish objects
		InputOutput jIO = new InputOutput();
		MorphologyAnalyzer morph = new MorphologyAnalyzer();
				
		// init variables
		int i;
		
		//construct image related objects
		ImagePlus soilSurface = new ImagePlus();
		
		//read base folder and number of 3D images
		String myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
		String[] myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
		
		//read soil top surface images
		String myTopSurfaceFolder = myBaseFolder + "SurfaceOfColumn";
		String[] myTSs = jIO.listTiffsInFolder(new File(myTopSurfaceFolder));
		
		//read gauges in the gauge folder
		InputOutput.MyMemory myMem = jIO.readBrainFile(myBaseFolder + "Path2Gauge.txt");
		String myGaugeFolder = myMem.gaugeFolder;
		String[] myGauges = jIO.listGaugeFiles(myGaugeFolder, myMem.filterTag);	
			
		//create output paths
		String myPreOutFolder;		
		myPreOutFolder = "WholeCol";		
		String myPreOutPath = myBaseFolder + myPreOutFolder;
		new File(myPreOutPath).mkdir();
		
		//init array for saving the bulk volumes
		int[][] myBulkVolumes = new int[myGauges.length][3]; 
		
		//loop over 3D images
		for (i = 0 ; i < myGauges.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();
			
			//select the correct gauge and surface files			
			int[] myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(myTiffs[i], myGauges, null);
					
			//load surfaces
			soilSurface = jIO.openTiff3D(myTopSurfaceFolder + "\\" + myTSs[myGandS[1]]);
							
			//apply segmentation	
			int[] nowVolume = new int[3];
			nowVolume = morph.calculateTheBulkSoilVolume(myPreOutPath, myGauges[myGandS[0]], soilSurface);
			for (int j = 0 ; j < 3 ; j++) myBulkVolumes[i][j] = nowVolume[j]; 			
			
		}
		
		//finally also save dimensions of cut-out volume..
		jIO.writeStringListIntoAsciiFile(myPreOutPath + "\\Stats\\BulkVolumes.txt", myTiffs, myBulkVolumes);
		
	}
}