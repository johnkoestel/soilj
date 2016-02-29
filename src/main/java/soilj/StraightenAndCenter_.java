package SoilJ;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.plugin.PlugIn;
import SoilJ.tools.ImageManipulator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.ObjectDetector;

import java.io.File;

public class StraightenAndCenter_ extends ImagePlus implements PlugIn  {

	// loads all 3D TIFFs in the specified folder and puts them straight
	// and aligns them as such that they stand 'perfectly' straight

	public void run(String arg) {

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
				ObjectDetector.ColumnCoordinates prelimCC = jOD.findOrientationOfPVCOrAluColumn(nowTiff, jCFS);
				
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