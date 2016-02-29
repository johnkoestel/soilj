package soilj;

import ij.ImagePlus;
import ij.plugin.PlugIn;
import soilj.tools.InputOutput;

import java.io.File;

public class CombineTiffStack2Tiff_ extends ImagePlus implements PlugIn {

	// loads a stack of tiff images from the Offer sampling campaign (June 2013)
	// and aligns them as such that they stand 'perfectly' straight

	public void run(String arg) {
		
		//init objects		
		InputOutput jIO = new InputOutput();
	
		ImagePlus imgStack;
		
		String myBaseFolder;
		String mySubBaseFolder;
		String[] myFolders;
		String myOutFolder = "3D";
		String myOutPath;
		int i;
		
		//read base folder and number of 3D images
		myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
		myFolders = jIO.listFoldersInFolder(new File(myBaseFolder));
		
		//create output folder
		mySubBaseFolder = jIO.getTheFolderAbove(myBaseFolder);
		myOutPath = mySubBaseFolder + "\\" + myOutFolder;
		new File(myOutPath).mkdir(); 
		
		//go through each tiffstack and save it as one 3D tiff in the output folder
		for (i = 0 ; i < myFolders.length ; i++) {
			
			File myFolder = new File(myBaseFolder + "\\" + myFolders[i]);
			imgStack = jIO.openStack(myFolder.getAbsolutePath());
			jIO.tiffSaver(myOutPath, myFolders[i] + ".tif", imgStack);						
					
		}	
	}
}