package soilj;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import soilj.tools.ArtificialPoreNetworkCreator;
import soilj.tools.InputOutput;
import soilj.tools.MenuWaiter;

import java.io.File;

public class X9_GenerateRandomPoreClusters_ extends ImagePlus implements PlugIn  {

	// loads all 3D TIFFs in the specified folder and puts them straight
	// and aligns them as such that they stand 'perfectly' straight

	public void run(String arg) {

		//construct biggish objects
		InputOutput jIO = new InputOutput();	
		MenuWaiter menu = new MenuWaiter();		
		ArtificialPoreNetworkCreator creator = new ArtificialPoreNetworkCreator();
		
		MenuWaiter.RandomClusterGenerator mRCG = menu.new RandomClusterGenerator();
				
		// init variables
		int i;
		
		//get information on what shall be done..
		mRCG = menu.showRandomClusterGeneratorMenu();
		
		//read base folder and number of 3D images
		String myBaseFolder = jIO.chooseAFolder("Please choose the folder into which you want to save your images");		
		String myOutFolder = myBaseFolder + mRCG.shape + mRCG.domainX + "x" + mRCG.domainY + "x" + mRCG.domainZ;
		new File(myOutFolder).mkdir();		
		
		//in case that the porosities should be taken from a list..
		if (mRCG.mode.equalsIgnoreCase("predefinedList")) {
			String myPorosityListFile = jIO.chooseAnAsciiFile("Please choose the file that contains the porosities..",myBaseFolder,"*.*");
			double[] myPorosities = jIO.readSingleColumnDoublesFromAscii(myPorosityListFile);	
			mRCG.porosityList = myPorosities;
		}
		
		//loop over 3D images
		for (i = 0 ; i < mRCG.numOfCopies ; i++) {  //myTiffs.length
			
			//apply segmentation			
			creator.createASetOfRandomFields(myOutFolder, i, mRCG);
			
			IJ.freeMemory();IJ.freeMemory();
			
		}
		
	}
}