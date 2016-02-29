package soilj.tools;

import java.awt.Color;
import java.util.Random;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.OvalRoi;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import soilj.tools.InputOutput;

public class ArtificialPoreNetworkCreator implements PlugIn {

	
	public void run(String arg) {
				//ok, this is not needed..
	}

	public void createASetOfRandomFields(String myPreOutPath, int randomSetNumber, MenuWaiter.RandomClusterGenerator mRCG) {
	
		InputOutput jIO = new InputOutput();
		Random randomNumberGenerator = new Random();	
		
		//init some variables
		int[] xyz = new int[]{mRCG.domainX, mRCG.domainY, mRCG.domainZ};		
		int numberOfPoreVoxels = 0;		
		int fieldsToCreate = 1;
		if (mRCG.mode.equalsIgnoreCase("predefinedList")) fieldsToCreate = mRCG.porosityList.length;
		int[] numberOfPoreVoxelList = new int[fieldsToCreate];
		
		//init random field generator following the selected mode
		double myPorosity = 0;
		if (mRCG.mode.equalsIgnoreCase("Gaussian")) {
			double pc = mRCG.pc;
			double relSTD = mRCG.standardDeviation;
			myPorosity = pc - relSTD * pc * randomNumberGenerator.nextGaussian();
			numberOfPoreVoxels = (int)Math.ceil(myPorosity * mRCG.domainX * mRCG.domainY * mRCG.domainZ);
		}
		if (mRCG.mode.equalsIgnoreCase("range")) {			
			double lBound = mRCG.porosityBounds[0];
			double uBound = mRCG.porosityBounds[1];
			myPorosity = lBound + (uBound - lBound) * randomNumberGenerator.nextDouble();
			numberOfPoreVoxels = (int)Math.ceil(myPorosity * mRCG.domainX * mRCG.domainY * mRCG.domainZ);
		}
		if (mRCG.mode.equalsIgnoreCase("predefinedList")) {		
			fieldsToCreate = mRCG.porosityList.length;
			for (int i = 0 ; i < fieldsToCreate ; i++) {
				myPorosity = mRCG.porosityList[i];				
				numberOfPoreVoxelList[i] = (int)Math.ceil(myPorosity * mRCG.domainX * mRCG.domainY * mRCG.domainZ);
			}
		}
		
		for (int k = 0 ; k < fieldsToCreate ; k++) {
			
			//create a new 3D image
			IJ.showStatus("Creating an empty 3D image ...");
			ImagePlus randomImage = new ImagePlus();
			ByteProcessor mySeedIP = new ByteProcessor(xyz[0], xyz[1]); 
			ImageStack randomStack = new ImageStack(xyz[0], xyz[1]);
			for (int z = 0 ; z < xyz[2] ; z++) randomStack.addSlice(mySeedIP.duplicate());
			randomImage.setStack(randomStack);		
				
			//create a random subset of 3D coordinates	
			int[] randomCoords = new int[3];
			if (mRCG.mode.equalsIgnoreCase("predefinedList")) {
				myPorosity = mRCG.porosityList[k];
				numberOfPoreVoxels = numberOfPoreVoxelList[k];			
			}
			IJ.showStatus("Creating " + numberOfPoreVoxels + " random pores ... (this may take a while ...)");
			for (int i = 0 ; i < numberOfPoreVoxels ; i++) {			
				
				boolean thisVoxelWasAlreadyAPore = true;
				
				while (thisVoxelWasAlreadyAPore == true) {
					for (int j = 0 ; j < 3 ; j++) randomCoords[j] = (int)Math.floor(randomNumberGenerator.nextDouble() * xyz[j]);
							
					//check if pixel is alreaday a pore
					randomImage.setSlice(randomCoords[2] + 1);
					ImageProcessor nowIP = randomImage.getProcessor();			
					
					int nowPix = nowIP.getPixel(randomCoords[0], randomCoords[1]);
					if (nowPix == 0) {
						nowIP.putPixel(randomCoords[0], randomCoords[1], 255);
						thisVoxelWasAlreadyAPore = false;
					}
				}
			}
			
			//if a cylinder should be considered, cut it out..
			if (mRCG.shape.equalsIgnoreCase("Cylindric") == true) {
				OvalRoi oRoi = new OvalRoi(0, 0, mRCG.domainX, mRCG.domainY);
				for (int z = 1 ; z <= randomImage.getNSlices() ; z++) {
					randomImage.setPosition(z);
					ImageProcessor rIP = randomImage.getProcessor();
					rIP.setRoi(oRoi);
					rIP.setColor(Color.BLACK);
					rIP.fillOutside(oRoi);
				}
			}
			
			//save image
			randomImage.updateAndDraw();
				
			String outName = "RF" + String.format("%04d",randomSetNumber + 1) + "_Poro_" + String.format("%1.4f", myPorosity) + ".tif";
			outName = outName.replace(",", ".");
				
			jIO.tiffSaver(myPreOutPath, outName, randomImage);		
		}
	}
	
}