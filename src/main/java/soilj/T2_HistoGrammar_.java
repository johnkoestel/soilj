package SoilJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import ij.process.AutoThresholder;
import SoilJ.tools.DisplayThings;
import SoilJ.tools.HistogramStuff;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;
import SoilJ.tools.RollerCaster;
import SoilJ.tools.TailoredMaths;

import java.io.File;

import org.apache.commons.math3.stat.StatUtils;

public class T2_HistoGrammar_ extends ImagePlus implements PlugIn  {

	// loads all 3D TIFFs in the specified folder and puts them straight
	// and aligns them as such that they stand 'perfectly' straight

	public void run(String arg) {

		//construct biggish objects
		InputOutput jIO = new InputOutput();
		MenuWaiter menu = new MenuWaiter();
		InputOutput.MyFolderCollection mFC = jIO.new MyFolderCollection();
		HistogramStuff hist = new HistogramStuff();
		AutoThresholder thresh = new AutoThresholder();
		TailoredMaths maths =  new TailoredMaths();
								
		// init variables
		int i;
		MenuWaiter.HistogramMenuReturn hMR; 
		
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();
		
		//ask for threshold choice
		hMR = menu.showHistogramDialog();
		
		//read base folder and number of 3D images
		String myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
		String[] myTiffs0 = jIO.listTiffsInFolder(new File(myBaseFolder));
		
		//if not all tiffs shall be thresholded
		String[] myTiffs;
		if (hMR.filterImages) {myTiffs = jIO.filterFolderList(myTiffs0, hMR.filterTag);} else myTiffs = myTiffs0;
		
		//read gauges in the gauge folder
		mFC = jIO.getAllMyNeededFolders(myBaseFolder, myTiffs, hMR.filterTag, false);
		
		//set output path			
		String outPath = mFC.myGaugeFolder + "\\Histograms.txt";
		if (hMR.filterImages) outPath = mFC.myGaugeFolder + "\\" + hMR.filterTag + "Histograms.txt";
				
		//init histogram vectors
		int[] myBulkHistogram = new int[256 * 256];
		int[] nowHist = new int[256 * 256];
						
		//loop over 3D images
		//for (i = 0 ; i < mFC.myGaugePaths.length ; i++) {  //myTiffs.length
		for (i = 0 ; i < 1 ; i++) {  

			//try to free up some memory			
			System.gc();
			
			//load image
			nowTiff = jIO.openTiff3D(myBaseFolder + "\\" + myTiffs[i]);	
			
			//select the correct gauge and surface files			
			int[] myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(myTiffs[i], mFC.myGaugePaths, null);
			
			//get 16-bit Histogram			
			nowHist = hist.sampleThisHistogram(nowTiff, myTiffs[i], mFC, hMR, myGandS);
						
			//add to bulk histogram of all images
			for (int j = 0 ; j < nowHist.length ; j++) myBulkHistogram[j] += nowHist[j];
				
			//save individual histogram			
			jIO.saveHistogram(mFC.myGaugeFolder, nowHist, outPath, mFC.myTiffs[i].substring(0, mFC.myTiffs[i].length() - 4));

		}
		
		//save bulk histogram		
		jIO.saveHistogram(mFC.myGaugeFolder, myBulkHistogram, outPath, "BulkHistogram");		
		
		//find thresholds
		SoilJ.tools.RollerCaster cast = new SoilJ.tools.RollerCaster();
		double[] medianFilteredBulkHistogram = cast.castInt2Double(myBulkHistogram);//, 11);  //maths.oneDMedianFilter
		double[] greyValues = new double[medianFilteredBulkHistogram.length];
		for (i = 0 ; i < greyValues.length ; i++) greyValues[i] = i;
		double otsu = thresh.getThreshold(AutoThresholder.Method.Otsu, cast.castDouble2Int(medianFilteredBulkHistogram));
		double[] myThreshX = {otsu, otsu};
		double[] myThreshY = {0, 1.1 * StatUtils.max(medianFilteredBulkHistogram)};
		
		DisplayThings disp = new DisplayThings();
		disp.plotXYXY(greyValues, medianFilteredBulkHistogram, myThreshX, myThreshY, "Otsu", "grey value", "frequency");
		
		
	}
	
}