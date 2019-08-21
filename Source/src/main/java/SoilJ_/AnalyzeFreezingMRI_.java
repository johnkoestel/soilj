package SoilJ_;

/**
 *SoilJ is a collection of ImageJ plugins for the semi-automatized processing of 3-D X-ray images of soil columns
 *Copyright 2014 2015 2016 2017 John Koestel
 *
 *This program is free software: you can redistribute it and/or modify
 *it under the terms of the GNU General Public License as published by
 *the Free Software Foundation, either version 3 of the License, or
 *(at your option) any later version.
 *
 *This program is distributed in the hope that it will be useful,
 *but WITHOUT ANY WARRANTY; without even the implied warranty of
 *MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *GNU General Public License for more details.
 *
 *You should have received a copy of the GNU General Public License
 *along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import ij.IJ;
import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import SoilJ.tools.DisplayThings;
import SoilJ.tools.HistogramStuff;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MorphologyAnalyzer;
import SoilJ.tools.MorphologyAnalyzer.ProfileStatistics;
import SoilJ.tools.RollerCaster;
import SoilJ.tools.InputOutput.MyFileCollection;

/** 
 * PlotVerticalProfile is a SoilJ plugin that extracts vertical profiles of the greyscale
 * statistics within horizontal slices of 3-D images. 
 * 
 * @author John Koestel
 *
 */

public class AnalyzeFreezingMRI_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "\\";
		
		/*
		 * // set the plugins.dir property to make the plugin appear in the Plugins menu
		 * Class<?> clazz = JointThresholdDetection_.class; String url =
		 * clazz.getResource("/" + clazz.getName().replace('.', '/') +
		 * ".class").toString(); String pluginsDir = url.substring(5, url.length() -
		 * clazz.getName().length() - 6); System.setProperty("plugins.dir", pluginsDir);
		 */
		
		InputOutput jIO = new InputOutput();
		
		InputOutput.MyFileCollection mFC = jIO.new MyFileCollection();
		mFC.nowTiffPath = "Z:\\FreezingMRI\\MRI\\Quantification\\SACoHSnIsoData.tif";
		mFC.nowWidth = 256;
		mFC.nowHeight = 248; //sand
		//mFC.nowHeight = 256; //glass
		int times = 49;  //sand
		//int times = 59;  // glass beads
		int range = 417;
		
		ImagePlus[] nowTiff = new ImagePlus[times];
		
		int[] nowRange = new int[range];
		
		for (int j = 0 ; j < range ; j++) nowRange[j] = j;
		
		for (int i = 0 ; i < times ; i++) {
			
			for (int j = 0 ; j < range ; j++) nowRange[j] = nowRange[j] + range;
		
			nowTiff[i] = jIO.openTiff3DSomeSlices(mFC, nowRange);
		}
		
		//Calculate statistics...		
		Roi myRoi = new OvalRoi(33,35,187,184);  //sand
		//Roi myRoi = new OvalRoi(29,43,190,188);  //glass
		
		//output structure
		double[][] mean = new double[times][range];
		double[][] median = new double[times][range];
		double[][] min = new double[times][range];
		double[][] max = new double[times][range];
		double[][] stdev = new double[times][range];
		double[][] skew = new double[times][range];
		double[][] kurt = new double[times][range];
		
		for (int i = 0 ; i < times ; i++) {		
							
			//lowerRef
			for (int z = 0 ; z < range ; z++) {
				
				nowTiff[i].setPosition(z+1);
				ImageProcessor nowIP = nowTiff[i].getProcessor();
								
				//ImagePlus nowTest = new ImagePlus("" + z, nowIP);
				//nowTest.show();
				
				nowIP.setRoi(myRoi);
				ImageStatistics jIS = ImageStatistics.getStatistics(nowIP);
								
				mean[i][z] = jIS.mean;
				median[i][z] = jIS.median;
				min[i][z] = jIS.min;
				max[i][z] = jIS.max;
				stdev[i][z] = jIS.stdDev;
				skew[i][z] = jIS.skewness;
				kurt[i][z] = jIS.kurtosis;						
				
			}
		}
		
		//save the results
		String path = "Z:\\FreezingMRI\\MRI\\Quantification";
		int top = 33;
		int bot = 333;
		jIO.writeMRIFreezingResults("binSA_means", path, mean, times, range, top, bot);
		jIO.writeMRIFreezingResults("binSA_medians", path, median, times, range, top, bot);
		jIO.writeMRIFreezingResults("binSA_mins", path, min, times, range, top, bot);
		jIO.writeMRIFreezingResults("binSA_maxs", path, max, times, range, top, bot);
		jIO.writeMRIFreezingResults("binSA_stdev", path, stdev, times, range, top, bot);
		jIO.writeMRIFreezingResults("binSA_skewnesses", path, skew, times, range, top, bot);
		jIO.writeMRIFreezingResults("binSA_kurtoses", path, kurt, times, range, top, bot);
		
	}
	
}