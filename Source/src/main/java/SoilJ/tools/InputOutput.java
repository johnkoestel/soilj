package SoilJ.tools;

import java.awt.Color;
import java.awt.Font;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.io.IOUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import SoilJ.tools.MorphologyAnalyzer.ProfileStatistics;

/**
 *SoilJ.tools is a collection of classes for SoilJ, 
 *a collection of ImageJ plugins for the semi-automatized processing of 3-D X-ray images of soil columns
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
import ij.ImageStack;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.io.DirectoryChooser;
import ij.io.FileInfo;
import ij.io.FileSaver;
import ij.io.OpenDialog;
import ij.io.Opener;
import ij.plugin.ContrastEnhancer;
import ij.plugin.PlugIn;
import ij.plugin.StackCombiner;
import ij.plugin.filter.RankFilters;
import ij.plugin.filter.ThresholdToSelection;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

// this file collects SoilJ classes that handle IO operations

/** 
 * InputOutput is a SoilJ class that collects IO functions
 * 
 * @author John Koestel
 *
 */

public class InputOutput extends ImagePlus implements PlugIn {
	
	public void run(String arg) {
		//ok, this is not needed..
	}
	
	public class SampleTiffWrapper {
		
		public ImagePlus samTiff;
		public int[] samTiffSliceNumbers;
		public int[] samSlices;	
		public boolean hasConverged;
		
	}
	
	public class Histograms {
	
		public int histogramClasses;
		public int numberOfHistograms;
		public int[][] histograms;
		public ArrayList<String> sampleNames;
		
	}
	
	public class MyFileCollection {
		
		public String pathSep;
			
		public String myBaseFolder;
		public String myPreOutFolder;
		public String myOutFolder;
		public String mySurfaceFolder;
		public String myCutSurfaceFolder;
		public String myResultsFolder;
		public String myHistogramFolder;
		public String myInnerCircleFolder;
		public String myGradientFolder;
		public String myPoreFolder;
		public String myPOMFolder;
		public String myMatrixFolder;
		public String myMineralFolder;
		public String myFloyWarshallFolder;
		
		public String[] myTiffs;
		public String[] mySurfaceFileNames;
		public String[] myInnerCircleFiles;
		
		public String nowTiffPath;
		public FileInfo[] fileInfo;
		public String fileName;   //filename
		public String colName;		//filename without ".tif";
		public int startSlice;
		public int stopSlice;
		public int nOfSlices;
		public int nowWidth;
		public int nowHeight;
		public int bitDepth;	
		public String nowInnerCirclePath;
		
		public boolean somethingIsWrong;
		public boolean imageHasBeenLoaded;
		public String eMsg;
				
		public ImagePlus nowTiff;
		
		public double caliZ;
		public double caliX;
		public double caliY;
		
	}
	
	public MyFileCollection addCurrentFileInfo8Bit(MyFileCollection mFC) {
		
		String pathSep = "/";
		
		mFC.somethingIsWrong = false;
		
		mFC.nowTiffPath = mFC.myBaseFolder + pathSep + mFC.fileName;	
		
		FileInfo[] fI = Opener.getTiffFileInfo(mFC.nowTiffPath);
		
		mFC.fileInfo = fI;
		mFC.fileName = fI[0].fileName;
		mFC.nOfSlices = fI[0].nImages;	
		mFC.nowWidth = fI[0].width;
		mFC.nowHeight = fI[0].height;		
		mFC.bitDepth = 8 * fI[0].getBytesPerPixel();
		
		String ending = mFC.fileName.substring(mFC.fileName.length() - 2, mFC.fileName.length());
		if (ending.equalsIgnoreCase("ff")) mFC.fileName.substring(0, mFC.fileName.length() - 5);
		else mFC.colName = mFC.fileName.substring(0, mFC.fileName.length() - 4);		
					
		/*if (mFC.nOfSlices == 0) {
			mFC.somethingIsWrong = true;
			mFC.eMsg = "I have just found a file with only one slice. I cannot handle this. \n";
			mFC.eMsg += "Please select a folder that only contains 16-bit 3-D TIFF images.";
			return null;
		}
		
		if (mFC.bitDepth != 1) {
			mFC.somethingIsWrong = true;
			if (mFC.bitDepth == 4) mFC.eMsg = "I have just found a " + (mFC.bitDepth * 8) + "-bit image. I cannot handle this. \n";
			else mFC.eMsg = "I have just found an " + (mFC.bitDepth * 8) + "-bit image. I cannot handle this. \n";
			mFC.eMsg += "Please select a folder that only contains a binary 8-bit 3-D TIFF images.";
			return null;
		}*/

		return mFC;
		
	}
		
	public MyFileCollection addCurrentFileInfo(MyFileCollection mFC) {
		
		String pathSep = "/";
		
		mFC.somethingIsWrong = false;
		
		mFC.nowTiffPath = mFC.myBaseFolder + pathSep + mFC.fileName;	
		
		FileInfo[] fI = Opener.getTiffFileInfo(mFC.nowTiffPath);	
		mFC.fileInfo = fI;
		mFC.fileName = fI[0].fileName;
		mFC.nOfSlices = fI[0].nImages;	
		mFC.nowWidth = fI[0].width;
		mFC.nowHeight = fI[0].height;		
		mFC.bitDepth = fI[0].getBytesPerPixel() * 8;
		
		String ending = mFC.fileName.substring(mFC.fileName.length()-1, mFC.fileName.length());
		if (ending.equalsIgnoreCase("ff")) mFC.fileName.substring(0, mFC.fileName.length() - 5);
		else mFC.colName = mFC.fileName.substring(0, mFC.fileName.length() - 4);
		
		return mFC;
		
	}
	
	public MyFileCollection getAllMyNeededFolders(String myBaseFolder, String[] myTiffs, String filterTag, boolean useInnerCircle, boolean useSurfaceFiles) {
		
		MyFileCollection mFC = new MyFileCollection();
		//boolean hasSurfaceFiles = false; 		
		
		//assign folders to folder collection
		mFC.myBaseFolder = myBaseFolder;
		mFC.myTiffs = myTiffs;
		
		/*
		 * //read soil top surface images //not used at the moment.. String[]
		 * filesInThisFolder = listFilesInFolder(new File(myBaseFolder)); for (int i = 0
		 * ; i < filesInThisFolder.length ; i++) { if
		 * (filesInThisFolder[i].contains("SurfaceOfColumn")) { hasSurfaceFiles = true;
		 * } }
		 */
		
		if (useSurfaceFiles == true) {
			String myTopSurfaceFolder = findTheSoilSurfaceFiles(myBaseFolder);			
			mFC.mySurfaceFolder = myTopSurfaceFolder;
			String[] myTSs = listTiffsInFolder(new File(myTopSurfaceFolder));
			mFC.mySurfaceFileNames = myTSs; 
		}
		else {
			mFC.mySurfaceFolder = null;
			mFC.mySurfaceFileNames = null; 
		}
		
		//read gauges in the gauge folder
		if (useInnerCircle) {
			String myGaugeFolder = findTheInnerCircle(myBaseFolder);
			String[] myGauges = listInnerCircleFiles(myGaugeFolder, filterTag);
			mFC.myInnerCircleFolder = myGaugeFolder;
			mFC.myInnerCircleFiles = myGauges;			
		}
		else {
			mFC.myInnerCircleFolder = null;
			mFC.myInnerCircleFiles = null;	
		}
		
		return mFC;
	}
	
	public MyFileCollection addInnerCircleFolder(MyFileCollection mFC) {
		
		String myGaugeFolder = findTheInnerCircle(mFC.myBaseFolder);
		String[] myGauges = listInnerCircleFiles(myGaugeFolder, "");
		mFC.myInnerCircleFolder = myGaugeFolder;
		mFC.myInnerCircleFiles = myGauges;		
		
		return mFC;
	}
			
	public MyFileCollection addSurfaceTopographyFolder(MyFileCollection mFC) {
		
		String myTopSurfaceFolder = findTheSoilSurfaceFiles(mFC.myBaseFolder);
		String[] myTSs = listTiffsInFolder(new File(myTopSurfaceFolder));
		mFC.mySurfaceFolder = myTopSurfaceFolder;
		mFC.mySurfaceFileNames = myTSs; 
		
		return mFC;
		
	}
	
	public SampleTiffWrapper assembleRepresentativeSample(MyFileCollection mFC) {
		
		MorphologyAnalyzer morpho = new MorphologyAnalyzer();
		SampleTiffWrapper sTW = new SampleTiffWrapper();
		RankFilters rF = new RankFilters();
						
		int numOfSampleSlices = 60; // just because..				
		
		//define the default case
		int startSlice = mFC.startSlice;
		int numOfSlices = mFC.stopSlice - mFC.startSlice; 
		double increment = (double)numOfSlices / (double)numOfSampleSlices;  
				    	    
		//assign current TIFF
		int[] sampleSlices = new int[numOfSampleSlices];
		for (int j = 0 ; j < sampleSlices.length ; j++)	sampleSlices[j] = (int)Math.floor(startSlice + (double)j * increment);
		
		//Load sample image
	    ImagePlus nowTiff = new ImagePlus();	    
	    if (mFC.imageHasBeenLoaded) {
	    	ImageStack newStack = new ImageStack(mFC.nowTiff.getWidth(), mFC.nowTiff.getHeight());
	    	for (int j = 0 ; j < sampleSlices.length ; j++) {
	    		
	    		mFC.nowTiff.setPosition(sampleSlices[j]);	    		
	    		ImageProcessor nowIP = mFC.nowTiff.getProcessor();
	    		newStack.addSlice(nowIP);	    		
	    	}
	    	nowTiff.setStack(newStack);
	    }
	    else nowTiff = openTiff3DSomeSlices(mFC, sampleSlices);			    
		
		//find approximate location of column top and bottom
		ProfileStatistics pS = morpho.findProfileStatistics(nowTiff, null);
		
		//calculate standard deviation of grey values in vertical middle of column
		double[] stdMiddle = new double[numOfSampleSlices/2];
		for (int j = numOfSampleSlices / 4 ; j < numOfSampleSlices / 4 * 3 ; j++) stdMiddle[j - numOfSampleSlices / 4] = pS.std[j];
		double medianStdMiddle = StatUtils.percentile(stdMiddle, 50);
		
		//calculate standard deviation of grey values in vertical fringes of column
		double[] stdFringes = new double[numOfSampleSlices/2];
		for (int j = 0 ; j < numOfSampleSlices / 4 ; j++) stdFringes[j] = pS.std[j];
		for (int j = numOfSampleSlices / 4 * 3 ; j < numOfSampleSlices ; j++) stdFringes[j - numOfSampleSlices / 4 * 3] = pS.std[j];
		double medianStdFringes = StatUtils.percentile(stdFringes, 50);
		
		//calculate delta std
		double deltaStd = medianStdMiddle - medianStdFringes;			
		
		//label all slices below half median std as probably not part of the column...
		int firstColumnSlice = 0;
		int lastColumnSlice = 0;
		int halfNumOfSamples = numOfSampleSlices / 2;
		int j = 0;
		while (firstColumnSlice == 0){
			if (pS.std[j] > medianStdFringes + 0.5 * deltaStd) firstColumnSlice = j;
			j++;
		}
		j = halfNumOfSamples ;
		while (lastColumnSlice == 0 & j < numOfSampleSlices){
			if (pS.std[j] < medianStdFringes + 0.5 * deltaStd) lastColumnSlice = j - 1;
			j++;
		}
		if (lastColumnSlice == 0) lastColumnSlice = numOfSampleSlices - 1;
		
		//remove top and bottom fringes from samples slice set
		ImagePlus samTiff = new ImagePlus();
		ImageStack samStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		int[] samSlices = new int[lastColumnSlice + 1 - firstColumnSlice];
		for (j = firstColumnSlice ; j <= lastColumnSlice ; j++) {
			nowTiff.setPosition(j + 1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			//apply a median filter to the slices
			double smoothing = (mFC.nowWidth + mFC.nowHeight) / 400;
			rF.rank(nowIP, smoothing, RankFilters.MEDIAN); 
			
			//add it to the stack..
			samStack.addSlice(nowIP);	
			samSlices[j - firstColumnSlice] = sampleSlices[j];		
		}
		samTiff.setStack(samStack);			

		//assign output variables
		sTW.samTiff = samTiff;
		sTW.samSlices = samSlices;	
			
		return sTW;
	}
	
	//this variant is used to assemble stacks for looking for the top and bottom of a soil column..
	public SampleTiffWrapper assembleRepresentativeSample(MyFileCollection mFC, int[] fringeAndFound) { 
				
		SampleTiffWrapper sTW = new SampleTiffWrapper();

		int[] sampleSlices = new int[3];
				
		//decide whether the top or the bottom of the column is sampled..		
		if (fringeAndFound[0] <= fringeAndFound[1]){ // the top is sought
					
			//define search range
			int fringe = fringeAndFound[0];
			int mid = (fringeAndFound[1] - fringeAndFound[0]) / 2 + fringeAndFound[0];
			
			sampleSlices[0] = mid;
			sampleSlices[1] = mid + 1;
			sampleSlices[2] = mid + 2;	
			
			//decide whether search has been converging already				
			if (mid - fringe <= 3) {
				sTW.hasConverged = true;
				if (fringe == 0) {  // if convergence occurs close to slice 1, column is probably there.. let's chekc it..
					sampleSlices[0] = 0;
					sampleSlices[1] = 1;
					sampleSlices[2] = 2;
				}
			}
			else sTW.hasConverged = false;
			
	    }	
		else {  // the bottom is sought
			
			//define search range
			int fringe = fringeAndFound[0];
			int mid = (fringeAndFound[0] - fringeAndFound[1]) / 2 + fringeAndFound[1];
			
			sampleSlices[0] = mid;
			sampleSlices[1] = mid + 1;
			sampleSlices[2] = mid + 2;
			
			if (fringe - mid <= 3) sTW.hasConverged = true;
			else sTW.hasConverged = false;
		}
		
		//Load sample image
		if (mFC.imageHasBeenLoaded) {
			ImageStack newStack = new ImageStack(mFC.nowTiff.getWidth(), mFC.nowTiff.getHeight());
	    	for (int j = 0 ; j < sampleSlices.length ; j++) {
	    		
	    		mFC.nowTiff.setPosition(sampleSlices[j]);	    		
	    		ImageProcessor nowIP = mFC.nowTiff.getProcessor();
	    		newStack.addSlice(nowIP);    		
	    		
	    	}
	    	ImagePlus outTiff = new ImagePlus("", newStack);
	    	sTW.samTiff = outTiff;
		}
		else sTW.samTiff = openTiff3DSomeSlices(mFC.nowTiffPath, mFC.nowWidth, mFC.nowHeight, sampleSlices);		    
	    sTW.samSlices = sampleSlices;
		    			
		return sTW;
	}
	
	public MyFileCollection fileSelector(String heading) {
		
		//selects a file or a folder (and thus all files in it) 
		
		MyFileCollection mFC = new MyFileCollection();
		MenuWaiter mW = new MenuWaiter();
		
		String myPath = chooseAFile(heading);
		
		File myFile = new File(myPath);
		myPath = myFile.getParent();
		myPath = myPath.replace('\\', '/');
		if (myPath.endsWith("/")) myPath = myPath.substring(0, myPath.length() - 1);
		
		//open menu asking whether selected files or all files should be processed. 
		MenuWaiter.SelectFiles mFS = mW.showFileSelectionMenu(myFile);	
		
		if (!mFS.selectChosen) {			
			
			File myPathAsFile = new File(myPath);
			String[] myTiffs = listTiffsInFolder(myPathAsFile);
			
			mFC = getAllMyNeededFolders(myPath, myTiffs, "", false, false);
		}
		
		else {
			
			String[] myTiff = new String[1]; 
			myTiff[0] = myFile.getName();
					
			mFC = getAllMyNeededFolders(myPath, myTiff, "", false, false);
			
		}
		
		return mFC;
	}	
	
	public String chooseAFolder(String heading) {
		
		DirectoryChooser od = new DirectoryChooser(heading);
		String dir = od.getDirectory();  
		if (null == dir) return null; // dialog was canceled
		
		return dir;
	}
	
	public String chooseAFile(String heading) {		
		
		OpenDialog od = new OpenDialog(heading);		
		String path = od.getPath();
		if (null == path) return null; // dialog was canceled
		
		return path;
	}
	
	public String chooseAFile(String heading, String dir, String file) {		
		
		OpenDialog od = new OpenDialog(heading, dir, file);		
		String path = od.getPath();
		if (null == path) return null; // dialog was canceled
		
		return path;
	}
	
	public String getTheFolderAbove(String myDir) {
		
		String pathSep = "/";
		
		int i;
		int cutoff = 0;
		
		char[] myDirAsChar = myDir.toCharArray(); 
		
		for (i = myDirAsChar.length-3 ; i > 0 ; i--) {
						
			if (myDirAsChar[i] == pathSep.charAt(0)) { 
				cutoff = i;
				break;
			}
			
		}
		
		return myDir.substring(0, cutoff);
		
	}
	
	public String[] listFoldersInFolder(final File folder) {

		String testFile;
		ArrayList<String> myFiles = new ArrayList<String>();
				
		for (final File fileEntry : folder.listFiles()) {
			testFile = fileEntry.getName();
			if (fileEntry.isDirectory()) {
				myFiles.add(testFile);				
			}
		}
		
		if (myFiles.size() == 0) return null;
		
		String[] myFolder = new String[myFiles.size()];
		for (int i = 0 ; i < myFiles.size() ; i++) {
			myFolder[i] = myFiles.get(i);
		}
		
		return myFolder;
	}
	
	public int[] getTheCorrectGaugeNSurfaceFiles(String myTiff, String[] myGauges, String[] myTSs) {
		
		String myName = myTiff.substring(0,myTiff.length() - 4);
		int[] myGandS = new int[2];	
				
		for (int i = 0 ; i < myGauges.length ; i++) {
			
			boolean foundName = false;
			for (int j = 0 ; j < myGauges[i].length() - myName.length() - 1 ; j++) {
				String myNextChar = myGauges[i].substring(j+myName.length(), j+myName.length()+1);
				String region2Compare = myGauges[i].substring(j, j + myName.length());
				if (myName.equalsIgnoreCase(region2Compare) & myNextChar.equalsIgnoreCase(".")) foundName = true;
				if (foundName) {
					myGandS[0] = i;
					break;
				}
			}
		}
				
		if (myTSs != null) {
			for (int i = 0 ; i < myTSs.length ; i++) {
				
				Boolean foundName = false;
				for (int j = 0 ; j < myTSs[i].length() ; j++) {
					if (myTSs[i].regionMatches(j, myName, 0, myName.length())) foundName = true;
					if (foundName) myGandS[1] = i;
				}
			}	
		}
			
		return myGandS;
	}
	
	public int[] getTheCorrectGaugeNSurfaceFiles(MyFileCollection mFC) {
		
		String myName = mFC.colName;	
		
		//remove "_Thickness" from file name 
		int lomN = myName.length(); 
		if (lomN > 10) {
			String endOfMyName = myName.substring(lomN - 10, lomN);
			if (endOfMyName.equalsIgnoreCase("_Thickness")) {
				myName = myName.substring(0, lomN - 10);
			}
		}
		
		int[] myGandS = new int[2];	
				
		for (int i = 0 ; i < mFC.myInnerCircleFiles.length ; i++) {
						
			String nowName = mFC.myInnerCircleFiles[i];
			int endOfName = -1;
			for (int j = nowName.length() - 1 ; j > 0 ; j--) {
				endOfName = j;
				if (nowName.substring(j-1,j).equalsIgnoreCase(".")) break;
			}
			nowName = nowName.substring(endOfName - myName.length() - 1, endOfName-1);
			
			//IJ.error(nowName + " --- " + myName);
			
			if (nowName.length() == myName.length()) {					
				if (myName.equalsIgnoreCase(nowName)) {
					myGandS[0] = i;
					break;					
				}
			}	
		}
				
		if (mFC.mySurfaceFolder != null) {
			for (int i = 0 ; i < mFC.mySurfaceFileNames.length ; i++) {
				
				Boolean foundName = false;
				int endIndex =  mFC.mySurfaceFileNames[i].length();
				for (int j = 0 ;  j < mFC.mySurfaceFileNames[i].length() ; j++) {
					String nowName = mFC.mySurfaceFileNames[i];
					if (nowName.substring(endIndex - j - 2, endIndex - j - 1).equalsIgnoreCase(".")) {
						break;
					}
					endIndex--;
				}
				String surName = mFC.mySurfaceFileNames[i].substring(0, endIndex - 3);				
				if (surName.equalsIgnoreCase(myName)) {
					foundName = true;
					myGandS[1] = i;	
				}		
			}	
		}
			
		return myGandS;
	}

	public boolean testIfFolderIsPresent(File folder, String folderName) {
				
		String testFile;
		boolean yes = false;
		
		if (folder.listFiles().length > 0) {
			for (final File fileEntry : folder.listFiles()) {
				testFile = fileEntry.getName();
				if (testFile.equalsIgnoreCase(folderName)) {
					yes = true;
					break;
				}
			}
		}
		
		return yes;
	}

	public String[] listTiffsInFolder(File folder) {

		String testFile, subString;
	
		ArrayList<String> myTiffs = new ArrayList<String>();
		
		if (folder.listFiles().length > 0) {
			for (final File fileEntry : folder.listFiles()) {
				testFile = fileEntry.getName();
				if (testFile.length() >= 3) {
					subString = testFile.substring(testFile.length() - 3,testFile.length());
					if (subString.equals("tif")) myTiffs.add(testFile);									
				}			
			}
			
			Collections.sort(myTiffs);
			
			String[] myFiles = new String[myTiffs.size()];
						
			for (int i = 0 ; i < myFiles.length ; i++) myFiles[i] = myTiffs.get(i); 
			
			return myFiles;
		}
		
		else return null;				
		
	}
	
	public String[] listFilesInFolder(File folder) {

		String testFile;
		
		int cc = 0;

		if (folder.listFiles().length > 0) {
			for (final File fileEntry : folder.listFiles()) {
				testFile = fileEntry.getName();
				if (testFile.length() >= 3) {
					cc++;									
				}			
			}
			
			String[] myFiles = new String[cc];
			cc = 0;
			for (final File fileEntry : folder.listFiles()) {
				testFile = fileEntry.getName();
				if (testFile.length() >= 3) {					
					myFiles[cc] = testFile;
					cc++;						
				}			
			}
			
			return myFiles;
		}
		
		else return null;				
		
	}
	
	public String findTheInnerCircle(String myBaseFolder){
		
		String pathSep = "/";
		
		File baseFolder = new File(myBaseFolder); 		//check if InnerCircle is in the base folder
		String thisParentFolder = null;
		String InnerCircle = "InnerCircle";				
		String myGaugeFolder = null;
		String testFile;
		boolean path2GaugeNeedsAChangeAtBase = false;
		boolean path2GaugeNeedsAChangeAtParent = false;
		
		for (final File fileEntry : baseFolder.listFiles()) {
			
			testFile = fileEntry.getName();			
			
			if (testFile.equals(InnerCircle)) {
				myGaugeFolder = myBaseFolder + pathSep + InnerCircle;			
			}		
			
			if (testFile.equals("Path2Gauge.txt")) {
				
				MyMemory myMem = readBrainFile(myBaseFolder + pathSep + "Path2Gauge.txt");
				
				//check if path is pointing at a Samba server
				if (myMem.gaugeFolder.substring(0,2).equalsIgnoreCase("smb")) {
					IJ.error("SoilJ does not support Samba yet. I am sorry.\n Please copy your files to your local Linux system.");
					return null;
				}
							
				//test if this really is the InnerCircle folder
				if (myMem.gaugeFolder == null | myMem.gaugeFolder.length() < 13) {
					IJ.error("This plugin requires InnerCircle files! \n\n "
							+ "If you have not yet created such files, please abort this program and first \n"
							+ "the FindColumnOutlines plugin.\n\n"
							+ "Otherwise please tell me the path to your InnerCircle files. Thank you.");
					myGaugeFolder=null;
				}
				
				else {
					int end = myMem.gaugeFolder.length();	
					String myCompareString = myMem.gaugeFolder.substring(end-11, end);
					
					if (myCompareString.equalsIgnoreCase("InnerCircle")) myGaugeFolder = myMem.gaugeFolder;
					
					//catch Linux different path syntax
					if (myCompareString.substring(myCompareString.length()-1).equalsIgnoreCase("/")) { 
						myCompareString = myMem.gaugeFolder.substring(end-12, end-1);
						if (myCompareString.equalsIgnoreCase("InnerCircle")) myGaugeFolder = myMem.gaugeFolder.substring(0,end-1);
					}
					
					else path2GaugeNeedsAChangeAtBase = true;
				}				
			}	
		}

		if (myGaugeFolder==null) {						//check if InnerCircle is in the folder above the base folder
			
			thisParentFolder = getTheFolderAbove(myBaseFolder);
			File myParentFolder = new File(thisParentFolder);
			
			for (final File fileEntry : myParentFolder.listFiles()) {
				
				testFile = fileEntry.getName();			
				
				if (testFile.equals(InnerCircle)) {
					myGaugeFolder = thisParentFolder + pathSep + InnerCircle;			
				}			
				
				if (testFile.equals("Path2Gauge.txt")) {
					
					MyMemory myMem = readBrainFile(thisParentFolder + pathSep + "Path2Gauge.txt");
					
					//check if path is pointing at a Samba server
					if (myMem.gaugeFolder.substring(0,2).equalsIgnoreCase("smb")) {
						IJ.error("SoilJ does not support Samba yet. I am sorry.\n Please copy your files to your local Linux system.");
						return null;
					}
					
					//test if this really is the InnerCircle folder
					int end = myMem.gaugeFolder.length();
					
					String myCompareString = myMem.gaugeFolder.substring(end-11, end);
					if (myCompareString.equalsIgnoreCase("InnerCircle")) myGaugeFolder = myMem.gaugeFolder;
					
					//catch Linux different path syntax
					if (myCompareString.substring(myCompareString.length()-1).equalsIgnoreCase("/")) { 
						myCompareString = myMem.gaugeFolder.substring(end-12, end-1);
						if (myCompareString.equalsIgnoreCase("InnerCircle")) myGaugeFolder = myMem.gaugeFolder.substring(0,end-1);
					}
					
					else path2GaugeNeedsAChangeAtParent = true;
				}	
			}			
		}	
		
		while (myGaugeFolder==null) {						//If still not found ask for location..
			
			myGaugeFolder = chooseAFolder("Could not find the location of the 'Inner Circle'.. please show me where it is! I need it!");
			
			//check if path is pointing at a Samba server
			if (myGaugeFolder.substring(0,2).equalsIgnoreCase("smb")) {
				IJ.error("SoilJ does not support Samba yet. I am sorry.\n Please copy your files to your local Linux system.");
				return null;
			}
			
			//test if this really is the InnerCircle folder
			int end = myGaugeFolder.length();
			
			if (end < 13) myGaugeFolder = null;
			else{
				
				String myCompareString = myGaugeFolder.substring(end-11, end);
				if (!myCompareString.equalsIgnoreCase("InnerCircle")) myGaugeFolder = null;
				else myGaugeFolder = myGaugeFolder;
				
				if (myCompareString.substring(myCompareString.length()-1).equalsIgnoreCase("/")) {
					myCompareString = myGaugeFolder.substring(end-12, end-1);
					if (!myCompareString.equalsIgnoreCase("InnerCircle")) myGaugeFolder = null;
					else myGaugeFolder = myGaugeFolder.substring(0, end-1);				
				}
			}
		}
		
		writeStringIntoAsciiFile(myBaseFolder + pathSep + "Path2Gauge.txt", myGaugeFolder);
		
		return myGaugeFolder;
		
	}
	
	public String findTheSoilSurfaceFiles(String myBaseFolder){
		
		String pathSep = "/";
		
		File baseFolder = new File(myBaseFolder); 		//check if SoC is in the base folder
		String SoC = "SurfaceOfColumn";				
		String mySurfaceFolder = null;
		String testFile;
		
		for (final File fileEntry : baseFolder.listFiles()) {
			
			testFile = fileEntry.getName();			
			
			if (testFile.equals(SoC)) {
				mySurfaceFolder = myBaseFolder + pathSep + SoC;			
			}		
		}
		
		int tester = 0;
		while (mySurfaceFolder==null) {						//If not found ask for location..
			
			if (tester == 0) mySurfaceFolder = chooseAFolder("Could not find the location of the 'SurfaceOfColumn'.. please show me where it is! I need it!");
			else mySurfaceFolder = chooseAFolder("At this locations were no surface files.. please show me where they really are or press cancel. Thank you.");
			
			tester++;
			File newFolder = new File(mySurfaceFolder);
			
			String newName = newFolder.getName();
			
			if (!newName.equals(SoC)) mySurfaceFolder = null;
						
			if (tester > 3) return null;
		}
		
		return mySurfaceFolder;
		
	}
	
	public String[] listInnerCircleFiles(String myGaugeFolder, String filterTag) {

		String pathSep = "/";
		
		File folder = new File(myGaugeFolder);
		String testFile, subString;
		ArrayList<String> myF = new ArrayList<String>();  
		
		for (final File fileEntry : folder.listFiles()) {
	
			testFile = fileEntry.getName();
			subString = testFile.substring(testFile.length() - 3,
					testFile.length());
			if (subString.equals("txt") && testFile.contains("Gauge")) {
				myF.add(testFile);			
			}			
			if (subString.equals("txt") && testFile.contains("Steel")) {
				myF.add(testFile);			
			}	
		}
		
		String[] myFiles = new String[myF.size()];		
		for (int i = 0 ; i < myF.size() ; i++) {			
			myFiles[i] = myGaugeFolder + pathSep + myF.get(i);			
		}
				
		return myFiles;
	}
	
	public String[] listSurfaceFiles(InputOutput.MyFileCollection mFC, String filterTag) {

		String pathSep = "/";
		
		File folder = new File(mFC.mySurfaceFolder);
		String testFile, subString;				
		String[] myFiles = listTiffsInFolder(folder);
		
		int cc = 0;

		for (final File fileEntry : folder.listFiles()) {
			testFile = fileEntry.getName();
			subString = testFile.substring(testFile.length() - 4,
					testFile.length());
			if (subString.contains("tif") & mFC.nOfSlices == 2) {				 
				cc++;			
			}			
		}
		
		String[] myFiles0 = new String[cc];
		cc = 0;
		for (final File fileEntry : folder.listFiles()) {
			testFile = fileEntry.getName();
			subString = testFile.substring(testFile.length() - 3,
					testFile.length());
			if (subString.contains("tif") & mFC.nOfSlices == 2) {
				myFiles0[cc] = mFC.mySurfaceFolder + pathSep + testFile;
				cc++;
			}
		}
		
		//filter gauge files in case it was desired
		if (filterTag!="") {myFiles = filterFolderList(myFiles0, filterTag);} else myFiles = myFiles0; 
				
		return myFiles;
	}
	
	public String[] listSteelGaugeFiles(String myGaugeFolder) {
		
		String pathSep = "/";

		File folder = new File(myGaugeFolder);
		String testFile, subString;	
		
		int cc = 0;

		for (final File fileEntry : folder.listFiles()) {
			testFile = fileEntry.getName();
			subString = testFile.substring(testFile.length() - 3,
					testFile.length());
			if (subString.equals("txt") && testFile.contains("Gauge")) {
				cc++;
				
			}			
		}
		
		String[] myFiles = new String[cc];
		cc = 0;
		for (final File fileEntry : folder.listFiles()) {
			testFile = fileEntry.getName();
			subString = testFile.substring(testFile.length() - 3,
					testFile.length());
			if (subString.equals("txt") && testFile.contains("Gauge")) {
				myFiles[cc] = myGaugeFolder + pathSep + testFile;
				cc++;
			}
		}
				
		return myFiles;
	}
	
	public String[] filterFolderList(String[] myTiffs0, String filterTag) {
		
		int cc = 0;
		int tiffListLength = myTiffs0.length;
		for (int i = 0 ; i < tiffListLength ; i++) {
			if (myTiffs0[i].contains(filterTag)) {
				cc++;
			}
		}
		String[] myTiffs = new String[cc];
		cc = 0;
		for (int i = 0 ; i < myTiffs0.length ; i++) if (myTiffs0[i].contains(filterTag)) {
			myTiffs[cc] = myTiffs0[i];
			cc++;
		}
		
		return myTiffs;		
	}
	
	public ImagePlus openStack(String myDir) {
		
		String pathSep = "/";

		int i;		
		int lengthOfTiffStack;
		ImagePlus newImg;
		ImagePlus imgStack = new ImagePlus();
		ImageProcessor ip;

		//check out how many entries there are in tiff stack		
		String[] myFiles;
		File myFolder = new File(myDir);
		myFiles = listTiffsInFolder(myFolder);		
		Object array = myFiles;
		lengthOfTiffStack = Array.getLength(array);
		
		// open first image to get the file info
		newImg = IJ.openImage(myDir + pathSep + myFiles[0]);

		// also create an image processor of it..
		ip = newImg.getProcessor();

		// create an image stack..
		ImageStack myStack = new ImageStack(newImg.getFileInfo().width,
				newImg.getFileInfo().height);
				
		// create the stack 
		for (i = 0; i < lengthOfTiffStack; i++) {
			// repeat the above for all images..
			IJ.showStatus("Reading: " + (i + 1) + "/" + lengthOfTiffStack);
			newImg = IJ.openImage(myDir + pathSep + myFiles[i]);
			ip = newImg.getProcessor();
			myStack.addSlice(myFiles[i], ip);

		}

		imgStack.setStack(myStack);

		return imgStack;
	}

	public ImagePlus openTiff3D(String nowTiffPath) {
		
		Opener oT3D = new Opener();
		ImagePlus nowTiff;
		
		nowTiff = oT3D.openImage(nowTiffPath);
		
		return nowTiff;
		
	}
	
	public ImagePlus openTiff3D(MyFileCollection mFC) {
		
		Opener oT3D = new Opener();
		ImagePlus nowTiff;
		
		nowTiff = oT3D.openImage(mFC.nowTiffPath);
		
		return nowTiff;
		
	}
	
	public ImagePlus openTiff3DSomeSlices(MyFileCollection mFC, int[] sampleSlices) {
		
		Opener oT3D = new Opener();
		ImagePlus nowTiff = new ImagePlus();		
				
		ImageStack outStack = new ImageStack(mFC.nowWidth, mFC.nowHeight);
		int numberOfSlides = sampleSlices.length;
			
		for (int i = 0 ; i < numberOfSlides ; i++) {
			
			IJ.showStatus("Opening sample slice " + (i + 1) + "/" + sampleSlices.length + " ...");
			
			try {
				ImageProcessor nowIP = oT3D.openImage(mFC.nowTiffPath, sampleSlices[i]).getProcessor();
				outStack.addSlice(nowIP);
			}	
			catch(Exception e){
				IJ.error("Something was fishy when trying to open "+ sampleSlices.length +" sample slices from " + mFC.colName);
			}		
			
		}
		
		nowTiff.setStack(outStack);
		
		return nowTiff;
		
	}
	
	public ImagePlus openTiff3DSomeSlices(String nowTiffPath, int width, int height, int[] sampleSlices) {
		
		Opener oT3D = new Opener();
		ImagePlus nowTiff = new ImagePlus();		
				
		ImageStack outStack = new ImageStack(width, height);
			
		for (int i = 0 ; i < sampleSlices.length ; i++) {
			
			IJ.showStatus("Opening sample slice " + (i + 1) + "/" + sampleSlices.length + " ...");
			
			try {
				ImageProcessor nowIP = oT3D.openImage(nowTiffPath, sampleSlices[i]).getProcessor();
				outStack.addSlice(nowIP);
			}	
			catch(Exception e){}		
			
		}
		
		nowTiff.setStack(outStack);
		
		return nowTiff;
		
	}
	
	public ImagePlus openTiff2D(String nowTiffPath) {
		
		Opener oT2D = new Opener();
		ImagePlus nowTiff;
		
		nowTiff = oT2D.openImage(nowTiffPath);
		
		return nowTiff;
		
	}
	
	public void tiffSaver(MyFileCollection mFC, ImagePlus outImg) {
		
		String pathSep = "/";
	
		FileSaver fs = new FileSaver(outImg);		
		String outPath = mFC.myOutFolder + pathSep + mFC.fileName;
		fs.saveAsTiffStack(outPath);
		
	}
	
	public void tiffSaver(String myDir, String outName, ImagePlus outImg) {
		
		String pathSep = "/";
		
		FileSaver fs = new FileSaver(outImg);		
		String outPath = myDir + pathSep + outName;
		fs.saveAsTiffStack(outPath);
				
	}
	
	public void save2DHistogramAsTiff(MyFileCollection mFC, double[][] histogram2D) {
		
		ImageProcessor histAsImage = new FloatProcessor(256, 256);		
				
		//create Image
		IJ.showStatus("Compiling histogram topography ...");
		for (int y = 0 ; y < 256 ; y++) {			
			for (int x = 0 ; x < 256 ; x++) {			
				double pix2Put = Math.log10(histogram2D[x][y]);
				if (pix2Put < -10) pix2Put = -10;
				pix2Put += 10;
				histAsImage.putPixelValue(x, y, pix2Put);								
			}			
		}
		
		//scale gray value between 0 and 1;
		histAsImage.multiply(1/histAsImage.getMax());
		
		//scale to one tenth of real gray value	
		IJ.showStatus("Scaling image to a nice size ...");
		histAsImage.setInterpolationMethod(ImageProcessor.BILINEAR);
		ImageProcessor scaledHistogram = histAsImage.resize(6554, 6554, true);
		scaledHistogram.setMinAndMax(0.25, 0.75);
		
		//save image
		ImagePlus outTiff = new ImagePlus();
		IJ.showStatus("Saving image ...");
		outTiff.setProcessor(scaledHistogram);
		outTiff.updateAndDraw();
		
		save2DTiff(mFC, outTiff);	
		
	}
	
	public void saveAsTiffStack4GeoDict(String myDir, String outName, ImagePlus outTiff) {
		
		String pathSep = "/";

		String subDir = outName.substring(0,outName.length()-4);
		String myOutPath = myDir + pathSep + subDir;
		new File(myOutPath).mkdir();
		String runCommand = "format=TIFF save=" + myOutPath + pathSep + subDir + "0000.tif";
		IJ.run(outTiff, "Image Sequence... ", runCommand);
		
	}
	
	public void save2DTiff(String myDir, String outName, ImagePlus outImg) {
		
		String pathSep = "/";
		
		FileSaver fs = new FileSaver(outImg);		
		String outPath = myDir + pathSep + outName;
		fs.saveAsTiff(outPath);
		
	}
	
	public void save2DTiff(MyFileCollection mFC, ImagePlus outImg) {
		
		FileSaver fs = new FileSaver(outImg);		
		String outPath = mFC.nowTiffPath;
		fs.saveAsTiff(outPath);
		
	}
	
	public int[][] startAndStopReader(String myFileLocation, int numberOfFiles) throws IOException {
		
		FileInputStream inputStream = new FileInputStream(myFileLocation);
		String startsAndStops;
		int[][] sns = new int[numberOfFiles][2]; 
		
		//read file
		try {
			startsAndStops = IOUtils.toString(inputStream);
		} 
		finally {		       
			inputStream.close();
		}
		
		//extract information from file
		int cc = 0;
		String stringBuffer = "";
		for (int i = 0 ; i < startsAndStops.length(); i++) {
			String nowString = startsAndStops.substring(i,i+1);			
			if (nowString.charAt(0) != '\t' && nowString.charAt(0) != '\r' && nowString.charAt(0) != '\n') {
				if (stringBuffer.isEmpty()) stringBuffer = nowString;
				else stringBuffer = stringBuffer + nowString;
			}
			if (nowString.charAt(0) == '\t') {
				sns[cc][0] = Integer.parseInt(stringBuffer);
				stringBuffer = "";
			}
			if (nowString.charAt(0) == '\n') {
				sns[cc][1] = Integer.parseInt(stringBuffer);
				stringBuffer = "";
				cc++;
			}			
		}
		return sns;
		
	}
	
	public boolean writeROIMorphoResults(String sampleName, String path, MorphologyAnalyzer.ROIMorphoProps rMP) {
       
		try{
            
			//open file
			FileOutputStream fos = new FileOutputStream(path);
            Writer w = new BufferedWriter(new OutputStreamWriter(fos));
            
            //write pre-header
            String myString = 	"sample name:                     " + sampleName + "\n";
           
            myString += "ROI bulk volume:                 "	+ String.format("%1.6e",rMP.roiBulkVolume) + " vx\n";
            myString += "phase volume:                    "	+ String.format("%1.6e",rMP.phaseVolume) + " vx\n"; 
            myString += "phase volume fraction:           "	+ String.format("%1.4e",rMP.phaseVolumeFraction) + "\n";
            myString += "surface area:                    "	+ String.format("%1.6e",rMP.surfaceArea) + " px\n";
           
            myString += "average phase thickness:         "	+ String.format("%3.3f",rMP.averagePhaseDiameter) + " res\n";
            myString += "average distance phase-boundary: "	+ String.format("%3.3f",rMP.averageDistance2PhaseBoundary) + " res\n";
            myString += "surface fractal dimension:       "	+ String.format("%1.3f",rMP.surfaceFractalDimension) + "\n";
            
            myString += "vertical/horizontal anisotropy:  "	+ String.format("%1.3f",rMP.vertHorzAnisotropy) + "\n";
            myString += "direction of anisotropy :        " + rMP.aniDirection + "\n";
            myString += "horizontal anisotropy:           "	+ String.format("%1.3f",rMP.horizontalAnisotropy) + "\n";
            
            myString += "gamma (connection probability):  "	+ String.format("%1.3f",rMP.gamma) + "\n";
            myString += "Euler number:                    "	+ String.format("%1.6e",rMP.eulerNumber) + "\n";
            myString += "phase is percolating:            "	+ rMP.phasePercolates + "\n";
            myString += "critical diameter:               "	+ String.format("%2.2f",rMP.criticalPhaseDiameter) + "\n";

            myString += "percolating volume fraction      " + String.format("%1.4e",rMP.percolatingVolumeFraction) + "\n";
            myString += "volume fraction connected to top " + String.format("%1.4e",rMP.volumeFractionConnected2Top) + "\n";
            myString += "KT87 volume fraction:            "	+ String.format("%1.4e",rMP.KT87VolumeFraction) + "\n";
            
            String cMyString = myString.replace(',', '.');
            w.write(cMyString);
            w.flush();
            
            //close file
            w.close();
            
        }catch(Exception e){
        
        	return true;
        	
        }
        
		return false;
        
    }
	
	public boolean writeNetworkAnalysesResults(MyFileCollection mFC, MorphologyAnalyzer.NetworkProps mNP) {
	       
		String path = mFC.myOutFolder + "/" + mFC.colName + ".net"; 
		
		try{
            
			//open file
			FileOutputStream fos = new FileOutputStream(path);
            Writer w = new BufferedWriter(new OutputStreamWriter(fos));
           
            //write pre-header
            String myString = 	"sample name:                               " + mFC.colName + "\n";
           
            myString += "\nglobal properties\n\n";
            
            myString += "global number of edges:                    " + mNP.globalEdges + "\n";
            myString += "global number of slabs:                    " + mNP.globalSlabs + "\n"; 
            myString += "global number of vertices:                 " + mNP.globalVertices + "\n";
            myString += "global number of dead-ends:                " + mNP.globalDeadEnds + "\n";
           
            myString += "number of single vertex graphs             " + mNP.numberOfSingleVertexGraphs + "\n";
            myString +=	"number of percolating top vertices         " + mNP.numberOfPercolatingTopPoints + "\n";
            myString +=	"number of pecolating bottom vertices       " + mNP.numberOfPercolatingBotPoints + "\n";
            myString +=	"number of vertices connected to top        " + mNP.numberOfCon2TopPoints + "\n";
            myString +=	"number of vertices connected to bottom     " + mNP.numberOfCon2BotPoints + "\n";

            myString += "global edge length                         " + String.format("%1.3f",mNP.globalEdgeLength) + "\n";
            myString += "global straight edge length                " + String.format("%1.3f",mNP.globalEuclidEdgeLength) + "\n";
            myString += "global tortuosity                          " + String.format("%1.3f",mNP.globalTortuosity) + "\n";
            myString += "mean local tortuosity                      " + String.format("%1.3f",mNP.globalMeanTortuosity) + "\n";
            myString += "mean coordination number                   " + String.format("%1.3f",mNP.globalMeanCN) + "\n";
            myString += "mean declination from Z axis               " + String.format("%1.3f",mNP.globalMeanAngleFromZ) + "\n";
            
            myString += "\nproperties of percolating graphs\n\n";

            myString += "number of percolating edges:               " + mNP.percolatingEdges + "\n";
            myString += "number of percolating slabs:               " + mNP.percolatingSlabs + "\n"; 
            myString += "number of percolating vertices:            " + mNP.percolatingVertices + "\n";
            myString += "number of percolating dead-ends:           " + mNP.percolatingDeadEnds + "\n";
            myString += "percolating mean local tortuosity          " + String.format("%1.3f",mNP.percolatingTortuosity) + "\n";
            myString += "percolating mean coordination number       " + String.format("%1.3f",mNP.percolatingMeanCN) + "\n";
            myString += "percolating mean declination from Z axis   " + String.format("%1.3f",mNP.percolatingMeanAngleFromZ) + "\n";

            myString += "\nproperties of graphs touching the top (but not the bottom)\n\n";
            
            myString += "number of top-connecting edges:            " + mNP.con2TopEdges + "\n";
            myString += "number of top-connecting  slabs:           " + mNP.con2TopSlabs + "\n"; 
            myString += "number of top-connecting  vertices:        " + mNP.con2TopVertices + "\n";
            myString += "number of top-connecting  dead-ends:       " + mNP.con2TopDeadEnds + "\n";
            myString += "top-connecting  mean local tortuosity      " + String.format("%1.3f",mNP.con2TopTortuosity) + "\n";
            myString += "top-connecting  mean coordination number   " + String.format("%1.3f",mNP.con2TopMeanCN) + "\n";
            myString += "top-connecting  mean declination from Z    " + String.format("%1.3f",mNP.con2TopMeanAngleFromZ) + "\n";
            
            myString += "\nproperties of graphs touching the bottom (but not the top)\n\n";
            
            myString += "number of bottom-connecting edges:         " + mNP.con2BotEdges + "\n";
            myString += "number of bottom-connecting slabs:         " + mNP.con2BotSlabs + "\n"; 
            myString += "number of bottom-connecting vertices:      " + mNP.con2BotVertices + "\n";
            myString += "number of bottom-connecting dead-ends:     " + mNP.con2BotDeadEnds + "\n";
            myString += "bottom-connecting mean local tortuosity    " + String.format("%1.3f",mNP.con2BotTortuosity) + "\n";
            myString += "bottom-connecting mean coordination number " + String.format("%1.3f",mNP.con2BotMeanCN) + "\n";
            myString += "bottom-connecting  mean declination from Z " + String.format("%1.3f",mNP.con2BotMeanAngleFromZ) + "\n";
            
            myString += "\nproperties of fastest percolating path\n\n";

            myString += "number of fastest-path slabs:              " + mNP.fastestPathSlabs + "\n"; 
            myString += "number of fastest-path vertices:           " + mNP.fastestPathVertices + "\n";
            myString += "fastest-path resistance:                   " + String.format("%1.3f",mNP.fastestPathResistance) + "\n";
            myString += "fastest-path bottleneck:                   " + String.format("%1.3f",mNP.fastestPathBottleneck) + "\n";
            myString += "fastest-path edge tortuosity:              " + String.format("%1.3f",mNP.fastestPathEdgeTortuosity) + "\n";
            myString += "fastest-path global tortuosity:            " + String.format("%1.3f",mNP.fastestPathGlobalTortuosity) + "\n";
            myString += "mean fastest-path coordination number:     " + String.format("%1.3f",mNP.meanFastestPathCN) + "\n";
    	
            String cMyString = myString.replace(',', '.');
            w.write(cMyString);
            w.flush();
            
            //close file
            w.close();
            
        }catch(Exception e){
        
        	return true;
        	
        }
        
		return false;
        
    }
	
	public boolean writeMLBHallePA3DResults(String path, MorphologyAnalyzer.MLJLabellingResults labelResults) {
	       
		try{
            
			//open file
			FileOutputStream fos = new FileOutputStream(path);
            Writer w = new BufferedWriter(new OutputStreamWriter(fos));
            
            //write 
            String myString = 	"volume\tsurface\teuler\n";           
            
            for (int i = 0 ; i < labelResults.volume.length ; i++) {
            	myString += String.format("%1.12e",labelResults.volume[i]) + "\t";            
            	myString += String.format("%1.12e",labelResults.surface[i]) + "\t"; 
            	myString += String.format("%1.12e",labelResults.euler[i]) + "\n";
            }
                      
            String cMyString = myString.replace(',', '.');
            w.write(cMyString);
            w.flush();
            
            //close file
            w.close();
            
        } catch(Exception e) {
        
        	return true;
        	
        }
        
		return false;
        
    }
	
	public boolean writeROIAnisoResults(String sampleName, String path, MorphologyAnalyzer.AnisotropyResults arreArre) {
		
		try{
            
			//open file
			FileOutputStream fos = new FileOutputStream(path);
            Writer w = new BufferedWriter(new OutputStreamWriter(fos));
            
            //write header
        	String myHeader = "direction\t" + "sum\t" + "mean\t" + "stdev\t" + "CV\n";
        	w.write(myHeader);
                        
            //write the data
            String outString = "z\t";
            for (int i = 0 ; i < 3 ; i++) outString += String.format("%1.6e", arreArre.z[i]) + "\t";
            outString += String.format("%1.6e", arreArre.z[3]) + "\n";
     
            outString += "x\t";
            for (int i = 0 ; i < 3 ; i++) outString += String.format("%1.6e", arreArre.x[i]) + "\t";
            outString += String.format("%1.6e", arreArre.x[3]) + "\n";
            
            outString += "y\t";
            for (int i = 0 ; i < 3 ; i++) outString += String.format("%1.6e", arreArre.y[i]) + "\t";
            outString += String.format("%1.6e", arreArre.y[3]) + "\n";
            
            outString += "xy\t";
            for (int i = 0 ; i < 3 ; i++) outString += String.format("%1.6e", arreArre.xy[i]) + "\t";
            outString += String.format("%1.6e", arreArre.xy[3]) + "\n";
            
            outString += "yx\t";
            for (int i = 0 ; i < 3 ; i++) outString += String.format("%1.6e", arreArre.yx[i]) + "\t";
            outString += String.format("%1.6e", arreArre.yx[3]) + "\n";
            
            outString += "dxz\t";
            for (int i = 0 ; i < 3 ; i++) outString += String.format("%1.6e", arreArre.dxz[i]) + "\t";
            outString += String.format("%1.6e", arreArre.dxz[3]) + "\n";
            
            outString += "dyz\t";
            for (int i = 0 ; i < 3 ; i++) outString += String.format("%1.6e", arreArre.dyz[i]) + "\t";
            outString += String.format("%1.6e", arreArre.dyz[3]) + "\n";
            
            outString += "dxyz\t";
            for (int i = 0 ; i < 3 ; i++) outString += String.format("%1.6e", arreArre.dxyz[i]) + "\t";
            outString += String.format("%1.6e", arreArre.dxyz[3]) + "\n";
            
            outString += "dyxz\t";
            for (int i = 0 ; i < 3 ; i++) outString += String.format("%1.6e", arreArre.dyxz[i]) + "\t";
            outString += String.format("%1.6e", arreArre.dyxz[3]) + "\n";
            
            outString += "uxz\t";
            for (int i = 0 ; i < 3 ; i++) outString += String.format("%1.6e", arreArre.uxz[i]) + "\t";
            outString += String.format("%1.6e", arreArre.uxz[3]) + "\n";
            
            outString += "uyz\t";
            for (int i = 0 ; i < 3 ; i++) outString += String.format("%1.6e", arreArre.uyz[i]) + "\t";
            outString += String.format("%1.6e", arreArre.uyz[3]) + "\n";
            
            outString += "uxyz\t";
            for (int i = 0 ; i < 3 ; i++) outString += String.format("%1.6e", arreArre.uxyz[i]) + "\t";
            outString += String.format("%1.6e", arreArre.uxyz[3]) + "\n";
            
            outString += "uyxz\t";
            for (int i = 0 ; i < 3 ; i++) outString += String.format("%1.6e", arreArre.uyxz[i]) + "\t";
            outString += String.format("%1.6e", arreArre.uyxz[3]) + "\n";
            
           	//replace comma by period
        	String cOutString = outString.replace(',', '.');            	
        	            	
        	//write n flush
        	w.write(cOutString);
        	w.flush();
            
            //close file
            w.close();
            
        }catch(Exception e){
        
        	return true;
        	
        }
        
		return false;
        
    }		
	
	public boolean writeClusterMorphoResults(String sampleName, String path, MorphologyAnalyzer.PoreClusterProps mPCP) {
	       
		try{
            
			//IJ.error(path);
			
			//open file
			FileOutputStream fos = new FileOutputStream(path);
            Writer w = new BufferedWriter(new OutputStreamWriter(fos));
            
            //write header
        	String myHeader = "id\t" + "volume\t" + "surfaceArea\t" + "sphericity\t" + "euler\t";
        	myHeader += "connects2Top\t" + "connects2Bottom\t" + "percolates\t";
        	myHeader += "xCenter\t" + "yCenter\t" + "zCenter\t"; 
        	myHeader += "xmin\t" + "ymin\t" + "zmin\t" + "xmax\t" + "ymax\t" + "zmax\n";     	    
            w.write(myHeader);
                        
            //write the data        
            for (int i = 0 ; i < mPCP.id.length ; i++) {
            	
            	String outString = "";
            	outString += mPCP.id[i] + "\t";
            	
            	outString += String.format("%1.6e", mPCP.volume[i]) + "\t";
            	outString += String.format("%1.6e", mPCP.surfaceArea[i]) + "\t"; 
            	outString += String.format("%1.6e", mPCP.sphericity[i]) + "\t"; 
            	outString += String.format("%1.6e", mPCP.euler[i]) + "\t";
            	
            	outString += String.valueOf(mPCP.touchesTop[i]) + "\t";
            	outString += String.valueOf(mPCP.touchesBot[i]) + "\t";
            	outString += String.valueOf(mPCP.isPercolating[i]) + "\t";
            	
            	outString += String.format("%4.2f", mPCP.xCenter[i]) + "\t";
            	outString += String.format("%4.2f", mPCP.yCenter[i]) + "\t";
            	outString += String.format("%4.2f", mPCP.zCenter[i]) + "\t";	  
            	
            	outString += mPCP.minX[i] + "\t";
            	outString += mPCP.minY[i] + "\t";
            	outString += mPCP.minZ[i] + "\t";
	
            	outString += mPCP.maxX[i] + "\t";
            	outString += mPCP.maxY[i] + "\t";
            	outString += mPCP.maxZ[i] + "\n";
            	
            	//replace comma by period
            	String cOutString = outString.replace(',', '.');            	
            	            	
            	//write n flush
            	w.write(cOutString);
            	w.flush();
                       
            }
            
            //close file
            w.close();
            
        }catch(Exception e){
        
        	IJ.error("something went wrong when writing the cluster statistic results!");
        	return true;
        	
        }
        
		return false;
        
    }
	
	public boolean writeMRIFreezingResults(String name, String path, double[][] results, int times, int layers, int top, int bot) {
	       
		try{
            
			//IJ.error(path);
			String fileName = name + "_profile.txt";
			
			//open file
			FileOutputStream fos = new FileOutputStream(path + "\\" + fileName);
            Writer w = new BufferedWriter(new OutputStreamWriter(fos));
            
            //write header
        	String myHeader = "layer";
        	for (int i = 0 ; i < times ; i++) {
        		myHeader += "\t" + i;        		
        	}
        	myHeader += "\n";
        	
            w.write(myHeader);
                        
            //write the data        
            for (int z = 0 ; z < layers ; z++) {
            	
            	String outString = "" + (z + 1) + "\t";
            	
            	outString += "" + String.format("%5.2f",results[0][z]);
            	for (int i = 1 ; i < times ; i++) {
            		
            		outString += "\t" + String.format("%5.2f",results[i][z]);
            		
            	}
            	outString += "\n";
            	
            	//replace comma by period
            	String cOutString = outString.replace(',', '.');            	
            	            	
            	//write n flush
            	w.write(cOutString);
            	w.flush();
                       
            }
            
            //close file
            w.close();
            
        }
		
		catch(Exception e) {
        
        	IJ.error("something went wrong when writing the MRI statistic results!");
        	return true;
        	
        }
		
		//Also write global time series
		try{
            
			String fileName = name + "_global_" + top +"_" + bot + ".txt";
			
			//open file
			FileOutputStream fos = new FileOutputStream(path + "\\" + fileName);
            Writer w = new BufferedWriter(new OutputStreamWriter(fos));
            
            //write header
        	String myHeader = "T0";
        	for (int i = 1 ; i < times ; i++) {
        		myHeader += "\t" +"T" + i;        		
        	}
        	myHeader += "\n";
        	
            w.write(myHeader);
                        
            //write the data     
            String outString = "";
            for (int i = 0 ; i < times ; i++) {
            	
            	double[] profile = new double[layers];            
            		
            	for (int z = top ; z < bot ; z++) profile[z] = results[i][z];
            	
            	outString += String.format("%5.2f",StatUtils.mean(profile));
            	
            	if (i < times - 1) outString += "\t";
            	else outString += "\n";
            	
            }
            
        	//replace comma by period
        	String cOutString = outString.replace(',', '.');            	
        	            	
        	//write n flush
        	w.write(cOutString);
        	w.flush();
                               
            //close file
            w.close();
            
        }
		
		catch(Exception e) {
        
        	IJ.error("something went wrong when writing the MRI statistic results!");
        	return true;
        	
        }
		
        
		return false;
        
    }
	
	public boolean writeWRCResultsInASCII(MyFileCollection mFC, MenuWaiter.WRCCalculatorMenu mWRC, int[] air, int[] water, double[] aSM, double[] wSM) {
	       
		try{
			
			String pathSep = "/";
            
			//open file
			String path = mFC.myOutFolder + pathSep + mFC.colName + ".asc"; 
			
			FileOutputStream fos = new FileOutputStream(path);
            Writer w = new BufferedWriter(new OutputStreamWriter(fos));
            
            //write header
            String myPreHeader = "tensionAtTop\t" + "tensionAtCenter\t" + "tensionAtBottom\t" +
            		"airVolume (vx)\t" + "waterVolume (vx)\t" + "airSaturation\t" + "waterSaturation\n";
            w.write(myPreHeader);
            w.flush();
            
            //write results
            for (int i = 0 ; i < air.length ; i++) {
            	String integralString = String.format("%3.2f", mWRC.tensionAtTop[i]) + "\t" +
            							String.format("%3.2f", mWRC.tensionAtCenter[i]) + "\t" +
            							String.format("%3.2f", mWRC.tensionAtBottom[i]) + "\t" +
            							air[i] + "\t" + water[i] + "\t" +
            							String.format("%1.4f", aSM[i]) + "\t" +
            							String.format("%1.4f", wSM[i]) + "\n";
            
            	String cIntegralString = integralString.replace(',', '.');  
            	w.write(cIntegralString);
            	w.flush();
            }
            
            //close file
            w.close();
            
        } 
		catch (Exception e) {
        
        	return true;
        	
        }
        
		return false;
        
    }
	
	public boolean writeCutEllipsoidMaskFile(int i, String path, MyFileCollection mFC, MenuWaiter.ClipperMenuReturn mSCM) {
	 
		//load correct gaugefile
		int[] myGandS = getTheCorrectGaugeNSurfaceFiles(mFC.myTiffs[i], mFC.myInnerCircleFiles, mFC.mySurfaceFileNames);				
		String nowGaugePath = mFC.myInnerCircleFiles[myGandS[0]];
		
		//read InnerCircle file
		ObjectDetector jOD = new ObjectDetector();
		ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();
		int versio = checkInnerCircleFileVersion(nowGaugePath);			
		if (versio == 0) jCO = readInnerCircleVer0(nowGaugePath);	
		else jCO = readInnerCircleVer1(nowGaugePath);
		
		//modify inner circle appropriately
		jCO.topOfColumn = mSCM.startAtSlice - 1;
		jCO.bottomOfColumn = mSCM.stopAtSlice - 1;
		jCO.heightOfColumn = mSCM.heightOfROI;
		for (int j = 0 ; j < jCO.innerMajorRadius.length ; i++) {
			jCO.innerMajorRadius[j] -=(mSCM.clipFromCanvasEdge + mSCM.clipFromInnerPerimeter - mSCM.canvasExceedsBy);					
			jCO.outerMajorRadius[j] -=(mSCM.clipFromCanvasEdge + mSCM.clipFromInnerPerimeter - mSCM.canvasExceedsBy);
			jCO.innerMinorRadius[j] -=(mSCM.clipFromCanvasEdge + mSCM.clipFromInnerPerimeter - mSCM.canvasExceedsBy);					
			jCO.outerMinorRadius[j] -=(mSCM.clipFromCanvasEdge + mSCM.clipFromInnerPerimeter - mSCM.canvasExceedsBy);
		}
		
		//save the new inner circle
		boolean writersSuccess = writeInnerCircleVer0(path, jCO);
		
		return writersSuccess;
	}
	
	public boolean writeInnerCircleVer0(String path, ObjectDetector.ColCoords3D jCO) {
	       
		try{
            //open file
			FileOutputStream fos = new FileOutputStream(path);
            Writer w = new BufferedWriter(new OutputStreamWriter(fos));
 
            //write header 1
        	String myHeader0 = "tiltInXZ\t" + "tiltInYZ\t" + "tiltTotal\t" + "wallThickness\t";        	
        	myHeader0 = myHeader0 + "heightOfColumn\t" + "numberOfImputedLayers\n";
        	w.write(myHeader0);
        	w.flush();
        	
    		//write data1        	
        	String scalars = "";
    		scalars += String.format("%1.4f\t",jCO.tiltInXZ);
    		scalars += String.format("%1.4f\t",jCO.tiltInYZ);
    		scalars += String.format("%1.4f\t",jCO.tiltTotal);  
    		scalars += String.format("%d\t",(int)Math.round(StatUtils.percentile(jCO.wallThickness, 50)));  
    		scalars += String.format("%d\t",jCO.heightOfColumn);    		
    		scalars += String.format("%d\n",jCO.numberOfImputedLayers);
    		w.write(scalars);
    		w.flush();
    	
            //write header2
        	String myHeader = "xOutMid\t" + "yOutMid\t" + "zOutMid\t";        	
        	myHeader += "xInnMid\t" + "yInnMid\t";  
        	myHeader += "outerMajorRadius\t" + "outerMinorRadius\t" + "innerMajorRadius\t" + "innerMinorRadius\t"; 
        	myHeader += "wallThickness\t" + "outerTheta\t" + "innerTheta\t";
        	myHeader += "tiltinXZ\t" + "tiltinXY\t" + "tiltTotal\t";
        	myHeader += "outerR2\t" + "innerR2\n";
        	w.write(myHeader);
        	w.flush();
            
            //write the data2
            for (int i = jCO.topOfColumn ; i < jCO.topOfColumn + jCO.heightOfColumn ; i++) {
            	 
            	String outString = "";
            	outString += String.format("%4.2f\t", jCO.xmid[i]);
            	outString += String.format("%4.2f\t", jCO.ymid[i]);
            	outString += String.format("%4.2f\t", jCO.zmid[i]);
            	outString += String.format("%4.2f\t", jCO.ixmid[i]);
            	outString += String.format("%4.2f\t", jCO.iymid[i]);
            	outString += String.format("%4.2f\t", jCO.outerMajorRadius[i]);
            	outString += String.format("%4.2f\t", jCO.outerMinorRadius[i]);
            	outString += String.format("%4.2f\t", jCO.innerMajorRadius[i]);
            	outString += String.format("%4.2f\t", jCO.innerMinorRadius[i]);
            	outString += String.format("%3.2f\t", jCO.wallThickness[i]);
            	outString += String.format("%3.4f\t", jCO.theta[i]);
            	outString += String.format("%3.4f\t", jCO.itheta[i]);
            	outString += String.format("%2.4f\t", jCO.outerR2[i]);
            	outString += String.format("%2.4f\n", jCO.innerR2[i]);
            	
            	//replace comma by period
            	String cOutString = outString.replace(',', '.');            	
            	
            	//write
            	w.write(cOutString);
            	w.flush();
            }            
            
            //close file
            w.close();  
            
        }catch(Exception e){
                
        	return true;
        	
        }
        
		return false;
        
    }
	
	public boolean writeInnerCircleSteel(ObjectDetector.EggShapedColCoords3D eggshapedCC, MyFileCollection mFC) {
	       
		try{
			
			//definePath
			String path = mFC.myInnerCircleFolder + mFC.pathSep + "Steel" + mFC.fileName + ".txt"; 
			
            //open file
			FileOutputStream fos = new FileOutputStream(path);
            Writer w = new BufferedWriter(new OutputStreamWriter(fos));
 
            //write header 1
        	String myHeader0 = "heightOfColumn\n";
        	w.write(myHeader0);
        	w.flush();
        	
    		//write data1        	
        	String scalars = "" + eggshapedCC.heightOfColumn + "\n";
        	w.write(scalars);
    		w.flush();
    		
    		//write header 2
    		myHeader0 = "Z\tXcenter\tYcenter\tWallThickness\t'then each line xOD, yOD, xID, yID, i.e. outer and inner perimeter coordinates\n";
        	w.write(myHeader0);
        	w.flush();
        	
        	int top = eggshapedCC.topOfColumn;
        	int bot = eggshapedCC.bottomOfColumn;        	      	
        	for (int i = top ; i <= bot ; i++) {
        		scalars = String.format("%1.2f\t",eggshapedCC.zmid[i]);;
        		scalars += String.format("%1.2f\t",eggshapedCC.xmid[i]);
        		scalars += String.format("%1.2f\t",eggshapedCC.ymid[i]);
        		scalars += String.format("%1.2f\n",eggshapedCC.wallThickness[i]);
        		w.write(scalars);
        		w.flush();
        		
        		int j = 0;
        		scalars = "";
        		for (j = 0 ; j < eggshapedCC.anglesChecked - 1 ; j++) {        			
        			scalars += String.format("%1.2f\t",eggshapedCC.xOD[i][j]);
        		}
        		scalars += String.format("%1.2f\n",eggshapedCC.xOD[i][j + 1]);
        		w.write(scalars);
        		w.flush();
        		
        		scalars = "";
           		for (j = 0 ; j < eggshapedCC.anglesChecked - 1 ; j++) {        			
        			scalars += String.format("%1.2f\t",eggshapedCC.yOD[i][j]);
        		}
           		scalars += String.format("%1.2f\n",eggshapedCC.yOD[i][j + 1]);
        		w.write(scalars);
        		w.flush();
        		
        		scalars = "";
           		for (j = 0 ; j < eggshapedCC.anglesChecked - 1 ; j++) {        			
        			scalars += String.format("%1.2f\t",eggshapedCC.xID[i][j]);
        		}
           		scalars += String.format("%1.2f\n",eggshapedCC.xID[i][j + 1]);
        		w.write(scalars);
        		w.flush();
        		
        		scalars = "";
           		for (j = 0 ; j < eggshapedCC.anglesChecked - 1 ; j++) {        			
        			scalars += String.format("%1.2f\t",eggshapedCC.yID[i][j]);
        		}
           		scalars += String.format("%1.2f\n\n",eggshapedCC.yID[i][j + 1]);
        		w.write(scalars);
        		w.flush();
        
            }            
            
            //close file
            w.close();  
            
        }catch(Exception e){
                
        	return true;
        	
        }
        
		return false;
        
    }
	
	public boolean writeInnerCircleVer1(String path, ObjectDetector.ColCoords3D jCO) {
	       
		try{
            //open file
			FileOutputStream fos = new FileOutputStream(path);
            Writer w = new BufferedWriter(new OutputStreamWriter(fos));
 
            //write header 1
        	String myHeader0 = "heightOfColumn\n";        	
        	w.write(myHeader0);
        	w.flush();
        	
    		//write data1        	
        	String scalars = "";
    		scalars += String.format("%d\n",jCO.heightOfColumn);  
    		scalars = scalars.replace(',', '.');
    		w.write(scalars);
    		w.flush();
    		    	
            //write header2
        	String myHeader = "number\t" + "wallThickness\t";          	
        	myHeader += "zOutMid\t" + "xOutMid\t" + "yOutMid\t" +"xInnMid\t" + "yInnMid\t";  
        	myHeader += "outerMajorRadius\t" + "outerMinorRadius\t" + "innerMajorRadius\t" + "innerMinorRadius\t"; 
        	myHeader +=  "outerTheta\t" + "innerTheta\t";
        	myHeader += "outerR2\t" + "innerR2\n";
        	w.write(myHeader);
        	w.flush();
            
            //write the data2
            for (int i = 0 ; i < jCO.xmid.length ; i++) {
            	 
            	String outString = "";
            	outString += String.format("%4.0f\t", (float)(i + 1));
            	outString += String.format("%3.2f\t", jCO.wallThickness[i]);
            	outString += String.format("%4.2f\t", jCO.zmid[i]);
            	outString += String.format("%4.2f\t", jCO.xmid[i]);
            	outString += String.format("%4.2f\t", jCO.ymid[i]);
            	outString += String.format("%4.2f\t", jCO.ixmid[i]);
            	outString += String.format("%4.2f\t", jCO.iymid[i]);
            	outString += String.format("%4.2f\t", jCO.outerMajorRadius[i]);
            	outString += String.format("%4.2f\t", jCO.outerMinorRadius[i]);
            	outString += String.format("%4.2f\t", jCO.innerMajorRadius[i]);
            	outString += String.format("%4.2f\t", jCO.innerMinorRadius[i]);            	
            	outString += String.format("%3.4f\t", jCO.theta[i]);
            	outString += String.format("%3.4f\t", jCO.itheta[i]);
            	outString += String.format("%1.4f\t", jCO.outerR2[i]);
            	outString += String.format("%1.4f\n", jCO.innerR2[i]);
            	
            	//replace comma by period
            	String cOutString = outString.replace(',', '.');            	
            	
            	//write
            	w.write(cOutString);
            	w.flush();
            }            
            
            //close file
            w.close();  
            
        } catch(Exception e){
        	
        	IJ.error("I am terribly sorry, but something went wrong with saving the Inner Circle file.");
                
        	return true;
        	
        }
        
		return false;
        
    }
	
	/*public boolean writeSteelMaskFile(String path, ObjectDetector.ColumnCoordinates jCO) {
	       
		try{
            //open file
			FileOutputStream fos = new FileOutputStream(path);
            Writer w = new BufferedWriter(new OutputStreamWriter(fos));
 
            //write header 0
        	String myHeader00 = "approximateColumnHeight\t" + "approximateColumnRadius\t" + "coningStartsApproximatelyBeforeBottom\n";
        	myHeader00 += String.format("%5.0f\t",jCO.approximateHeight) + "\t" + String.format("%4.0f\t",jCO.approximateRadius) + "\t" + String.format("%3.0f\t",jCO.coningStartsBeforeEnd) + "\n";
        	w.write(myHeader00);
        	w.flush();
            
            //write header 1
        	String myHeader0 = "tiltInXZ\t" + "tiltInYZ\t" + "tiltTotal\t" + "topOfColumn\t" + "bottomOfColumn\t";        	
        	myHeader0 += "heightOfColumn\t" + "wallThickness\t" + "wallGreyValue\t" + "#ofAnglesChecked\n";
        	w.write(myHeader0);
        	w.flush();
        	
    		//write data1        	
        	String scalars = "";
    		scalars += String.format("%1.4f\t",jCO.tiltInXZ);
    		scalars += String.format("%1.4f\t",jCO.tiltInYZ);
    		scalars += String.format("%1.4f\t",jCO.tiltTotal);
    		scalars += String.format("%d\t",jCO.topOfColumn);
    		scalars += String.format("%d\t",jCO.bottomOfColumn);
       		scalars += String.format("%d\t",jCO.heightOfColumn); 
    		scalars += String.format("%3.0f\t",jCO.wallThickness);
    		scalars += String.format("%5.0f\t",jCO.steelGreyValue);
    		scalars += jCO.anglesChecked + "\n";
    		w.write(scalars);
    		w.flush();
    	
            //write header2
        	String myHeader = "layer\t" + "xOutMid\t" + "yOutMid\t" + "zOutMid\t" + "medianOfWallGreyValue\n";
        	w.write(myHeader);
        	w.flush();
            
            //write the data2
            for (int i = 0 ; i < jCO.xmid.length ; i++) {
            	 
            	String outString = i + "\t";
            	outString += String.format("%4.2f\t", jCO.xmid[i]);
            	outString += String.format("%4.2f\t", jCO.ymid[i]);
            	outString += String.format("%4.2f\t", jCO.zmid[i]);
            	outString += String.format("%5.0f\n", jCO.wallGreyValues[i]);
            	
            	w.write(outString);
            	w.flush();
            }            
            
            //write ROI coordinates
        	for (int j = 0 ; j < jCO.xmid.length ; j++) {
        		String outString = "\n";	//write one empty line first
        		outString += "Level " + j + "\n";
        		
        		//xOD
        		for (int i = 0 ; i < jCO.xOD[j].length - 1; i++) {
        			outString += String.format("%4.0f\t", jCO.xOD[j][i]);
        		}
        		outString += String.format("%4.0f\n", jCO.xOD[j][jCO.xOD[j].length - 1]);
        		
        		//yOD
        		for (int i = 0 ; i < jCO.yOD[j].length - 1; i++) {
        			outString += String.format("%4.0f\t", jCO.yOD[j][i]);
        		}
        		outString += String.format("%4.0f\n", jCO.yOD[j][jCO.yOD[j].length - 1]);
        		
        		//xID
        		for (int i = 0 ; i < jCO.xID[j].length - 1; i++) {
        			outString += String.format("%4.0f\t", jCO.xID[j][i]);
        		}
        		outString += String.format("%4.0f\n", jCO.xID[j][jCO.xID[j].length - 1]);
        		
        		//yID
        		for (int i = 0 ; i < jCO.yID[j].length - 1; i++) {
        			outString += String.format("%4.0f\t", jCO.yID[j][i]);
        		}
        		outString += String.format("%4.0f\n", jCO.yID[j][jCO.yID[j].length - 1]);     
        		
            	//replace comma by period
            	String cOutString = outString.replace(',', '.');            	
            	
            	//write        		
            	w.write(cOutString);
            	w.flush();
        	}
            
            //close file
            w.close();  
            
        }catch(Exception e){
                
        	return true;
        	
        }
        
		return false;
        
    }*/
	
	public int checkInnerCircleFileVersion(String nowGaugePath) {
		
		int cc = 0;
		String line;
		
		try{
            //open file
			FileReader fr = new FileReader(nowGaugePath);		
            BufferedReader br = new BufferedReader(fr);
	            
            //read header line 1
            line = br.readLine();
            
		} catch(Exception e) {return -1;}
		
		String checkString = line.substring(0, 6);
        if (checkString.equalsIgnoreCase("height")) return 1;
        else return 0;
		
	}
	
	public ObjectDetector.EggShapedColCoords3D readInnerCircleSteel(MyFileCollection mFC) {
		
		ObjectDetector jOD = new ObjectDetector();
		ObjectDetector.EggShapedColCoords3D jCO = jOD.new EggShapedColCoords3D();
		
		String scalar;
		String vectors;
		int cc = 0;
		
		try{
            //open file
			FileReader fr = new FileReader(mFC.nowInnerCirclePath);		
            BufferedReader br = new BufferedReader(fr);
		
            //init all data string
            StringBuilder sb = new StringBuilder();
            
            //skip header1
            String line = br.readLine();
            
            //read scalars
            scalar = br.readLine();
            
            //skip header2
            line = br.readLine();
            
            //read data
            while (line != null) {
            	if (cc > 0) {
            		sb.append(line + "\n");
            	}            	
            	line = br.readLine();
            	cc++;
            }
            vectors = sb.toString();
                 
            br.close();
            
		} catch(Exception e) {return null;}
		
		//parse data
		String corrScalars = scalar.replace(',','.');
		jCO.heightOfColumn = Integer.parseInt(corrScalars);
		
		//parse data2			
		int dataNumber = (cc - 1) / 6; 
		double[] xmid = new double[dataNumber];			//x midpoint
		double[] ymid = new double[dataNumber];			//y midpoint
		double[] zmid = new double[dataNumber];			//y midpoint
		double[] wallThickness = new double[dataNumber];
		double[][] xOD = new double[dataNumber][72];  //assuming that there are 72 checked angles..
		double[][] yOD = new double[dataNumber][72];  
		double[][] xID = new double[dataNumber][72];  
		double[][] yID = new double[dataNumber][72];  
		
		int[] lineBreaks = findLineBreaks(vectors, cc);
		cc = 0;
		for (int i = 0 ; i < lineBreaks.length - 1 ; i++) {
			
			String myLine0 = vectors.substring(lineBreaks[i], lineBreaks[i + 1]);
			String myLine1 = myLine0.replace('\n', ' ');
			String myLine = myLine1.replace(',', '.');
			
			int[] tabPos2 = findTabPositions(myLine);
			
			if (Math.floorMod(i, 6) == 0) {
				zmid[cc] = Double.parseDouble(myLine.substring(0, tabPos2[0]));
				xmid[cc] = Double.parseDouble(myLine.substring(tabPos2[0] + 1, tabPos2[1]));
				ymid[cc] = Double.parseDouble(myLine.substring(tabPos2[1] + 1, tabPos2[2]));
				wallThickness[cc] = Double.parseDouble(myLine.substring(tabPos2[2] + 1, myLine.length() - 1));
			}

			if (Math.floorMod(i, 6) == 1) {
				xOD[i][0] = Double.parseDouble(myLine.substring(0, tabPos2[0]));
				for (int j = 1 ; j < 71 ; j++) {
					xOD[i][j] = Double.parseDouble(myLine.substring(tabPos2[j] + 1, tabPos2[j + 1]));
				}
				xOD[i][72] = Double.parseDouble(myLine.substring(tabPos2[71] + 1, myLine.length() - 1));
			}

			if (Math.floorMod(i, 6) == 1) {
				yOD[i][0] = Double.parseDouble(myLine.substring(0, tabPos2[0]));
				for (int j = 1 ; j < 71 ; j++) {
					yOD[i][j] = Double.parseDouble(myLine.substring(tabPos2[j] + 1, tabPos2[j + 1]));
				}
				yOD[i][72] = Double.parseDouble(myLine.substring(tabPos2[71] + 1, myLine.length() - 1));
			}
			
			if (Math.floorMod(i, 6) == 1) {
				xID[i][0] = Double.parseDouble(myLine.substring(0, tabPos2[0]));
				for (int j = 1 ; j < 71 ; j++) {
					xID[i][j] = Double.parseDouble(myLine.substring(tabPos2[j] + 1, tabPos2[j + 1]));
				}
				xID[i][72] = Double.parseDouble(myLine.substring(tabPos2[71] + 1, myLine.length() - 1));
			}
			
			if (Math.floorMod(i, 6) == 1) {
				yID[i][0] = Double.parseDouble(myLine.substring(0, tabPos2[0]));
				for (int j = 1 ; j < 71 ; j++) {
					yID[i][j] = Double.parseDouble(myLine.substring(tabPos2[j] + 1, tabPos2[j + 1]));
				}
				yID[i][72] = Double.parseDouble(myLine.substring(tabPos2[71] + 1, myLine.length() - 1));
			}
					
			cc++;
		}

		//transfer readout to output structure
		jCO.xmid = xmid;
		jCO.ymid = ymid;
		jCO.zmid = zmid;
		jCO.wallThickness = wallThickness;
		
		jCO.xOD = xOD;
		jCO.yOD = yOD;
		jCO.xID = xID;
		jCO.yID = yID;
		
		return jCO;
		
	}
	
	public ObjectDetector.ColCoords3D readInnerCircleVer0(String nowGaugePath) {
		
		ObjectDetector jOD = new ObjectDetector();
		ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();
		
		int cc = 0;
		String scalars;
		String vectors;
		
		try{
            //open file
			FileReader fr = new FileReader(nowGaugePath);		
            BufferedReader br = new BufferedReader(fr);
		
            //init all data string
            StringBuilder sb = new StringBuilder();
            
            //skip header1
            String line = br.readLine();
            
            //read scalars
            scalars = br.readLine();
            
            //skip header2
            line = br.readLine();
            
            //read data    
            cc = 0;
            while (line != null) {
            	if (cc > 0) {
            		sb.append(line + "\n");
            	}            	
            	line = br.readLine();
            	cc++;
            }
            vectors = sb.toString();
                 
            br.close();
            
		} catch(Exception e) {return null;}
		
		//parse data
		String corrScalars = scalars.replace(',','.');
		int[] tabPos = findTabPositions(corrScalars);
		jCO.tiltInXZ = Double.parseDouble(corrScalars.substring(0, tabPos[0]));
		jCO.tiltInYZ = Double.parseDouble(corrScalars.substring(tabPos[0] + 1, tabPos[1]));
		jCO.tiltTotal = Double.parseDouble(corrScalars.substring(tabPos[1] + 1, tabPos[2]));	
		Integer.parseInt(corrScalars.substring(tabPos[2] + 1, tabPos[3]));	//average wall Thickness dummy
		jCO.heightOfColumn = Integer.parseInt(corrScalars.substring(tabPos[3] + 1, tabPos[4]));		
		jCO.numberOfImputedLayers = Integer.parseInt(corrScalars.substring(tabPos[4] + 1));
		
		//parse data2			
		double[] xmid = new double[cc - 1];			//x midpoint
		double[] ymid = new double[cc - 1];			//y midpoint
		double[] zmid = new double[cc - 1];			//y midpoint		
		double[] ixmid = new double[cc - 1];			//x midpoint (inner circle)
		double[] iymid = new double[cc - 1];			//y midpoint (inner circle)		
		double[] outerMajorRadius = new double[cc - 1];
		double[] innerMajorRadius = new double[cc - 1];
		double[] outerMinorRadius = new double[cc - 1];
		double[] innerMinorRadius = new double[cc - 1];
		double[] wallThickness = new double[cc - 1];
		double[] theta = new double[cc - 1];  //angle of major ellipse axis
		double[] itheta = new double[cc - 1];  //angle of major ellipse axis (inner circle)
		double[] outerR2 = new double[cc - 1];
		double[] innerR2 = new double[cc - 1];
		
		int[] lineBreaks = findLineBreaks(vectors, cc);
		cc = 0;
		for (int i = 0 ; i < lineBreaks.length - 1 ; i++) {
			
			String myLine0 = vectors.substring(lineBreaks[i], lineBreaks[i + 1]);
			String myLine1 = myLine0.replace('\n', ' ');
			String myLine = myLine1.replace(',', '.');
			
			int[] tabPos2 = findTabPositions(myLine);
			
			xmid[cc] = Double.parseDouble(myLine.substring(0, tabPos2[0]));
			ymid[cc] = Double.parseDouble(myLine.substring(tabPos2[0] + 1, tabPos2[1]));		
			zmid[cc] = Double.parseDouble(myLine.substring(tabPos2[1] + 1, tabPos2[2]));	
			ixmid[cc] = Double.parseDouble(myLine.substring(tabPos2[2] + 1, tabPos2[3]));	
			iymid[cc] = Double.parseDouble(myLine.substring(tabPos2[3] + 1, tabPos2[4]));	
			outerMajorRadius[cc] = Double.parseDouble(myLine.substring(tabPos2[4] + 1, tabPos2[5]));	
			outerMinorRadius[cc] = Double.parseDouble(myLine.substring(tabPos2[5] + 1, tabPos2[6]));	
			innerMajorRadius[cc] = Double.parseDouble(myLine.substring(tabPos2[6] + 1, tabPos2[7]));
			innerMinorRadius[cc] = Double.parseDouble(myLine.substring(tabPos2[7] + 1, tabPos2[8]));	
			wallThickness[cc] = Double.parseDouble(myLine.substring(tabPos2[8] + 1, tabPos2[9]));	
			theta[cc] = Double.parseDouble(myLine.substring(tabPos2[9] + 1, tabPos2[10]));
			itheta[cc] = Double.parseDouble(myLine.substring(tabPos2[10] + 1, tabPos2[11]));	
			outerR2[cc] = Double.parseDouble(myLine.substring(tabPos2[11] + 1, tabPos2[12]));	
			innerR2[cc] = Double.parseDouble(myLine.substring(tabPos2[12] + 1));
					
			cc++;
		}
			
		jCO.xmid = xmid;
		jCO.ymid = ymid;
		jCO.zmid = zmid;
		jCO.ixmid = ixmid;
		jCO.iymid = iymid;
		jCO.outerMajorRadius = outerMajorRadius;
		jCO.innerMajorRadius = innerMajorRadius;
		jCO.outerMinorRadius = outerMinorRadius;
		jCO.innerMinorRadius = innerMinorRadius;
		jCO.wallThickness = wallThickness;
		jCO.theta = theta;
		jCO.itheta = itheta;
		jCO.outerR2 = outerR2;
		jCO.innerR2 = innerR2;
		
		return jCO;
		
	}
	
	public ObjectDetector.ColCoords3D readInnerCircleVer1(String nowGaugePath) {
		
		ObjectDetector jOD = new ObjectDetector();
		ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();
		
		int cc = 0;
		String scalars;
		String vectors;
				
		try{
            //open file
			FileReader fr = new FileReader(nowGaugePath);		
            BufferedReader br = new BufferedReader(fr);
		
            //init all data string
            StringBuilder sb = new StringBuilder();
            
            //skip header1
            String line = br.readLine();
            
            //read scalars
            scalars = br.readLine();
            
            //skip header2
            line = br.readLine();
            
            //read data    
            cc = 0;
            while (line != null) {
            	if (cc > 0) {
            		sb.append(line + "\n");
            	}            	
            	line = br.readLine();
            	cc++;
            }
            vectors = sb.toString();
                 
            br.close();
            
		} catch(Exception e) {return null;}
		
		//parse data
		String corrScalars = scalars.replace(',','.');		
		jCO.heightOfColumn = Integer.parseInt(corrScalars);		
		
		//parse data2			
		double[] xmid = new double[cc - 1];			//x midpoint
		double[] ymid = new double[cc - 1];			//y midpoint
		double[] zmid = new double[cc - 1];			//y midpoint		
		double[] ixmid = new double[cc - 1];			//x midpoint (inner circle)
		double[] iymid = new double[cc - 1];			//y midpoint (inner circle)		
		double[] outerMajorRadius = new double[cc - 1];
		double[] innerMajorRadius = new double[cc - 1];
		double[] outerMinorRadius = new double[cc - 1];
		double[] innerMinorRadius = new double[cc - 1];
		double[] wallThickness = new double[cc - 1];
		double[] theta = new double[cc - 1];  //angle of major ellipse axis
		double[] itheta = new double[cc - 1];  //angle of major ellipse axis (inner circle)
		double[] outerR2 = new double[cc - 1];
		double[] innerR2 = new double[cc - 1];
		
		int[] lineBreaks = findLineBreaks(vectors, cc);
		cc = 0;
		for (int i = 0 ; i < lineBreaks.length - 1 ; i++) {
			
			String myLine0 = vectors.substring(lineBreaks[i], lineBreaks[i + 1]);
			String myLine1 = myLine0.replace('\n', ' ');
			String myLine = myLine1.replace(',', '.');
			
			int[] tabPos2 = findTabPositions(myLine);
	
			wallThickness[cc] = Double.parseDouble(myLine.substring(tabPos2[0] + 1, tabPos2[1]));	
			
			zmid[cc] = Double.parseDouble(myLine.substring(tabPos2[1] + 1, tabPos2[2]));
			xmid[cc] = Double.parseDouble(myLine.substring(tabPos2[2] + 1, tabPos2[3]));
			ymid[cc] = Double.parseDouble(myLine.substring(tabPos2[3] + 1, tabPos2[4]));
			ixmid[cc] = Double.parseDouble(myLine.substring(tabPos2[4] + 1, tabPos2[5]));	
			iymid[cc] = Double.parseDouble(myLine.substring(tabPos2[5] + 1, tabPos2[6]));		
					
			outerMajorRadius[cc] = Double.parseDouble(myLine.substring(tabPos2[6] + 1, tabPos2[7]));	
			outerMinorRadius[cc] = Double.parseDouble(myLine.substring(tabPos2[7] + 1, tabPos2[8]));	
			innerMajorRadius[cc] = Double.parseDouble(myLine.substring(tabPos2[8] + 1, tabPos2[9]));
			innerMinorRadius[cc] = Double.parseDouble(myLine.substring(tabPos2[9] + 1, tabPos2[10]));	
			
			theta[cc] = Double.parseDouble(myLine.substring(tabPos2[10] + 1, tabPos2[11]));
			itheta[cc] = Double.parseDouble(myLine.substring(tabPos2[11] + 1, tabPos2[12]));	
			outerR2[cc] = Double.parseDouble(myLine.substring(tabPos2[12] + 1, tabPos2[13]));	
			innerR2[cc] = Double.parseDouble(myLine.substring(tabPos2[13] + 1));
					
			cc++;
		}
			
		jCO.xmid = xmid;
		jCO.ymid = ymid;
		jCO.zmid = zmid;
		jCO.ixmid = ixmid;
		jCO.iymid = iymid;
		jCO.outerMajorRadius = outerMajorRadius;
		jCO.innerMajorRadius = innerMajorRadius;
		jCO.outerMinorRadius = outerMinorRadius;
		jCO.innerMinorRadius = innerMinorRadius;
		jCO.wallThickness = wallThickness;
		jCO.theta = theta;
		jCO.itheta = itheta;
		jCO.outerR2 = outerR2;
		jCO.innerR2 = innerR2;
		
		return jCO;
		
	}
	
	/*public ObjectDetector.ColumnCoordinates readSteelGaugeFile(String nowGaugePath) {
		
		ObjectDetector jOD = new ObjectDetector();
		ObjectDetector.ColumnCoordinates jCO = jOD.new ColumnCoordinates();
		
		int cc = 0;
		int cc1 = 0;
		int ccAll = 0;
		String scalars;
		String vectors;	
		String walls;
		
		try{
            //open file
			FileReader fr = new FileReader(nowGaugePath);		
            BufferedReader br = new BufferedReader(fr);
		
            //init data string 4 basic data
            StringBuilder sb = new StringBuilder();
            
            //skip header1
            br.readLine();
            
            //skip scalars0
            br.readLine();
            
            //skip header2
            br.readLine();
            
            //read scalars1
            scalars = br.readLine();
            
            //skip header2
            String line = br.readLine();            
            
            //read basic data    
            cc = 0;
            while (cc1 == 0) {
            	if (cc > 0) {
            		sb.append(line + "\n");         
            	}            	
            	line = br.readLine();
             	
            	if (cc1 == 0 & line.length() == 0) cc1 = cc; //remember when block one is over..
            	
            	cc++;            	
            }
            
            //save collected string for analysis and free sb space for reuse 
            vectors = sb.toString();
            sb.delete(0, sb.toString().length());
            sb.trimToSize();
            
            //parseBasicData
            jCO = parseBasicData(scalars, vectors, cc, cc1);
            
            //init wall position matrices
            float[][] xOD = new float[jCO.xmid.length][jCO.anglesChecked];
            float[][] xID = new float[jCO.xmid.length][jCO.anglesChecked];
            float[][] yOD = new float[jCO.xmid.length][jCO.anglesChecked];
            float[][] yID = new float[jCO.xmid.length][jCO.anglesChecked];
                        
            //read steel wall coordinates
            cc=0;
            while (line != null) {            	
            	if (!line.isEmpty()) {            		           		
            		sb.append(line + "\n");  
            		            		           		
            		if (cc == 4) {         
            			
            			//init parsing of next depth
            			cc = 0;
            			walls = sb.toString();
            			walls.replace(' ','0');
            			sb.delete(0,sb.toString().length());
            			sb.trimToSize();
            			
            			//parse wall coords
            			float[][] wallCoords = parseWallCoords(walls, jCO);
            			
            			//transfer wallCoords to jCO transferable output            			
            			for (int j = 0 ; j < jCO.anglesChecked ; j++) {            			
            				xOD[ccAll][j] = wallCoords[0][j];
            				yOD[ccAll][j] = wallCoords[1][j];
            				xID[ccAll][j] = wallCoords[2][j];
            				yID[ccAll][j] = wallCoords[3][j];
            			}
            			
            			ccAll++;
            			
            		}
            		else cc++;
            	}

            	line = br.readLine();
            	
            }
                             
            br.close();
            
            //transfer wall choords to jCO
    		jCO.xOD = xOD;
    		jCO.yOD = yOD;
    		jCO.xID = xID;
    		jCO.yID = yID;		
    		
    		return jCO;
            
            
		} catch(Exception e) {			
			return null;
		}
		
	}*/
	
	public void saveSurfaceStatistics(String myOutPath, String[] myTiffs, MorphologyAnalyzer.SurfaceStatistics[] mSS) {
		
		String pathSep = "/";
		
		try{
			
			//open file
			String path = myOutPath + pathSep + "SurfaceStatistics.txt";
	    	FileOutputStream fos = new FileOutputStream(path);
	    	Writer w = new BufferedWriter(new OutputStreamWriter(fos));
			
			//writeHeader
			String outString = "ColumnName\t" + "maxTopSurface\t" + "medianTopSurface\t" + "meanTopSurface\t" + "minTopSurface\t" +
					"maxBottomSurface\t" + "medianBottomSurface\t" + "meanBottomSurface\t" + "minBottompSurface\n";
			w.write(outString);
			
			//write data
			//int i = 1;
			for (int i = 0 ; i < mSS.length ; i++) {
				
				String name = myTiffs[i].substring(0, myTiffs[i].length() - 5);
				
				String dataString = "";
				dataString += name + "\t";
				dataString += mSS[i].highestElevation  + "\t";
				dataString += mSS[i].medianElevation  + "\t";
				dataString += mSS[i].meanElevation  + "\t";
				dataString += mSS[i].lowestElevation  + "\t";
				
				dataString += mSS[i].highestIntrusion  + "\t";
				dataString += mSS[i].medianIntrusion  + "\t";
				dataString += mSS[i].meanIntrusion  + "\t";
				dataString += mSS[i].lowestIntrusion  + "\n";
				
				w.write(dataString);
		    	w.flush();
				
			}
				
            //close file
            w.close();  
            
            
        }catch(Exception e){        	
        	return;        	
        }
		
		
	}
	
	public void writeRadialGreyValues(double[][][] radialGreyValues,  String nowGaugePath, int slices, int anglesChecked, int standardRadius) { 
				
		int maxAlpha = 360;
		int dAlpha = maxAlpha / anglesChecked;
		Percentile jP = new Percentile();
		
		//break down grey-values to percentiles..
		int[][][] q = new int[slices][standardRadius][4];
		 
		for (int i = 0 ; i < slices; i++) {
			for (int j = 0 ; j < standardRadius ; j++) {
				
				IJ.showStatus("Calculating percentiles #" + (i + 1) + "/" + slices);
				
				int angleCounter = 0;
				double[] thisRadius = new double[anglesChecked];
				for (double angle = 0 ; angle < 2 * Math.PI - Math.PI/400 ; angle = angle + 2 * Math.PI / (maxAlpha/dAlpha)) {	
					thisRadius[angleCounter] = radialGreyValues[i][angleCounter][j];
					angleCounter++;
				}
				int q01 = (int)Math.round(jP.evaluate(thisRadius, 1	));
				int q40 = (int)Math.round(jP.evaluate(thisRadius, 40));
				int q60 = (int)Math.round(jP.evaluate(thisRadius, 60));
				int q80 = (int)Math.round(jP.evaluate(thisRadius, 80));
				
				q[i][j][0] = q01;
				q[i][j][1] = q40;
				q[i][j][2] = q60;
				q[i][j][3] = q80;
			}
		}
		
		try{			
			
            //write the data2
            for (int k = 0 ; k < 4 ; k++)  {
            	
            	IJ.showStatus("Writing File #" + (k + 1) + "/" + 4);
            	
            	//init output path
            	String badReggae = "InnerCircle";
            	int startPosition = 0;
            	for (int i = 0 ; i < nowGaugePath.length() - badReggae.length() ; i++) {
            		String nowSub = nowGaugePath.substring(i,i+11);
            		if (nowSub.equals(badReggae)) {
            			startPosition = i + 11;
            		}
            	}
            	String path = nowGaugePath.substring(0, startPosition+1) + "q01_" + nowGaugePath.substring(startPosition + 7);

				if (k==1) path = nowGaugePath.substring(0, startPosition+1) + "q40_" + nowGaugePath.substring(startPosition + 7);
				if (k==2) path = nowGaugePath.substring(0, startPosition+1) + "q60_" + nowGaugePath.substring(startPosition + 7);
				if (k==3) path = nowGaugePath.substring(0, startPosition+1) + "q80_" + nowGaugePath.substring(startPosition + 7);
				
            	//open file
            	FileOutputStream fos = new FileOutputStream(path);
            	Writer w = new BufferedWriter(new OutputStreamWriter(fos));
                        	
            	for (int i = 0 ; i < slices ; i++) {
            	 
            		String outString = String.format("%d\t",q[i][0][k]);
            	            	            	
            		int j;
            		for (j = 1 ; j < standardRadius - 2 ; j++) {
            			outString += String.format("%d\t",q[i][j][k]);
            		}
            		outString += String.format("%d\n",q[i][j+1][k]);
            	
                	//replace comma by period
                	String cOutString = outString.replace(',', '.');            	
                	
                	//write            		
                	w.write(cOutString);
            		w.flush();
            	}
            	
            	//close file
                w.close();  
            }            
            
            
        }catch(Exception e){        	
        	return;        	
        }
	}
	
	public double[][][] readSteelBeamHardeningCorrectionParameters(String nowGaugePath, int slices) {
		
		double[][][] bhc = new double[slices][6][4];
		String outString;		
				
		try{			
			
            //write the data2
            for (int k = 0 ; k < 4 ; k += 3)  {
            	
            	//init output path
            	String badReggae = "InnerCircle";
            	int startPosition = 0;
            	for (int i = 0 ; i < nowGaugePath.length() - badReggae.length() ; i++) {
            		String nowSub = nowGaugePath.substring(i,i+11);
            		if (nowSub.equals(badReggae)) {
            			startPosition = i + 11;
            		}
            	}
            	
            	//init output path
            	String path = nowGaugePath.substring(0, startPosition + 1) + "c01_" + nowGaugePath.substring(startPosition + 7);

				if (k==1) path = nowGaugePath.substring(0, startPosition + 1) + "c40_" + nowGaugePath.substring(startPosition + 7);
				if (k==2) path = nowGaugePath.substring(0, startPosition + 1) + "c60_" + nowGaugePath.substring(startPosition + 7);
				if (k==3) path = nowGaugePath.substring(0, startPosition + 1) + "c80_" + nowGaugePath.substring(startPosition + 7);
				
	            //open file
				FileReader fr = new FileReader(path);		
	            BufferedReader br = new BufferedReader(fr);
			
	            //init data string 4 basic data
	            StringBuilder sb = new StringBuilder();
                    
	            int i = 0;
            	for (i = 0 ; i < slices ; i++) {            		
            		String line = br.readLine();
            		sb.append(line + "\n");            		
            	}
            	
            	outString = sb.toString();
                sb.delete(0, sb.toString().length());
                sb.trimToSize();
                br.close();   
                
                //parse data
        		double a, b, c, r, dy, R2;
        		int[] lineBreaks = findLineBreaks(outString, slices);   
        		for (int j = 0 ; j < slices ; j++) {
        			
        			String myLine0 = null;
        			if (j == slices - 1) myLine0 = outString.substring(lineBreaks[j]);
        			else myLine0 = outString.substring(lineBreaks[j], lineBreaks[j + 1]);
        			String myLine1 = myLine0.replace('\n', ' ');
        			String myLine = myLine1.replace(',', '.');
        			int[] tabPos2 = findTabPositions(myLine);
        			
        			a = Double.parseDouble(myLine.substring(0, tabPos2[0]));
        			b = Double.parseDouble(myLine.substring(tabPos2[0], tabPos2[1]));
        			c = Double.parseDouble(myLine.substring(tabPos2[1], tabPos2[2]));
        			r = Double.parseDouble(myLine.substring(tabPos2[2] + 1, tabPos2[3]));		
        			dy = Double.parseDouble(myLine.substring(tabPos2[3] + 1, tabPos2[4]));
        			R2 = Double.parseDouble(myLine.substring(tabPos2[4] + 1));
        			
        			bhc[j][0][k]=a;
        			bhc[j][1][k]=b;
        			bhc[j][2][k]=c;
        			bhc[j][3][k]=r;
        			bhc[j][4][k]=dy;
        			bhc[j][5][k]=R2;
        		}
                        		
     		}
            
		} catch(Exception e) {return null;}
		
		return bhc;
		
	}
	
	public double[][][] readPVCBeamHardeningCorrectionParameters(String nowGaugePath, int slices) {
		
		double[][][] bhc = new double[slices][7][4];
		String outString;		
				
		try{			
			
            //write the data2
            for (int k = 2 ; k < 4 ; k++)  {
            	
            	//init output path
            	String badReggae = "InnerCircle";
            	int startPosition = 0;
            	for (int i = 0 ; i < nowGaugePath.length() - badReggae.length() ; i++) {
            		String nowSub = nowGaugePath.substring(i,i+11);
            		if (nowSub.equals(badReggae)) {
            			startPosition = i + 11;
            		}
            	}
            	
            	//init output path
            	String path = nowGaugePath.substring(0, startPosition + 1) + "c01_" + nowGaugePath.substring(startPosition + 7);

				if (k==1) path = nowGaugePath.substring(0, startPosition + 1) + "c40_" + nowGaugePath.substring(startPosition + 7);
				if (k==2) path = nowGaugePath.substring(0, startPosition + 1) + "c60_" + nowGaugePath.substring(startPosition + 7);
				if (k==3) path = nowGaugePath.substring(0, startPosition + 1) + "c80_" + nowGaugePath.substring(startPosition + 7);
				
	            //open file
				FileReader fr = new FileReader(path);		
	            BufferedReader br = new BufferedReader(fr);
			
	            //init data string 4 basic data
	            StringBuilder sb = new StringBuilder();
                    
	            int i = 0;
            	for (i = 0 ; i < slices ; i++) {            		
            		String line = br.readLine();
            		sb.append(line + "\n");           		
            	}
            	
            	outString = sb.toString();
                sb.delete(0, sb.toString().length());
                sb.trimToSize();
                br.close();   
                
                //parse data
        		double my, sd, thresh, skiplast, r, dy, R2;
        		int[] lineBreaks = findLineBreaks(outString, slices);   
        		for (int j = 0 ; j < slices ; j++) {
        			
        			String myLine0 = null;
        			if (j == slices - 1) myLine0 = outString.substring(lineBreaks[j]);
        			else myLine0 = outString.substring(lineBreaks[j], lineBreaks[j + 1]);
        			String myLine1 = myLine0.replace('\n', ' ');
        			String myLine = myLine1.replace(',', '.');
        			int[] tabPos2 = findTabPositions(myLine);
        			
        			my = Double.parseDouble(myLine.substring(0, tabPos2[0]));
        			sd = Double.parseDouble(myLine.substring(tabPos2[0], tabPos2[1]));
        			thresh = Double.parseDouble(myLine.substring(tabPos2[1], tabPos2[2]));
        			skiplast = Double.parseDouble(myLine.substring(tabPos2[2], tabPos2[3]));
        			r = Double.parseDouble(myLine.substring(tabPos2[3] + 1, tabPos2[4]));		
        			dy = Double.parseDouble(myLine.substring(tabPos2[4] + 1, tabPos2[5]));
        			R2 = Double.parseDouble(myLine.substring(tabPos2[5] + 1));
        			
        			bhc[j][0][k]=my;
        			bhc[j][1][k]=sd;
        			bhc[j][2][k]=thresh;
        			bhc[j][3][k]=skiplast;
        			bhc[j][4][k]=r;
        			bhc[j][5][k]=dy;
        			bhc[j][6][k]=R2;
        		}
                        		
     		}
            
		} catch(Exception e) {return null;}
		
		return bhc;
		
	}
	
	/*public MorphologyAnalyzer.PoreClusterProperties readPoreClusterProperties(String poreClusterPropertiesFilePath) {
		
		MorphologyAnalyzer morph = new MorphologyAnalyzer();
		MorphologyAnalyzer.PoreClusterProperties mCP = morph.new PoreClusterProperties();
				
		int i = 0, j = 0;
		int cc = -1;
		String errString = null;
		
		//find out how long the file is		
		try {		
			BufferedReader pr = new BufferedReader(new FileReader(poreClusterPropertiesFilePath));
			for(String line; (line = pr.readLine()) != null; ) {
		      if (!line.startsWith("0")) cc++;
		    }
		    pr.close();		    
		} catch(Exception e) {			
			IJ.error("While reading the pore cluster properties..", "Could not find out the number of lines in the ascii file!");
			return null;
		}
		
		//cluster properties
		int[] id = new int[cc];
		double[] volume = new 		double[cc];	
		double[] xCenter = new double[cc];	
		double[] yCenter = new double[cc];		
		double[] zCenter = new double[cc];		
		double[] momentOfInertiaShortestAxis = new double[cc];		
		double[] momentOfInertiamiddleAxis = new double[cc];		
		double[] momentOfInertiaLongestAxis = new double[cc];		
		double[] unitVectorInXDirection = new double[cc];		
		double[] unitVectorInYDirection = new double[cc];		
		double[] unitVectorInZDirection = new double[cc];		
		double[] euler = new double[cc];
		double[] holes = new double[cc];
		double[] cavities = new double[cc];
		double[] thickness = new double[cc];	
		double[] SDThickness = new double[cc];		
		double[] maxThickness = new double[cc];			
		
		try {
			 //open file
			FileReader fr = new FileReader(poreClusterPropertiesFilePath);		
            BufferedReader br = new BufferedReader(fr);
            
            //read header line
        	String line = br.readLine();
        	
        	IJ.showStatus("Reading properties of pore clusters ..");
        	
        	//read all the rest
        	for (i = 0 ; i < cc ; i++) {
        		
        		line = br.readLine();
        		String corrLine = line.replace(',','.');
        		int[] tabPos = findTabPositions(corrLine);	
        		
        		if (!line.startsWith("0")) {
        			for (j = 0 ; j < tabPos.length ; j++) {
	        		
	        			id[i] = Integer.parseInt(corrLine.substring(0, tabPos[0]));
	        			volume[i] = Double.parseDouble(corrLine.substring(tabPos[0], tabPos[1]));	
	        			xCenter[i] = Double.parseDouble(corrLine.substring(tabPos[1], tabPos[2]));		
	        			yCenter[i] = Double.parseDouble(corrLine.substring(tabPos[2], tabPos[3]));		
	        			zCenter[i] = Double.parseDouble(corrLine.substring(tabPos[3], tabPos[4]));		
	        			momentOfInertiaShortestAxis[i] = Double.parseDouble(corrLine.substring(tabPos[4], tabPos[5]));		
	        			momentOfInertiamiddleAxis[i] = Double.parseDouble(corrLine.substring(tabPos[5], tabPos[6]));		
	        			momentOfInertiaLongestAxis[i] = Double.parseDouble(corrLine.substring(tabPos[6], tabPos[7]));		
	        			unitVectorInXDirection[i] = Double.parseDouble(corrLine.substring(tabPos[7], tabPos[8]));		
	        			unitVectorInYDirection[i] = Double.parseDouble(corrLine.substring(tabPos[8], tabPos[9]));		
	        			unitVectorInZDirection[i] = Double.parseDouble(corrLine.substring(tabPos[9], tabPos[10]));		
	        			euler[i] = Double.parseDouble(corrLine.substring(tabPos[10], tabPos[11]));
	        			holes[i] = Double.parseDouble(corrLine.substring(tabPos[11], tabPos[12]));
	        			cavities[i] = Double.parseDouble(corrLine.substring(tabPos[12], tabPos[13]));
	        			thickness[i] = Double.parseDouble(corrLine.substring(tabPos[13], tabPos[14]));	
	        			SDThickness[i] = Double.parseDouble(corrLine.substring(tabPos[14], tabPos[15]));		
	        			maxThickness[i] = Double.parseDouble(corrLine.substring(tabPos[15]));	
	        			
	        		}        	
        		}
        	}
			
            br.close();            
		} catch(Exception e) {		
			errString = e.getMessage();
			return null;
		} finally {
			if (i != cc) IJ.error("While reading the pore cluster properties..", "Error reading line " + (i + 1) + "\n" 
						+ "Computer says '" + errString + "'!");
		}
		
		
		mCP.id = id;
		mCP.volume = volume;	
		mCP.xCenter = xCenter;	
		mCP.yCenter = yCenter;	
		mCP.zCenter = zCenter;	
		mCP.momentOfInertiaShortestAxis = momentOfInertiaShortestAxis;	
		mCP.momentOfInertiamiddleAxis = momentOfInertiamiddleAxis;	
		mCP.momentOfInertiaLongestAxis = momentOfInertiaLongestAxis;	
		mCP.unitVectorInXDirection = unitVectorInXDirection;	
		mCP.unitVectorInYDirection = unitVectorInYDirection;	
		mCP.unitVectorInZDirection = unitVectorInZDirection;	
		mCP.euler = euler;	
		mCP.holes = holes;
		mCP.cavities = cavities;
		mCP.thickness = thickness;
		mCP.sdThickness = SDThickness;	
		mCP.maxThickness = maxThickness;
				
		return mCP;
	}*/
	
	public int findNumberOfTabsInString(String myString) {
		
		int numberOfTabPositions = 0;
		for (int i = 0 ; i < myString.length() ; i++) {
			if (myString.charAt(i) == '\t') {				
				numberOfTabPositions++;
			}
		}
		
		return numberOfTabPositions;
	}
	
	public int[] findTabPositions(String myString) {
		
		int numberOfTabPositions = findNumberOfTabsInString(myString);		
		
		int[] tabPos = new int[numberOfTabPositions];
		
		int cc = 0;
		for (int i = 0 ; i < myString.length() ; i++) {
			if (myString.charAt(i) == '\t') {
				tabPos[cc] = i;
				cc++;
			}
		}
		
		return tabPos;
	}
	
	/*public ObjectDetector.ColumnCoordinates parseBasicData(String scalars, String vectors, int cc, int cc1) {		
		
		ObjectDetector jOD = new ObjectDetector();
		ObjectDetector.ColumnCoordinates jCO = jOD.new ColumnCoordinates();
		
		//parse data1
		String corrScalars = scalars.replace(',','.');
		int[] tabPos = findTabPositions(corrScalars);		
		jCO.tiltInXZ = Double.parseDouble(corrScalars.substring(0, tabPos[0]));
		jCO.tiltInYZ = Double.parseDouble(corrScalars.substring(tabPos[0] + 1, tabPos[1]));
		jCO.tiltTotal = Double.parseDouble(corrScalars.substring(tabPos[1] + 1, tabPos[2]));
		jCO.topOfColumn = Integer.parseInt(corrScalars.substring(tabPos[2] + 1, tabPos[3]));
		jCO.bottomOfColumn = Integer.parseInt(corrScalars.substring(tabPos[3] + 1, tabPos[4]));
		jCO.heightOfColumn = Integer.parseInt(corrScalars.substring(tabPos[4] + 1, tabPos[5]));
		jCO.steelWallThickness = Double.parseDouble(corrScalars.substring(tabPos[5] + 1, tabPos[6]));
		jCO.steelGreyValue = Double.parseDouble(corrScalars.substring(tabPos[6] + 1, tabPos[7])); 
		jCO.anglesChecked = Integer.parseInt(corrScalars.substring(tabPos[7] + 1));
			  		
		//parse data2
		int[] level = new int[cc1];
		double[] xmid = new double[cc1];			//x midpoint
		double[] ymid = new double[cc1];			//y midpoint
		double[] zmid = new double[cc1];			//y midpoint
		double[] wallGrey = new double[cc];
		int[] wallGreyValues = new int[cc1];
		
		int[] lineBreaks = findLineBreaks(vectors, cc);
		cc = 0;
		for (int i = 0 ; i < cc1 ; i++) {
			
			String myLine0 = vectors.substring(lineBreaks[i], lineBreaks[i + 1]);
			String myLine1 = myLine0.replace('\n', ' ');
			String myLine = myLine1.replace(',', '.');
			if (myLine.startsWith(" ")) myLine = myLine .substring(1);
			
			int[] tabPos2 = findTabPositions(myLine);
			
			level[cc] = Integer.parseInt(myLine.substring(0, tabPos2[0]));
			xmid[cc] = Double.parseDouble(myLine.substring(tabPos2[0], tabPos2[1]));
			ymid[cc] = Double.parseDouble(myLine.substring(tabPos2[1] + 1, tabPos2[2]));		
			zmid[cc] = Double.parseDouble(myLine.substring(tabPos2[2] + 1, tabPos2[3]));
			String theLastBit = myLine.substring(tabPos2[3] + 1);
			if (theLastBit.startsWith(" ")) theLastBit = theLastBit.substring(1);
			if (theLastBit.endsWith("\r")) theLastBit = theLastBit.substring(0,theLastBit.length() - 1);
			wallGreyValues[cc] = Integer.parseInt(theLastBit);
					
			cc++;
		}
	
		jCO.xmid = xmid;
		jCO.ymid = ymid;
		jCO.zmid = zmid;
		for (int fuck = 0 ; fuck < wallGreyValues.length ; fuck++) wallGrey[fuck] = wallGreyValues[fuck]; //Java sucks!!!!!!!!!!!!!!!!!
		jCO.wallGreyValues = wallGrey;
		
		return jCO;
	
	}
	*/
	public float[][] parseWallCoords(String vectors, ObjectDetector.ColCoords3D jCO) {
		
		//read level		
		int[] lineBreaks = findLineBreaks(vectors, 5);
		String myLine0 = vectors.substring(0, lineBreaks[1]);
		
		//parse data3	
		float[][] wallPositions = new float[4][jCO.anglesChecked];
		
		//parse steel wall coordinates		
		//XOD
		myLine0 = vectors.substring(lineBreaks[1], lineBreaks[2]);			
		String myLine1 = myLine0.replace('\n', ' ');
		String myLine = myLine1.replace(',', '.');
		myLine = myLine.replace(' ', '0');	
		int[] tabPos2 = findTabPositions(myLine);				
		wallPositions[0][0] = Integer.parseInt(myLine.substring(0, tabPos2[0]));
		for (int k = 1 ; k < jCO.anglesChecked; k++ ) {
			String nowLine = myLine.substring(tabPos2[k-1] + 1, tabPos2[k]);
			wallPositions[0][k] = Integer.parseInt(nowLine);
		}

		//YOD
		myLine0 = vectors.substring(lineBreaks[2], lineBreaks[3]);			
		myLine1 = myLine0.replace('\n', ' ');
		myLine = myLine1.replace(',', '.');
		myLine = myLine.replace(' ', '0');
		tabPos2 = findTabPositions(myLine);				
		wallPositions[1][0] = Integer.parseInt(myLine.substring(0, tabPos2[0]));
		for (int k = 1 ; k < jCO.anglesChecked; k++ ) wallPositions[1][k] = Integer.parseInt(myLine.substring(tabPos2[k-1] + 1, tabPos2[k]));

		//XID
		myLine0 = vectors.substring(lineBreaks[3], lineBreaks[4]);			
		myLine1 = myLine0.replace('\n', ' ');
		myLine = myLine1.replace(',', '.');		
		myLine = myLine.replace(' ', '0');
		tabPos2 = findTabPositions(myLine);				
		wallPositions[2][0] = Integer.parseInt(myLine.substring(0, tabPos2[0]));
		for (int k = 1 ; k < jCO.anglesChecked; k++ ) wallPositions[2][k] = Integer.parseInt(myLine.substring(tabPos2[k-1] + 1, tabPos2[k]));

		//YID
		myLine0 = vectors.substring(lineBreaks[4]);			
		myLine1 = myLine0.replace('\n', ' ');
		myLine = myLine1.replace(',', '.');
		myLine = myLine.replace(' ', '0');
		tabPos2 = findTabPositions(myLine);				
		wallPositions[3][0] = Integer.parseInt(myLine.substring(0, tabPos2[0]));
		for (int k = 1 ; k < jCO.anglesChecked; k++ ) wallPositions[3][k] = Integer.parseInt(myLine.substring(tabPos2[k-1] + 1, tabPos2[k]));			

		return wallPositions;
		
	}
	
	public int[] findLineBreaks(String myString, int numberOfLineBreaks) {
		
		int[] lineBreak = new int[numberOfLineBreaks];
		
		int cc = 1;		
		lineBreak[0] = 0;
		for (int i = 0 ; i < myString.length() ; i++) {
			if (myString.charAt(i) == '\n') {
				lineBreak[cc] = i;
				cc++;
			}
			if (cc == numberOfLineBreaks) break;
		}
		
		return lineBreak;
	}
	
	public boolean writeStringIntoAsciiFile(String path, String myString) {
	       
		try{
            //open file
			FileOutputStream fos = new FileOutputStream(path);
            Writer w = new BufferedWriter(new OutputStreamWriter(fos)); 
           
    		w.write(myString);
    		w.flush();
            w.close();  
            
        }catch(Exception e){
        
        	return true;
        	
        }
        
		return false;
        
    }
	
	public boolean writeStringListIntoAsciiFile(String path, String[] myNames, int[][] myBulkVolumes) {
		
		try{
            
			//open file
			FileOutputStream fos = new FileOutputStream(path);
            Writer w = new BufferedWriter(new OutputStreamWriter(fos));
            
            //write header
        	String myHeader = "FileName\tBulkVolume (vx)\tVolumeAboveTheColumn (vx)\tVolumeBelowTheColumn (vx)\n";
            w.write(myHeader);
                        
            //write the data
            for (int i = 0 ; i < myNames.length ; i++) {
            	
            	String outString = myNames[i] + "\t" + String.format("%d", myBulkVolumes[i][0]) + "\t";
            	outString += String.format("%d", myBulkVolumes[i][1]) + "\t" + String.format("%d", myBulkVolumes[i][2]) + "\n";
            	           	
            	//replace comma by period
            	String cOutString = outString.replace(',', '.');            	
            	
            	//write n flush
            	w.write(cOutString);
            	w.flush();
            }            
            
            //close file
            w.close();
            
        }catch(Exception e){
        
        	return true;
        	
        }
        
		return false;
		
	}
	
	public MyMemory readBrainFile(String pathOfGauges) {
		
		MyMemory myMem = new MyMemory();
		
		try{
            //open file
			FileReader fr = new FileReader(pathOfGauges);		
            BufferedReader br = new BufferedReader(fr);
           
            StringBuilder sb = new StringBuilder();
            String line = "nix";
            int cc = 0;
            
            //read data            
            while (line != null) {
            	
        		sb.append(line + "\n");  
        		
            	line = br.readLine();            
            	
            	if (line == null) return myMem;
            	
            	if ((cc == 0) && (line != null)) myMem.gaugeFolder = line;
            	if ((cc == 1) && (line != null)) myMem.cutAwayFromTop = Integer.parseInt(line);
            	if ((cc == 2) && (line != null)) myMem.filterTag = line;
            	
            	cc++;
            }
            
            myMem.filterImages = true;
            if (cc == 2) myMem.filterImages = false;
            
            br.close();
            
        }catch(Exception e){
        
        	return myMem;        	
        }
		
		return myMem;
	}
	
	public class MyMemory{
		
		public String gaugeFolder;
		public Boolean filterImages;
		public String filterTag;
		public int cutAwayFromTop;
		
	}
	
	public void writeThreshold(String thresholdSaverPath, String myFileName, int myThresh) {
		
		try{
			
            //open file
			BufferedWriter w = new BufferedWriter(new FileWriter(thresholdSaverPath,true)) ;
           
    		w.write(myFileName + "\t" + String.format("%d", myThresh) + "\n");
    		w.flush();
            w.close();  
            
        }
		
		catch(Exception e){}
		
	}
	
	public void writeSnapshots4Comparison(MyFileCollection mFC, ImagePlus nowTiff, ImagePlus binTiff, PolygonRoi[] pRoi, int[] myZ) {
		
		String pathSep = "/";
		
		ContrastEnhancer jCE = new ContrastEnhancer(); 
		
		//create the root output folder for segmentation result comparisons	 
		String myCompDir = mFC.myOutFolder + pathSep + "4Comp" ;		
		new File(myCompDir).mkdir();	

		//find correct filename
		String myFileName = mFC.colName;
		
		//init plotting of hor. cross-sections ..
		ImageStack rgbStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImageStack outStack = new ImageStack(binTiff.getWidth(), binTiff.getHeight());
		
		//check if there is the original image of just some samples
		int[] checkedZ = new int[myZ.length];
		if (nowTiff.getNSlices() == myZ.length) for (int i = 0 ; i < checkedZ.length ; i++) checkedZ[i] = i + 1;
		else checkedZ = myZ;
		
		//really do it..
		for (int i = 0 ; i < checkedZ.length ; i++) {
			
			//set original image and segmented image to the same slice
			nowTiff.setPosition(checkedZ[i]);
			ImageProcessor myIP = nowTiff.getProcessor();
			
			binTiff.setPosition(checkedZ[i]);
			ImageProcessor binIP = binTiff.getProcessor();
			
			//create overlay from binary image
			Overlay myO = new Overlay();
			ThresholdToSelection t2S = new ThresholdToSelection();
			binIP.setBinaryThreshold();
			Roi binRoi = t2S.convert(binIP);			
			myO.add(binRoi);
			
			//create two copies of this IP .. one with and one without the threshold superimposed
			jCE.stretchHistogram(myIP, 0.4);
			ImageProcessor rgbIP = myIP.convertToRGB();
			rgbIP.setColor(Color.YELLOW);			
			
			//ImagePlus test = new ImagePlus("",rgbIP);
			//test.updateAndDraw();
			//test.show();
			
			binRoi.drawPixels(rgbIP);  //superimpose bin Overlay to image
			
			ImageProcessor outIP = myIP.duplicate();
			outIP.setColor(Color.WHITE);
			
			//add label
			String myLabel = "slice " + myZ[i] + "/" + myZ[myZ.length - 1];			
			Font font = new Font("Verdana", Font.PLAIN, 40);
			TextRoi tRoi = new TextRoi((int)(0.01f * binTiff.getWidth()), (int)(0.01f * binTiff.getHeight()),myLabel, font);
			tRoi.drawPixels(rgbIP);
			tRoi.drawPixels(outIP);
			
			
			//ImagePlus sTiff = new ImagePlus("", rgbIP);
			//sTiff.setOverlay(myO);
			//sTiff.updateAndDraw();
			//sTiff.show();
			
			//prep Tiffs		
			rgbStack.addSlice(rgbIP);
			outStack.addSlice(outIP);
		}
			
		//save it
		ImagePlus rgbTiff = new ImagePlus();
		rgbTiff.setStack(rgbStack);		
		jCE.equalize(rgbTiff);		
		rgbTiff.updateAndDraw();		
		String rgbName = myFileName + "Threshed.tif";
		tiffSaver(myCompDir, rgbName, rgbTiff);	
		
		//save it
		ImagePlus outTiff = new ImagePlus();
		outTiff.setStack(outStack);
		jCE.equalize(outTiff);		
		outTiff.updateAndDraw();	
		String outName = myFileName + "Original.tif";		
		tiffSaver(myCompDir, outName , outTiff);		
	
		//combine the two, tune to correct contrast and convert to RGB
/*		StackCombiner jSC = new StackCombiner(); 
		ImageStack combinedImage = jSC.combineHorizontally(oriStack, binStack);
		ImagePlus cI = new ImagePlus();
		cI.setStack(combinedImage);*/
		
		//cI.draw();cI.show();
		
		/*for (int i = 0 ; i < combinedImage.getSize() ; i++) {	
			cI.setPosition(i + 1);
			ImageProcessor toBeSaved = cI.getProcessor().convertToRGB();			
			ImagePlus outout = new ImagePlus("Depth " + i,toBeSaved);
			//outout.draw();outout.show();
			singleTiffSaver(myCompDir+pathSep + eDC[i] + "%_Depth", myFileName, outout);
		}*/
	}
	
	public void writeSnapshots4ComparisonMega(String myOutPath, String myFileName, ImagePlus nowTiff, ImagePlus[] segmentedTiff, 
			double illCF, int largestVal, int lowestVal, float valSpan, PolygonRoi[] pRoi) {
		
		String pathSep = "/";

		ImageProcessor blankIP = new ByteProcessor(nowTiff.getWidth(), nowTiff.getHeight());
		
		//create the root output folder for segmentation result comparisons	 
		String myCompDir = myOutPath + pathSep + "4Comp" ;		
		new File(myCompDir).mkdir();	
		
		//create folder for each horizontal cross-section for checking
		int[] eDC = {20, 40, 60, 80};  //eDC stands for evaluation depth choices 
		for (int i = 0 ; i < eDC.length ; i++) new File(myCompDir+pathSep + eDC[i] + "%_Depth").mkdir();
		
		//define evaluation depths
		int[] myZ = new int[eDC.length];
		for (int i = 0 ; i < eDC.length ; i++) {
			double checkZ = (double)eDC[i] * nowTiff.getNSlices() / 100;
			myZ[i] = (int)Math.round(checkZ);
		}
		
		//init plotting of hor. cross-sections ..
		ImageStack oriStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImageStack airStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImageStack waterStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImageStack stoneStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		
		//really do it..
		for (int i = 0 ; i < myZ.length ; i++) {
			
			//prep orig Tiff
			double stretcher = 2.25;
			nowTiff.setPosition(myZ[i]);
			ImageProcessor myIP = nowTiff.getProcessor();
			myIP.multiply(illCF);
			myIP.subtract(lowestVal);
			myIP.max(stretcher * valSpan);			
			myIP.multiply(256 / (stretcher * valSpan));
			myIP.min(0);
			myIP.max(255);
			ImageProcessor eightIP = myIP.convertToByte(false);	
			if (pRoi != null) {
				eightIP.setRoi(pRoi[myZ[i] - 1]);			
				eightIP.setColor(0);
				eightIP.fillOutside(pRoi[myZ[i] - 1]);
			}
			ContrastEnhancer jCE = new ContrastEnhancer(); 
			jCE.equalize(eightIP);
			oriStack.addSlice(eightIP);
			
			//prep segmented Tiffs
			ImageProcessor outIP = blankIP;
			if (segmentedTiff[0] != null) {
				segmentedTiff[0].setPosition(myZ[i]);			
				outIP = segmentedTiff[0].getProcessor();
			}
			airStack.addSlice(outIP);
			
			//prep segmented Tiffs
			outIP = blankIP;
			if (segmentedTiff[1] != null) {
				segmentedTiff[1].setPosition(myZ[i]);
				outIP = segmentedTiff[1].getProcessor();
			}				
			waterStack.addSlice(outIP);
			
			
			//prep segmented Tiffs
			outIP = blankIP;
			if (segmentedTiff[2] != null) {
				segmentedTiff[2].setPosition(myZ[i]);			
				outIP = segmentedTiff[2].getProcessor();
			}
			stoneStack.addSlice(outIP);
			
		}
			
		//combine the two, tune to correct contrast and convert to RGB
		StackCombiner jSC = new StackCombiner(); 
		ImageStack combinedImage1 = jSC.combineHorizontally(oriStack, airStack);
		ImageStack combinedImage2 = jSC.combineHorizontally(waterStack, stoneStack);		
		
		//cI.draw();cI.show();
		
		for (int i = 0 ; i < combinedImage1.getSize() ; i++) {	
			ImageProcessor cI1 = combinedImage1.getProcessor(i + 1);
			ImageProcessor cI2 = combinedImage2.getProcessor(i + 1);
			
			ImageProcessor combo = new ByteProcessor(combinedImage1.getWidth(), 2 * combinedImage1.getHeight());
			for (int x = 0 ; x < combo.getWidth() ; x++) {
				for (int y = 0 ; y < combo.getHeight() ; y++) {
					if (y < combinedImage1.getHeight()) combo.putPixelValue(x, y, cI1.getPixelValue(x, y));
					else combo.putPixelValue(x, y, cI2.getPixelValue(x, y - combinedImage1.getHeight()));
				}
			}

			ImagePlus outout = new ImagePlus("Depth " + i, combo.convertToRGB());
			
			
			save2DTiff(myCompDir+ pathSep + eDC[i] + "%_Depth", myFileName + ".tif", outout);
		}
    }
	
	public MyFileCollection createFolders4SubROIData(MenuWaiter.ROISelectionOptions mRSO) {
		
		String pathSep = "/";		
		
		//select file or files
	  	MyFileCollection mFC = fileSelector("Please choose a file or folder with your image data");	
					
		//create output paths
		String myPreOutFolder = "";		
		if (mRSO.choiceOfRoi.equals("RealSample")) {			
			
			if (mRSO.heightOfRoi > 0) myPreOutFolder = "Top" + mRSO.cutAwayFromTop + "Height" + mRSO.heightOfRoi + "Wall" + mRSO.cutAwayFromWall;
			boolean isCut = false;
			if (mRSO.cutAwayFromBottom > 0 | mRSO.cutAwayFromTop > 0 | mRSO.cutAwayFromWall > 0 | mRSO.cutAwayFromCenter > 0) isCut = true; 
			if (isCut & mRSO.heightOfRoi > 0) myPreOutFolder = "Tv" + mRSO.cutAwayFromTop + "Height" + mRSO.heightOfRoi + "Wall" + mRSO.cutAwayFromWall;
			if (isCut & mRSO.cutZPercent & mRSO.cutXYPercent & mRSO.heightOfRoi == 0) myPreOutFolder = 
					"Tp" + mRSO.cutAwayFromTop + "Bp" + mRSO.cutAwayFromBottom +	"Wp" + mRSO.cutAwayFromWall + "Cp" + mRSO.cutAwayFromCenter;
			if (isCut & !mRSO.cutZPercent & mRSO.cutXYPercent & mRSO.heightOfRoi == 0) myPreOutFolder = 
					"Tv" + mRSO.cutAwayFromTop + "Bv" + mRSO.cutAwayFromBottom +	"Wp" + mRSO.cutAwayFromWall + "Cp" + mRSO.cutAwayFromCenter;
			if (isCut & mRSO.cutZPercent & !mRSO.cutXYPercent & mRSO.heightOfRoi == 0) myPreOutFolder = 
					"Tp" + mRSO.cutAwayFromTop + "Bp" + mRSO.cutAwayFromBottom +	"Wv" + mRSO.cutAwayFromWall + "Cv" + mRSO.cutAwayFromCenter;
			if (isCut & !mRSO.cutZPercent & !mRSO.cutXYPercent & mRSO.heightOfRoi == 0) myPreOutFolder = 
					"Tv" + mRSO.cutAwayFromTop + "Bv" + mRSO.cutAwayFromBottom +	"Wv" + mRSO.cutAwayFromWall + "Cv" + mRSO.cutAwayFromCenter;
			if (!isCut) myPreOutFolder = "InnerCircleColumn";	

			if (mRSO.includeSurfaceTopography) myPreOutFolder = "S_"+ myPreOutFolder;
			else myPreOutFolder = "C_"+ myPreOutFolder;
			mFC = getAllMyNeededFolders(mFC.myBaseFolder, mFC.myTiffs, "", mRSO.useInnerCircleFiles, mRSO.includeSurfaceTopography);   //"" put because no possibility for a filterTag implemented yet
			
		} else {
						
			if (mRSO.choiceOfRoi.equals("Cuboid")) myPreOutFolder = "Cuboid_XL" + mRSO.cubeX1 + "XR" + mRSO.cubeX2 + "YL" + mRSO.cubeY1 + "YR" + mRSO.cubeY2 + "ZT" + mRSO.cubeZ1 + "ZB" + mRSO.cubeZ2;
			if (mRSO.choiceOfRoi.equals("Cylinder")) myPreOutFolder = "Cyl_X" + mRSO.cylX + "Y" + mRSO.cylY+ "ZT" + mRSO.cylZ1 + "ZB" + mRSO.cylZ2 + "R" + mRSO.cylRadius;
			if (mRSO.choiceOfRoi.equals("Everything!")) myPreOutFolder = "EntireImage";
			if (mRSO.choiceOfRoi.equals("TopOfEveryThing")) myPreOutFolder = "TopOfEntireImage";
			if (mRSO.choiceOfRoi.equals("BottomOfEveryThing")) myPreOutFolder = "BottomOfEntireImage";
			
		}
		
		String myPreOutPath = mFC.myBaseFolder + pathSep + "Histo_" + myPreOutFolder;
		new File(myPreOutPath).mkdir();		
		mFC.myOutFolder = myPreOutPath; 	//also add it to the folder collection
		
		//save pathSep
		mFC.pathSep = pathSep;
		
		return mFC;
		
	}
	
	public MyFileCollection createFolders4SubROIData(MenuWaiter.PoreSpaceAnalyzerOptions mPSA) {
		
		String pathSep = "/";	
		
		MyFileCollection mFC = fileSelector("Please choose a file or folder with your image data");
					
		//create output paths
		String myPreOutFolder = "";		
		if (mPSA.mRSO.choiceOfRoi.equals("RealSample")) {			
			
			if (mPSA.mRSO.heightOfRoi > 0) myPreOutFolder = "Top" + mPSA.mRSO.cutAwayFromTop + "Height" + mPSA.mRSO.heightOfRoi + "Wall" + mPSA.mRSO.cutAwayFromWall;
			boolean isCut = false;
			if (mPSA.mRSO.cutAwayFromBottom > 0 | mPSA.mRSO.cutAwayFromTop > 0 | mPSA.mRSO.cutAwayFromWall > 0 | mPSA.mRSO.cutAwayFromCenter > 0) isCut = true; 
			if (isCut & mPSA.mRSO.heightOfRoi > 0) myPreOutFolder = "Tv" + mPSA.mRSO.cutAwayFromTop + "Height" + mPSA.mRSO.heightOfRoi + "Wall" + mPSA.mRSO.cutAwayFromWall;
			if (isCut & mPSA.mRSO.cutZPercent & mPSA.mRSO.cutXYPercent & mPSA.mRSO.heightOfRoi == 0) myPreOutFolder = 
					"Tp" + mPSA.mRSO.cutAwayFromTop + "Bp" + mPSA.mRSO.cutAwayFromBottom +	"Wp" + mPSA.mRSO.cutAwayFromWall + "Cp" + mPSA.mRSO.cutAwayFromCenter;
			if (isCut & !mPSA.mRSO.cutZPercent & mPSA.mRSO.cutXYPercent & mPSA.mRSO.heightOfRoi == 0) myPreOutFolder = 
					"Tv" + mPSA.mRSO.cutAwayFromTop + "Bv" + mPSA.mRSO.cutAwayFromBottom +	"Wp" + mPSA.mRSO.cutAwayFromWall + "Cp" + mPSA.mRSO.cutAwayFromCenter;
			if (isCut & mPSA.mRSO.cutZPercent & !mPSA.mRSO.cutXYPercent & mPSA.mRSO.heightOfRoi == 0) myPreOutFolder = 
					"Tp" + mPSA.mRSO.cutAwayFromTop + "Bp" + mPSA.mRSO.cutAwayFromBottom +	"Wv" + mPSA.mRSO.cutAwayFromWall + "Cv" + mPSA.mRSO.cutAwayFromCenter;
			if (isCut & !mPSA.mRSO.cutZPercent & !mPSA.mRSO.cutXYPercent & mPSA.mRSO.heightOfRoi == 0) myPreOutFolder = 
					"Tv" + mPSA.mRSO.cutAwayFromTop + "Bv" + mPSA.mRSO.cutAwayFromBottom +	"Wv" + mPSA.mRSO.cutAwayFromWall + "Cv" + mPSA.mRSO.cutAwayFromCenter;
			if (!isCut) myPreOutFolder = "InnerCircleColumn";	
						
			if (mPSA.mRSO.includeSurfaceTopography) myPreOutFolder = "S_"+ myPreOutFolder;
			else myPreOutFolder = "C_"+ myPreOutFolder;
			mFC = getAllMyNeededFolders(mFC.myBaseFolder, mFC.myTiffs, "", mPSA.mRSO.useInnerCircleFiles, mPSA.mRSO.includeSurfaceTopography);   //"" put because no possibility for a filterTag implemented yet
			
		} else {
			
			if (mPSA.mRSO.choiceOfRoi.equals("Cuboid")) myPreOutFolder = "Cuboid_XL" + mPSA.mRSO.cubeX1 + "XR" + mPSA.mRSO.cubeX2 + "YL" + mPSA.mRSO.cubeY1 + "YR" + mPSA.mRSO.cubeY2 + "ZT" + mPSA.mRSO.cubeZ1 + "ZB" + mPSA.mRSO.cubeZ2;
			if (mPSA.mRSO.choiceOfRoi.equals("Cylinder")) myPreOutFolder = "Cyl_X" + mPSA.mRSO.cylX + "Y" + mPSA.mRSO.cylY+ "ZT" + mPSA.mRSO.cylZ1 + "ZB" + mPSA.mRSO.cylZ2 + "R" + mPSA.mRSO.cylRadius;
			if (mPSA.mRSO.choiceOfRoi.equals("Everything!")) myPreOutFolder = "EntireImage";
			if (mPSA.mRSO.choiceOfRoi.equals("TopOfEveryThing")) myPreOutFolder = "TopOfEntireImage";
			if (mPSA.mRSO.choiceOfRoi.equals("BottomOfEveryThing")) myPreOutFolder = "BottomOfEntireImage";

		}
		
		String myPreOutPath = mFC.myBaseFolder + pathSep + mPSA.imagePhase2BeAnalyzed + "_" + mPSA.nameOfAnalyzedPhase;
		new File(myPreOutPath).mkdir();	
		myPreOutPath += pathSep + myPreOutFolder;
		new File(myPreOutPath).mkdir();
		mFC.myPreOutFolder = myPreOutPath; 	//also add it to the folder collection
		
		//also create a folder for the cut surfaces
		if (mPSA.mRSO.includeSurfaceTopography) {
			String myCutSurfacePath = mFC.myPreOutFolder + pathSep + "CutSurfaceFiles";
			new File(myCutSurfacePath).mkdir();
			mFC.myCutSurfaceFolder = myCutSurfacePath;
		}
					
		//save pathSep
		mFC.pathSep = pathSep;
		
		return mFC;
		
	}
	
	
	public int[] findStartAndStopSlices(MyFileCollection mFC, MenuWaiter.PoreSpaceAnalyzerOptions mPSA) {
		
		ObjectDetector jOD = new ObjectDetector();
		
		//check if soil surface should be taken into account
		int topSurface = 0;
		int botSurface = mFC.nOfSlices;
		if (mPSA.mRSO.includeSurfaceTopography) {			
			MorphologyAnalyzer.SurfaceStatistics mST = jOD.extractSurfaceStatistics(mFC);
			topSurface = mST.medianElevation;
			botSurface = mST.medianIntrusion;
		}
		
		//define start and stop
		int startSlice = topSurface + mPSA.mRSO.cutAwayFromTop;
		int stopSlice = botSurface - mPSA.mRSO.cutAwayFromBottom;
		
		//check if there is something to cut away
		if (mPSA.mRSO.choiceOfRoi.equalsIgnoreCase("Cuboid")) mPSA.mRSO.cutAwayFromTop = mPSA.mRSO.cubeZ1 - 1; 
		
		if (mPSA.mRSO.choiceOfRoi.equalsIgnoreCase("Cylinder")) {			
			mPSA.mRSO.cutAwayFromTop = mPSA.mRSO.cylZ1; 
			startSlice = mPSA.mRSO.cylZ1;
			stopSlice = mPSA.mRSO.cylZ2;
			if (stopSlice == 0) stopSlice = mFC.nOfSlices;
		}
		if (mPSA.mRSO.cutAwayFromTop < 0) mPSA.mRSO.cutAwayFromTop = 0;
		if (startSlice == 0) startSlice = 1; 			
		
		if (mPSA.mRSO.cutZPercent) {
			startSlice = (int)Math.round((float)mPSA.mRSO.cutAwayFromTop / 100f * (float)mFC.nOfSlices);
			if (startSlice == 0) startSlice = 1; 
			if (mPSA.mRSO.cutAwayFromBottom == 0) stopSlice = mFC.nOfSlices;
			else stopSlice = mFC.nOfSlices - (int)Math.round((float)mPSA.mRSO.cutAwayFromBottom / 100f * (float)mFC.nOfSlices);
		}
		if (mPSA.mRSO.heightOfRoi > 0) stopSlice = startSlice + mPSA.mRSO.heightOfRoi;
				
		int[] startStop = {startSlice, stopSlice};
		
		return startStop;
		
	}
	
	public int[] findStartAndStopSlices(MyFileCollection mFC, MenuWaiter.ROISelectionOptions mRSO) {
		
		ObjectDetector jOD = new ObjectDetector();
		
		//check if soil surface should be taken into account
		int topSurface = 1;
		int botSurface = mFC.nOfSlices;
		if (mRSO.includeSurfaceTopography) {			
			MorphologyAnalyzer.SurfaceStatistics mST = jOD.extractSurfaceStatistics(mFC);
			topSurface = mST.medianElevation;
			botSurface = mST.medianIntrusion;
		}
		
		//check if there is something to cut away
		if (mRSO.choiceOfRoi.equalsIgnoreCase("Cuboid")) mRSO.cutAwayFromTop = mRSO.cubeZ1 - 1; 
		if (mRSO.choiceOfRoi.equalsIgnoreCase("Cylinder")) mRSO.cutAwayFromTop = mRSO.cylZ1 - 1; 
		if (mRSO.cutAwayFromTop < 0) mRSO.cutAwayFromTop = 0;
		int startSlice = topSurface + mRSO.cutAwayFromTop;
		mFC.startSlice = startSlice;
		int stopSlice = botSurface - mRSO.cutAwayFromBottom;
		if (mRSO.cutZPercent) {
			startSlice = (int)Math.round((float)mRSO.cutAwayFromTop / 100f * (float)mFC.nOfSlices);
			stopSlice = mFC.nOfSlices - (int)Math.round((float)mRSO.cutAwayFromBottom / 100f * (float)mFC.nOfSlices);
		}
		if (mRSO.heightOfRoi > 0) stopSlice = startSlice + mRSO.heightOfRoi;
				
		int[] startStop = {startSlice, stopSlice};
		
		return startStop;
		
	}
	
	public void writeIlluminationCorrectionIntoAsciiFile(MenuWaiter.CalibrationReferences myNR, String outDir, String fileName) {
		
		//Set
		String myFilePath = outDir + "//" + "Illumi_" + fileName + ".txt";		
	
		try{
			
            //open file
			BufferedWriter w = new BufferedWriter(new FileWriter(myFilePath,false)) ;
           	
			String targetValues = String.format("%5.0f\t", myNR.lowerTarget) + "\t" + String.format("%5.0f\t", myNR.upperTarget) + "\n\n";
			w.append(targetValues);
			w.flush();
			
			for (int i = 0 ; i < myNR.originalLower.length ; i++) {
				String outString = String.format("%5.0f\t", myNR.originalLower[i]) + "\t" + String.format("%5.0f\t", myNR.originalUpper[i]) + "\n";
				w.append(outString);
				w.flush();
			}			
            w.close();
        }	
		
		catch(Exception e){}
		
	}
	
	public void writeHistogram16(MyFileCollection mFC, int[][] myHist) {
	
		String savePath = mFC.myOutFolder + mFC.pathSep + "Histograms.16b"; 
		
		try {
			
            //open file
			BufferedWriter w = null;			
			w = new BufferedWriter(new FileWriter(savePath,false));
		
			for (int line = 0 ; line < myHist.length ; line++) {
				
				String colName = mFC.myTiffs[line].substring(0,mFC.myTiffs[line].length()-4);
				String nowLine = colName + "\t";
				
				IJ.showStatus("Saving histogram of file " + (line + 1) + "/" + mFC.myTiffs.length + " into histogram folder..");
				
				for (int column = 0 ; column < (int)Math.round(Math.pow(2, 16)) - 1 ; column++) {
					
					nowLine += Integer.toString(myHist[line][column]) + "\t";
					
				}
				
				nowLine += Integer.toString(myHist[line][(int)Math.round(Math.pow(2, 16)) - 1 ]) + "\n";
				
				w.write(nowLine);
				w.flush();
				
			}
			
			w.close();
            
        }
		
		catch(Exception e){
			IJ.error("Exception " + e.toString() + "\nSaving path: " + savePath);
		}
	    
	}
	
	public void writeHistogram8(MyFileCollection mFC, int[][] myHist) {
		
		String savePath = mFC.myOutFolder + mFC.pathSep + "Histograms.8b"; 
		
		try {
			
            //open file
			BufferedWriter w = null;			
			w = new BufferedWriter(new FileWriter(savePath,false));
		
			for (int line = 0 ; line < myHist.length ; line++) {
				
				String colName = mFC.myTiffs[line].substring(0,mFC.myTiffs[line].length()-4);
				String nowLine = colName + "\t";
				
				IJ.showStatus("Saving histogram of file " + (line + 1) + "/" + mFC.myTiffs.length);
				
				for (int column = 0 ; column < (int)Math.round(Math.pow(2, 8)) - 1 ; column++) {
					
					nowLine += Integer.toString(myHist[line][column]) + "\t";
					
				}
				
				nowLine += Integer.toString(myHist[line][(int)Math.round(Math.pow(2, 8)) - 1 ]) + "\n";
				
				w.write(nowLine);
				w.flush();
				
			}
			
			w.close();
            
        }
		
		catch(Exception e){
			IJ.error("Exception " + e.toString() + "\nSaving path: " + savePath);
		}
	    
	}
	
	public void writeMultiThresholds(MyFileCollection mFC, int[][] myThresholds, MenuWaiter.MultiRegionSegmentation mRS) {
		
		String savePath = mFC.myOutFolder + mFC.pathSep + mFC.colName;
		
		//re-calculate to 16 bit if necessary
		if (mFC.bitDepth == 16) {
			for (int i = 0 ; i < myThresholds.length ; i++) { 
				for (int j = 0 ; j < myThresholds[i].length ; j++) {
					if (mRS.lowerBoundary > 0 | mRS.upperBoundary > 0) {
						myThresholds[i][j] = (int)Math.round((float)myThresholds[i][j] / 256 * (mRS.upperBoundary - mRS.lowerBoundary) + mRS.lowerBoundary);
					}
					else myThresholds[i][j] = myThresholds[i][j] * 256;
				}
			}
		}
		
		try {
			
            //open file
			BufferedWriter w = null;			
			w = new BufferedWriter(new FileWriter(savePath,false));
		
			String line2Write = "Multi Otsu thresholds:       \t";
			for (int i = 0 ; i < mRS.numberOfThresholds - 1 ; i++) line2Write += myThresholds[0][i] + "\t";
			line2Write += myThresholds[0][mRS.numberOfThresholds - 1] + "\n";
			
			w.write(line2Write);
			w.flush();		
		
			line2Write = "Multi max entropy thresholds:\t";
			for (int i = 0 ; i < mRS.numberOfThresholds - 1 ; i++) line2Write += myThresholds[1][i] + "\t";
			line2Write += myThresholds[1][mRS.numberOfThresholds - 1] + "\n";
			
			w.write(line2Write);
			w.flush();		
			
			line2Write = "Multi minimum error:         \t";
			for (int i = 0 ; i < mRS.numberOfThresholds - 1 ; i++) line2Write += myThresholds[2][i] + "\t";
			line2Write += myThresholds[2][mRS.numberOfThresholds - 1] + "\n";
			
			w.write(line2Write);
			w.flush();		
			
			w.close();
        
		}
	
		catch(Exception e){
			IJ.error("Exception " + e.toString() + "\nSaving path: " + savePath);
		}
	}
	
	public void write2DHistogram16(String path, double[][] out2DHist) {
		
		int b8 = 256;
		
		try {
			
			 FileWriter file = new FileWriter(path);
			    
			 file.write("-1\t");
			 for (int x = 0 ; x < b8 - 1 ; x++) file.write((x * 256 + 128) + "\t");
			 file.write((255 * 256 + 128) + "\n");
			 
			 for(int y = 0 ; y < b8 ; y++) {
				 file.write(String.format("%1.4e\t", (double)(y) * 25.6d + 128d/25.6d));		
				 for (int x = 0 ; x <  b8 - 1; x++) {				 
				      file.write(String.format("%1.4e\t", out2DHist[x][y]));			 
				 }
				 file.write(String.format("%1.4e\n", out2DHist[255][y]));	
			 }
			 
			 file.close();			
		}

		catch (IOException e) {
		      System.out.println("Error - " + e.toString());
	    }

	}
	
	public ImagePlus loadCarvedTiff(MyFileCollection mFC, MenuWaiter.ROISelectionOptions mRSO) {
	
		ObjectDetector jOD = new ObjectDetector();
		
		//check if soil surface should be taken into account
		int topSurface = 0;
		int botSurface = mFC.nOfSlices;
		if (mRSO.includeSurfaceTopography) {			
			MorphologyAnalyzer.SurfaceStatistics mST = jOD.extractSurfaceStatistics(mFC);
			topSurface = mST.medianElevation;
			botSurface = mST.medianIntrusion;
		}
		
		//check if there is something to cut away
		int startSlice = topSurface + mRSO.cutAwayFromTop;
		mFC.startSlice = startSlice;
		int stopSlice = botSurface - mRSO.cutAwayFromBottom;
		if (mRSO.cutZPercent) {
			startSlice = (int)Math.round((float)mRSO.cutAwayFromTop / 100f * (float)mFC.nOfSlices);
			stopSlice = mFC.nOfSlices - (int)Math.round((float)mRSO.cutAwayFromBottom / 100f * (float)mFC.nOfSlices);
		}
		if (mRSO.heightOfRoi > 0) stopSlice = startSlice + mRSO.heightOfRoi;
				
		//load file
		int[] colSlices = new int[stopSlice - startSlice];
		for (int j = 0 ; j < colSlices.length ; j++) colSlices[j] = startSlice + j;
		ImagePlus nowTiff = openTiff3DSomeSlices(mFC, colSlices);
		
		return nowTiff;
		
	}
	
	public int returnLinesInAsciiFile(String myPath) {
		
		try {		
		
			FileReader fr = new FileReader(myPath);		
	        BufferedReader br = new BufferedReader(fr);
	        
	        String line = "something";
	        int cc = 0;
	        while (line != null) {
	        	line = br.readLine();
	        	cc++;	        
	        }	     
	        int histogramNumber = cc - 1;
	        
	        //move file pointer back to start of the file
	        br.close();
	        fr.close();
	        
	        return histogramNumber;
		}
		catch(Exception e) {
			
			IJ.error("Something went wrong when reading from the histogram file...");
			return 0;
			
        }
	}
	
	public int returnHistogramClassesFromFile(String myPath) {
		
		try {		
		
			FileReader fr = new FileReader(myPath);		
	        BufferedReader br = new BufferedReader(fr);
	        
	        String line = br.readLine();			
			String[] wordsArray = line.split("\t");
			
			int bitDepth = wordsArray.length - 1;
	        
	        //move file pointer back to start of the file
	        br.close();
	        fr.close();
	        
	        return bitDepth;
		}
		catch(Exception e) {
			
			IJ.error("Something went wrong when reading from the histogram file...");
			return 0;
			
        }
	}
	
	public Histograms readHistogram(String myPath) {
		
		Histograms hists = new Histograms();
		
		int histogramNumber = returnLinesInAsciiFile(myPath);
		int histogramClasses = returnHistogramClassesFromFile(myPath);
		
		try {
			
			FileReader fr = new FileReader(myPath);		
	        BufferedReader br = new BufferedReader(fr);
	        
			int[][] myHists = new int[histogramNumber][histogramClasses]; 
			
			//read data
			ArrayList<String> samples = new ArrayList<>();
			
			for (int i = 0 ; i < histogramNumber ; i++) {
				
				String line = br.readLine();			
				String[] wordsArray = line.split("\t");
				ArrayList<String> words = new ArrayList<>(); 
				
				for(String each : wordsArray){
					if(!"".equals(each)){
						words.add(each);
					}					
				}
				
				//only start at line 1 to skip the header
				samples.add(words.get(0));
				for (int j = 1 ; j < histogramClasses ; j++) myHists[i][j] = Integer.parseInt(words.get(j));
				
			}
		
			br.close();
			fr.close();
			
			hists.histogramClasses = histogramClasses;
			hists.numberOfHistograms = histogramNumber;
			hists.histograms = myHists;
			hists.sampleNames = samples;
			
			return hists;
			
		}
		catch(Exception e) {
        			
			IJ.error("Something went wrong when reading from the histogram file...");
			return null;
			
        }
		
	}
	
	public int[][] readHistogram(String nowGaugePath, ObjectDetector.ColCoords3D jCO) {
		
		String pathSep = "/";
		
		int cc = 0;
		int[][] myHist = new int[jCO.heightOfColumn][256 * 256];
		
		//get filename
		int lastBackSlash = 0;
		
		for (int i = 0 ; i < nowGaugePath.length() ; i++) if (nowGaugePath.charAt(i) == pathSep.charAt(0)) lastBackSlash = i; 
		String myFilePath = nowGaugePath.substring(0, lastBackSlash) + pathSep + "Histogram_" + nowGaugePath.substring(lastBackSlash + 7, nowGaugePath.length() - 4) + ".txt";		
		
		String histAsString0 = "";
		
		try{
            //open file
			FileReader fr = new FileReader(myFilePath);		
            BufferedReader br = new BufferedReader(fr);
           
            StringBuilder sb = new StringBuilder();
            String line = "";            
            
            //read data            
            while (line != null) {
            	
        		sb.append(line + "\n");  
            	line = br.readLine();
            	cc++;
            }            
            
            histAsString0 = sb.toString();
         
            br.close();
            
        } catch(Exception e) {return null;}
		
		String histAsString = histAsString0.replace('\r', ' ');
		int[] mLB = findLineBreaks(histAsString, cc + 1);
		for (int i = 1 ; i < mLB.length - 1 ; i++) {
			
			IJ.showStatus("Reading histogram of slice #" + (i + 1) + "/" + cc);
			
			String myLine0 = histAsString.substring(mLB[i], mLB[i + 1]);
			String myLine = myLine0.replace('\n', ' ');			
			
			int[] tabPos = findTabPositions(myLine);
			
			for (int j = 0 ; j < tabPos.length - 1; j++) myHist[i][j] = Integer.parseInt(myLine.substring(tabPos[j] + 1, tabPos[j + 1]));;						
		}
		
		return myHist;		
	}

	public void saveCriticalGreyValues(String myGauge, int[] onePercQuantile, double[] wallOfGrey) {
		
		String pathSep = "/";
		
		//get filename
		int lastBackSlash = 0;
		for (int i = 0 ; i < myGauge.length() ; i++) if (myGauge.charAt(i) == pathSep.charAt(0)) lastBackSlash = i; 
		String myFilePath = myGauge.substring(0, lastBackSlash) + pathSep + "Q001PlusWall_" + myGauge.substring(lastBackSlash + 7, myGauge.length() - 4) + ".txt";		
	
		try{
			
            //open file
			BufferedWriter w = new BufferedWriter(new FileWriter(myFilePath,false)) ;
           			    		
			for (int i = 0 ; i < onePercQuantile.length ; i++) {
				String outString = onePercQuantile[i] + "\t" + (int)Math.round(wallOfGrey[i]) + "\n";
				w.append(outString);
				w.flush();
			}			
            w.close();
        }		
		catch(Exception e){}		
	}
	
	public double[] readSingleColumnDoublesFromAscii(String myPorosityListFile) {

		ArrayList<Double> myList = new ArrayList<Double>();
		
		try{
			
            //open file
			BufferedReader r = new BufferedReader(new FileReader(myPorosityListFile)) ;
           	
			//read first line
			String line = r.readLine();			
				
			//go through all the lines
			while (line != null) {
				myList.add(Double.parseDouble(line));
				line = r.readLine();	
			}			
            r.close();
        }		
		catch(Exception e){}		
		
		//store myList in an array
		double[] myDoubles = new double[myList.size()];
		for (int i = 0 ; i < myDoubles.length ; i++) myDoubles[i] = myList.get(i);
		
		return myDoubles;
	}
		
}
