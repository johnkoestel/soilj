package SoilJ.tools;

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


import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.StatUtils;
import org.doube.bonej.Thickness;
import org.doube.jama.EigenvalueDecomposition;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.PolygonRoi;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import SoilJ.copiedTools.JAnisotropy;
import SoilJ.copiedTools.JFractalCount;
import SoilJ.copiedTools.JParticleCounter;
import SoilJ.tools.InputOutput;
import SoilJ.tools.ObjectDetector;
import SoilJ.tools.ObjectDetector.ColCoords3D;

/** 
 * is a SoilJ class
 * 
 * @author John Koestel
 *
 */

public class MorphologyAnalyzer implements PlugIn {

	public void run(String arg) {
		//ok, this is not needed..
	}

	public class SurfaceStatistics {
		
		public int highestElevation;
		public int medianElevation;
		public int meanElevation;
		public int lowestElevation;
		
		public int highestIntrusion;
		public int medianIntrusion;
		public int meanIntrusion;
		public int lowestIntrusion;		
	}		
	
	public class ProfileStatistics {
		
		public int[] numberOfNonZeros;
		public double[] mean;
		public double[] geomean;	
		public double[] logMean;
		public double[] median;
		public double[] mode;
		public double[] std;
		public double[] mini;
		public double[] maxi;		
		
	}
	
	public void extractPoresizeDistro(int i, InputOutput.MyFolderCollection mFC, ImagePlus nowTiff) {
		
		InputOutput jIO = new InputOutput();
		
		int[] psd = new int[256];
		
		//get histogram
		for (int j = 0 ; j < nowTiff.getNSlices() ; j++) { 
		
			nowTiff.setPosition(j+1);
			ImageProcessor nowIP = nowTiff.getProcessor();
						 
			for (int x = 0 ; x < nowTiff.getWidth() ; x++) {
				for (int y = 0 ; y < nowTiff.getHeight() ; y++) {
					int nowPix = (int)Math.round(nowIP.getPixelValue(x, y));
					if (nowPix > 255) psd[255]++;
					else psd[nowPix]++;
				}
			}				
		}		
		
		//convert histogram into string
		String myHist = "" + psd[0] + "\t"; 
		for (int j = 1 ; j < psd.length ; j++) myHist += psd[j] + "\t";
		
		//save histogram
		jIO.writeStringIntoAsciiFile(mFC.myPreOutFolder + "\\" + mFC.myTiffs[i].substring(0, mFC.myTiffs[i].length() - 4), myHist);
	}
	
	public ImagePlus extractPercolatingPorosity(ImagePlus idImp, ArrayList<Integer> percolatingClusters) {
		
		ImageStack outStack = new ImageStack(idImp.getWidth(), idImp.getHeight());
				
		for (int z = 1 ; z <= idImp.getNSlices() ; z++) {
			
			IJ.showStatus("Extracting percolating clusters in slice " + z + "/" + idImp.getNSlices());
		
			ImageProcessor binIP = new ByteProcessor(idImp.getWidth(), idImp.getHeight());
			idImp.setPosition(z);
			ImageProcessor idIP = idImp.getProcessor();
			
			for (int x = 0 ; x < binIP.getWidth() ; x++) {
				
				for (int y = 0 ; y < binIP.getWidth() ; y++) {
					
					int id = (int)Math.round(idIP.getPixelValue(x, y));
					if (percolatingClusters.contains(id)) binIP.putPixel(x, y, 255);
					else binIP.putPixel(x, y, 0);
					
				}
				
			}
			
			outStack.addSlice(binIP);
			
		}
		
		ImagePlus outTiff = new ImagePlus();
		outTiff.setStack(outStack);
		
		return outTiff;
		
	}
	
	public PoreClusterProperties calculateChiOfPercolation(ImagePlus nowTiff, InputOutput.MyFolderCollection mFC, PoreClusterProperties mCP) {
			
			AndAllTheRest aa = new AndAllTheRest(); 
		
			int i,j;
			int x,y,z;
			
			double dSum = 0;
			double gSumSum = 0;
			
			double[] connectedCorrelationLength = new double[mCP.id.length];
			double numerator = 0;
			double gSum = 0;
			double[] chi = new double[mCP.id.length];
			
			IJ.showStatus("Calculating Chi of percolation ...");
			IJ.freeMemory();IJ.freeMemory();
			
			//sweep through the clusters
			int numberOfConsideredClusters = mCP.id.length;  //all clusters would be  mCP.id.length			
			for (i = 0 ; i < numberOfConsideredClusters ; i++) {			 
				
				int nowCluster = mCP.id[i];
					
				if (mCP.volume[i] == 1) {
					
					connectedCorrelationLength[i] = 0;
					numerator = 0;
					gSum = 0;
					
				} else if (!mCP.isPercolating[i]) {
					
					//long startTime = System.currentTimeMillis();		
					
					double rmax = 0;		//maximally possible distance between two voxels		
					
					//find all voxels belonging to the cluster	
					List<Integer> nX = new ArrayList<Integer>();
					List<Integer> nY = new ArrayList<Integer>();
					List<Integer> nZ = new ArrayList<Integer>();
					
					//define bounding box around the cluster 
					int x0 = mCP.minX[i];
					int y0 = mCP.minY[i];
					int z0 = mCP.minZ[i];
					int xe = mCP.maxX[i];
					int ye = mCP.maxY[i];
					int ze = mCP.maxZ[i];
					
					double xy = Math.sqrt((xe - x0) * (xe - x0) + (ye - y0) * (ye - y0));
					rmax = Math.sqrt((ze - z0) * (ze - z0) + xy * xy);
											
					for (z = z0 ; z <= ze ; z++) { 
						nowTiff.setPosition(z + 1);
						ImageProcessor nowIP = nowTiff.getProcessor();
						
						for (x = x0 ; x <= xe ; x++) {
							for (y = y0 ; y <= ye ; y++) {
								int nowPixel = (int)nowIP.getPixelValue(x, y);
								if (nowPixel == nowCluster) { 
									nX.add(x);
									nY.add(y);
									nZ.add(z);
								}
							}
						}					
					}				
							
					//give a warning if the cluster size is too large
					if (nX.size() > Integer.MAX_VALUE) IJ.error("cluster size exceeds Java Integer range.. arrgghhh...");
					
					//rewrite listed coordinates as array
					int nowX[] = new int[nX.size()];
					int nowY[] = new int[nX.size()];
					int nowZ[] = new int[nX.size()];
					for(j = 0 ; j < nX.size() ; j++) {
						nowX[j] = nX.get(j);
						nowY[j] = nY.get(j);
						nowZ[j] = nZ.get(j);
					}
					
					//if cluster size is larger than w * rmax take a random subset of all voxels below this number
					int numOfConsideredVoxels = nowX.length;
					int requiredMultipleOfRMax = 100;
					ArrayList<Integer> retain = new ArrayList<Integer>();
					if (numOfConsideredVoxels > requiredMultipleOfRMax * rmax) {
						Random randomGenerator = new Random();					
						while (retain.size() < (int)Math.floor(requiredMultipleOfRMax * rmax) + 1) {						
							int newNum = randomGenerator.nextInt(nowX.length);
							if (!retain.contains(newNum)) retain.add(newNum);
						}
						numOfConsideredVoxels = (int)Math.floor(requiredMultipleOfRMax * rmax);										
					}
					
					//calculate distances	
					int iRmax = (int)Math.ceil(rmax) + 1;
					int[] myHist = new int[iRmax];
					for (j = 0 ; j < myHist.length ; j++) myHist[j] = 0;						
					for (int a = 0 ; a < numOfConsideredVoxels ; a++) {
						for (int b = 0 ; b < numOfConsideredVoxels ; b++) {
							int ia = a; int ib = b;
							if (!retain.isEmpty()) {
								ia = retain.get(a);
								ib = retain.get(b);
							}
							int dx = nowX[ia] - nowX[ib];
							int dy = nowY[ia] - nowY[ib];
							int dz = nowZ[ia] - nowZ[ib];	
							double dxy = Math.sqrt(dx * dx + dy *dy);
							double dist = Math.sqrt(dxy * dxy + dz * dz);
							myHist[(int)Math.round(dist)]++;											
						}					
					}
									
					//normalize everything
					double[] r = new double[iRmax];
					for (j = 0 ; j < iRmax ; j++) r[j] = (j + 0.5); 
					double[] g = new double[iRmax];
					for (j = 0 ; j < iRmax; j++) g[j] = myHist[j] / (mCP.volume[i] - 1);
					
					//calculate pseudo chi
					double numer = 0;
					double denom = 0;
					for (j = 0 ; j < myHist.length ; j++) numer += r[j] * r[j] * g[j];
					for (j = 0 ; j < myHist.length ; j++) denom += g[j];
					connectedCorrelationLength[i] = Math.sqrt(numer / denom);	
					numerator = numer;
					gSum = denom;
											
	/*				float interval = (System.currentTimeMillis() - startTime) / 1000;		
					DecimalFormat df = new DecimalFormat();
					df.setMaximumFractionDigits(2);		
					IJ.error("Calculation time for " + numOfConsideredVoxels + " voxels was " + df.format(interval) + " seconds ..");*/
					
				} else {
					
					connectedCorrelationLength[i] = -1;
					numerator = 0;
					gSum = 0;
					
				}
				
				//calculate chi			
				dSum += numerator;
				gSumSum += gSum;			
				chi[i] = Math.sqrt(dSum / gSumSum);		
			}	
			
			IJ.freeMemory();IJ.freeMemory();
			
			//assign values to PoreClusterProperties
			mCP.connectedCorrelationLength = connectedCorrelationLength;
			mCP.chiOfPercolation = chi;
			
			return mCP;		
		}

	public double[] calculateTheBulkSoilVolume(String nowGaugePath, ImagePlus soilSurface, MenuWaiter.PoreSpaceAnalyzerOptions mPSA) {
		
		InputOutput jIO = new InputOutput();
		RoiHandler rH  = new RoiHandler();
		
		double[] bulkVolume = new double[4];
		
		//read InnerCircle file
		ObjectDetector jOD = new ObjectDetector();
		ObjectDetector.ColCoords3D jCO = jOD.new ColCoords3D();
		int versio = jIO.checkInnerCircleFileVersion(nowGaugePath);			
		if (versio == 0) jCO = jIO.readInnerCircleVer0(nowGaugePath);	
		else jCO = jIO.readInnerCircleVer1(nowGaugePath);
		
		PolygonRoi[] pRoi = rH.makeMeAPolygonRoiStack("inner", "tight", jCO, -mPSA.cutAwayFromWall);
		
		//assign soil top and bottom surfaces		
		ImageStack surStack = soilSurface.getStack();
		ImageProcessor topIP = surStack.getProcessor(1);		
		ImageProcessor botIP = surStack.getProcessor(2);
		
		//cut away outside of column of top ...
		topIP.setRoi(pRoi[0]);
		topIP.setColor(0);
		topIP.fillOutside(pRoi[0]);
		
		// ... and bottom
		botIP.setRoi(pRoi[0]);
		botIP.setColor(0);
		botIP.fillOutside(pRoi[pRoi.length - 1]);
				
		// count volume above the soil surface ...
		double volumeAboveColumn = 0;
		for (int x = 0 ; x < topIP.getWidth() ; x++) {
			for (int y = 0 ; y < topIP.getHeight() ; y++) {
				int surPix = topIP.getPixel(x, y);
				volumeAboveColumn += surPix;
			}
		}
		
		// .. and below..
		double volumeBelowColumn = 0;
		for (int x = 0 ; x < botIP.getWidth() ; x++) {
			for (int y = 0 ; y < botIP.getHeight() ; y++) {
				int surPix = botIP.getPixel(x, y);
				volumeBelowColumn += surPix;
			}
		}
		
		//calculate volume of cut-out column (in case the parts close to the walls have been cut away)
		double columnVolume = 0;
		for (int i = 0 ; i < jCO.heightOfColumn ; i++) {
			columnVolume += Math.round(Math.PI * (jCO.innerMajorRadius[i] - mPSA.cutAwayFromWall) * (jCO.innerMinorRadius[i]) - mPSA.cutAwayFromWall);						
		}
		
		//calculate volume of complete column
		double volume2Walls = 0;
		double originalVolume = 0;
		for (int i = 0 ; i < jCO.heightOfColumn ; i++) {
			originalVolume += Math.round(Math.PI * jCO.innerMajorRadius[i] * jCO.innerMinorRadius[i]);						
		}
		volume2Walls = originalVolume - columnVolume;
		
		bulkVolume[0] = columnVolume - volumeAboveColumn - volumeBelowColumn;
		bulkVolume[1] = volumeAboveColumn;
		bulkVolume[2] = volumeBelowColumn;
		bulkVolume[3] = volume2Walls;
		
		return bulkVolume;
		
	}

	public ArrayList<Integer> check4TouchingTheBottom(int imageNumber, InputOutput.MyFolderCollection mFC, ImagePlus nowTiff, boolean isRealColumn) {
	
		ArrayList<Integer> conBot = new ArrayList<Integer>();
	
		//get cluster image			
		int i,x,y;
		int iW = nowTiff.getWidth();
		int iH = nowTiff.getHeight();		
		int stackHeight = nowTiff.getNSlices();
		
		//define neighborhood variables
		int n;
		
		if (isRealColumn == true) {
			
			InputOutput jIO = new InputOutput(); 
			
			ImagePlus soilSurface = jIO.openTiff3D(mFC.mySurfaceFolder + "//" + mFC.mySurfaceFileNames[imageNumber]);
			soilSurface.setPosition(2);
			ImageProcessor botIP = soilSurface.getProcessor();
			int cc = 0;
			
			for (x = 0 ; x < iW ; x++) {
				for (y = 0 ; y < iH ; y++) {
					
					cc++;
					IJ.showStatus("Finding pore clusters connected to the bottom surface in kilo-pixel " + (cc / 1000) + "/" + (iW * iH / 1000));
					
					int myVox = botIP.getPixel(x, y);
					if (myVox > 0) {  //check if pixel is within the soil column outlines					
						for (int ix = -1 ; ix < 2 ; ix++) {  
							for (int iy = -1 ; iy < 2 ; iy++) {
								if (ix != 0 & iy != 0) {
									n = botIP.getPixel(x + ix, y + iy);  //get values of neighborhood 
									if (n > myVox) {
										for (i = 0 ; i < n - myVox ; i++) {
											nowTiff.setPosition(stackHeight - (myVox + i));
											ImageProcessor nowIP = nowTiff.getProcessor();
											int nPixel = (int)nowIP.getPixelValue(x + ix, y + iy);										
											if (nPixel > 0) if (!conBot.contains(nPixel)) conBot.add(nPixel);										
										}
									} else if (n > 0 & n < myVox) {
										for (i = 0 ; i < myVox - n ; i++) {
											nowTiff.setPosition(stackHeight - (myVox - i));
											ImageProcessor nowIP = nowTiff.getProcessor();
											int nPixel = (int)nowIP.getPixelValue(x + ix, y + iy);								
											if (nPixel > 0) if (!conBot.contains(nPixel)) conBot.add(nPixel);										
										}
									}
								}
							}
						}
					}
				}
			}	
		}
		
		if (isRealColumn == false) {
			
			//check bottom of sub-sample
			nowTiff.setPosition(stackHeight);
			ImageProcessor nowIP = nowTiff.getProcessor();
			int cc = 0;		
			for (x = 0 ; x < iW ; x++) {
				for (y = 0 ; y < iH ; y++) {				
					cc++;
					IJ.showStatus("Finding pore clusters connected to the bottom surface in kilo-pixel " + (cc / 1000) + "/" + (iW * iH / 1000));
					
					int nPixel = (int)nowIP.getPixelValue(x, y);										
					if (nPixel > 0) if (!conBot.contains(nPixel)) conBot.add(nPixel);	
				}
			}	
		}
		
		return conBot;
	
	}

	public ArrayList<Integer> check4TouchingTheTop(int imageNumber, InputOutput.MyFolderCollection mFC, ImagePlus nowTiff, boolean isRealColumn) {
	
		ArrayList<Integer> conTop = new ArrayList<Integer>();
		
		//get cluster image			
		int i,x,y;
		int iW = nowTiff.getWidth();
		int iH = nowTiff.getHeight();		
		
		//define neighborhood variables
		int n;
		
		if (isRealColumn == true) {
			
			InputOutput jIO = new InputOutput(); 
			
			ImagePlus soilSurface = jIO.openTiff3D(mFC.mySurfaceFolder + "//" + mFC.mySurfaceFileNames[imageNumber]);
			
			soilSurface.setPosition(1);
			ImageProcessor surIP = soilSurface.getProcessor();
			int cc = 0;
			
			for (x = 0 ; x < iW ; x++) {
				for (y = 0 ; y < iH ; y++) {
				
					cc++;
					IJ.showStatus("Finding pore clusters connected to the top surface in kilo-pixel " + (cc / 1000) + "/" + (iW * iH / 1000));
					
					int myVox = surIP.getPixel(x, y);
					if (myVox > 0) {  //check if pixel is within the soil column outlines					
						for (int ix = -1 ; ix < 2 ; ix++) {  
							for (int iy = -1 ; iy < 2 ; iy++) {
								if (ix != 0 & iy != 0) {
									n = surIP.getPixel(x + ix, y + iy);  //get values of neighborhood 
									if (n > myVox) {
										for (i = 0 ; i < n - myVox ; i++) {
											nowTiff.setPosition(myVox + i);
											ImageProcessor nowIP = nowTiff.getProcessor();
											int nPixel = (int)nowIP.getPixelValue(x + ix, y + iy);							
											if (nPixel > 0) if (!conTop.contains(nPixel)) conTop.add(nPixel);										
										}
									} else if (n > 0 & n < myVox) {
										for (i = 0 ; i < myVox - n ; i++) {
											nowTiff.setPosition(myVox - i);
											ImageProcessor nowIP = nowTiff.getProcessor();
											int nPixel = (int)nowIP.getPixelValue(x + ix, y + iy);								
											if (nPixel > 0) if (!conTop.contains(nPixel)) conTop.add(nPixel);										
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		if (isRealColumn == false) {
			
			int cc = 0;
			
			//check top of sub sample
			nowTiff.setPosition(1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			for (x = 0 ; x < iW ; x++) {
				for (y = 0 ; y < iH ; y++) {
					cc++;
					IJ.showStatus("Finding pore clusters connected to the top surface in kilo-pixel " + (cc / 1000) + "/" + (iW * iH / 1000));
				
					int nPixel = (int)nowIP.getPixelValue(x, y);							
					if (nPixel > 0) if (!conTop.contains(nPixel)) conTop.add(nPixel);
				}
			}
		}
		
		return conTop;
		
	}

	public boolean[] checkForPercolatingClusters(int nParticles, ArrayList<Integer> conTop, ArrayList<Integer> conBot) {
			
			int i;		
							
			//sort Lists
			Collections.sort(conTop);
			Collections.sort(conBot);
			
			//checking if a cluster is percolating..
			IJ.showStatus("Pinning down percolating clusters ...");
					
			//find percolating clusters
			ArrayList<Integer> conPerk= new ArrayList<Integer>();
	        for (int temp : conTop) {
	        	if (conBot.contains(temp)) conPerk.add(temp);
	        }
					
			//writing touches to cluster properties
			int numOfCluster = nParticles;
			boolean[] isPercolating = new boolean[numOfCluster];
			
			//init touching information
			for (i = 0 ; i < isPercolating.length ; i++) {
				if (conPerk.contains(i)) isPercolating[i] = true;
				else isPercolating[i] = false;
			}
		
			return isPercolating;	
		}

	public double[] calculateFractalDimension(ImagePlus nowTiff, InputOutput.MyFolderCollection mFC) {
				 
		JFractalCount myFBC = new JFractalCount();		
		double[] myFractalResults = myFBC.count(nowTiff);
		
		return myFractalResults;
		
	}
	
	public ProfileStatistics findProfileStatistics(ImagePlus nowTiff, PolygonRoi[] pRoi) {
		
		ProfileStatistics pS = new ProfileStatistics();		
		int[] numberOfNonZeros = new int[nowTiff.getNSlices()];
		double[] mean = new double[nowTiff.getNSlices()];	
		double[] geomean = new double[nowTiff.getNSlices()];	
		double[] logMean = new double[nowTiff.getNSlices()];
		double[] median = new double[nowTiff.getNSlices()];
		double[] mode = new double[nowTiff.getNSlices()];
		double[] std = new double[nowTiff.getNSlices()];
		double[] minVal = new double[nowTiff.getNSlices()];
		double[] maxVal = new double[nowTiff.getNSlices()];	
		
		for (int i = 1 ; i <= nowTiff.getNSlices() ; i++) {
			
			IJ.showStatus("Finding statistics of along vertical profile " + i + "/" + nowTiff.getNSlices() + "...");
		
			//move to next slice
			nowTiff.setPosition(i);
			ImageProcessor nowIP = nowTiff.getProcessor();
		
			//sample all grey values in a float array list
			ArrayList<Float> myValues = new ArrayList<Float>();			
			for (int x = 0 ; x < nowIP.getWidth() ; x++) {
				for (int y = 0 ; y < nowIP.getHeight() ; y++) {
					if (nowIP.getPixelValue(x, y) > 0) myValues.add(nowIP.getPixelValue(x, y));										
				}
			}
		
			//convert the float array list to an array
			double[] greyValues = new double[myValues.size()];
			for (int j = 0 ; j < myValues.size(); j++) greyValues[j] = myValues.get(j); 
			
			//calculate statistical values from the float array
			numberOfNonZeros[i-1] = myValues.size();
			mean[i-1] = StatUtils.mean(greyValues);
			logMean[i-1] = StatUtils.sumLog(greyValues)/numberOfNonZeros[i-1];
			median[i-1] = StatUtils.percentile(greyValues,50);			
			geomean[i-1] = StatUtils.geometricMean(greyValues);			
			std[i-1] = Math.sqrt(StatUtils.variance(greyValues));
			maxVal[i-1] = StatUtils.max(greyValues);	
			minVal[i-1] = StatUtils.min(greyValues);	
		}	
		
		pS.numberOfNonZeros = numberOfNonZeros;
		pS.mean = mean;
		pS.logMean = logMean;
		pS.median = median;
		pS.geomean = geomean;
		pS.mode = mode;
		pS.std = std;
		pS.maxi = maxVal;
		pS.mini = minVal;
		
		return pS;
	}

	public void tailoredPoreSpaceAnalyses(int imageNumber, InputOutput.MyFolderCollection mFC, ImagePlus[] nowTiffs, MenuWaiter.PoreSpaceAnalyzerOptions mPSA) {
				
		InputOutput jIO = new InputOutput();
		JParticleCounter jJPA = new JParticleCounter();
		ImageManipulator jIM = new ImageManipulator();
		PoreClusterProperties myP = new PoreClusterProperties();		
		AndAllTheRest aa = new AndAllTheRest();
		
		//init some very basic variables..		
		String nowImageName = mFC.colName;
		
		//create output folders and save statistics		
		IJ.showStatus("Creating directories for the pore-cluster properties ...");
		String myOutFolder = "Stats";
		String myOutPath = mFC.myPreOutFolder + "\\" + myOutFolder;
		new File(myOutPath).mkdir();		
		String outPath = mFC.myPreOutFolder + "\\Stats\\" + nowImageName + "_BJResults.txt";
				
		//input parameters for particle analyzer
		int FORE = -1;
		int slicesPerChunk = 2; //input parameter 4 particle analyzer
		double minVol = 0; //input parameter 4 particle analyzer  ..364 corresponds to 100,000,000 cubic micrometers or 0.1 cubic millimeter
		double maxVol = Double.POSITIVE_INFINITY; //input parameter 4 particle analyzer
		Boolean doExclude = false;		
		
		//check is it is a real column with non-plane top and bottom surface
		boolean isRealColumn = false;
		if (nowTiffs[1]!=null) isRealColumn = true;
		
		//check if the clusters need to be identified
		myP.containsAParticleAnalysis = false;
		if (mPSA.performParticleAnalyses == true) {
		
			myP.containsAParticleAnalysis = true;
			
			// get the particles	
			IJ.showStatus("Identifying connected pore-clusters ...");
			Object[] result = jJPA.getParticles(nowTiffs[0], slicesPerChunk, minVol, maxVol, FORE, doExclude);
			IJ.freeMemory();IJ.freeMemory();
			
			//extract the basic information
			int[][] particleLabels = (int[][]) result[1];
			long[] particleSizes = jJPA.getParticleSizes(particleLabels);
			final int nParticles = particleSizes.length;
			double[][] centroids = jJPA.getCentroids(nowTiffs[0], particleLabels, particleSizes);
			int[][] limits = jJPA.getParticleLimits(nowTiffs[0], particleLabels, nParticles);
			IJ.freeMemory();IJ.freeMemory();
			
			//volumes		 
			double[] volumes = jJPA.getVolumes(nowTiffs[0], particleSizes);
			IJ.freeMemory();IJ.freeMemory();
			
			/*// surfaces
			ArrayList<List<Point3f>> surfacePoints = new ArrayList<List<Point3f>>();
			double[] surfaceAreas = new double[nParticles];
			int resampling = 6;
			if (mPSA.calcSurface == true) {
				//surfacePoints = jJPA.getSurfacePoints(nowTiffs[0], particleLabels, limits, resampling, nParticles);
				//surfaceAreas = jJPA.getSurfaceArea(surfacePoints);				
			}		
			IJ.freeMemory();IJ.freeMemory();*/
			
						
			//Euler
			double[][] eulerCharacters = new double[nParticles][3];		
			if (mPSA.calcEuler == true) eulerCharacters = jJPA.getEulerCharacter(nowTiffs[0], particleLabels, limits, nParticles);
			IJ.freeMemory();IJ.freeMemory();
			
			//calculate thicknesses
			double[][] thick = new double[nParticles][2];		
			Thickness th = new Thickness();
			ImagePlus thickImp = null;
			if (mPSA.plotThickness == true | mPSA.plotPercolation == true | mPSA.calcThickness == true | mPSA.calcCriticalPoreDiameter == true) {			
				
				double max = 0;
				
				//create thickness image
				ImagePlus rawThickImp = new ImagePlus();
				if (nowTiffs[1] != null) rawThickImp = th.getLocalThickness(nowTiffs[1], false);
				else rawThickImp = th.getLocalThickness(nowTiffs[0], false);
				
				//cut out pore space again.. to circumvent thickness bug..
			    thickImp = jIM.cutOutRealThickness(nowTiffs[0], rawThickImp);
				thickImp.getProcessor().setMinAndMax(0, max);
				thickImp.setTitle(nowTiffs[0].getShortTitle() + "_thickness");
				
				//calculate average std and max thicknesses..
				thick = jJPA.getMeanStdDev(rawThickImp, particleLabels, particleSizes, 0);			
				for (int i = 1; i < nParticles; i++) max = Math.max(max, thick[i][2]);
				
				rawThickImp.killStack();
				IJ.freeMemory();IJ.freeMemory();
			}
			
			//calculate percolating clusters ... needs the other parameters.. and is needed for following parameters
			ImagePlus particleLabelTiff = jJPA.displayParticleLabels(particleLabels, nowTiffs[0]);
			ArrayList<Integer> conTop = check4TouchingTheTop(imageNumber, mFC, particleLabelTiff, isRealColumn);	
			ArrayList<Integer> conBot = check4TouchingTheBottom(imageNumber, mFC, particleLabelTiff, isRealColumn);
			boolean[] cTop = new boolean[nParticles];
			for (int i = 0 ; i < nParticles ; i++) if (conTop.contains(i)) cTop[i] = true;
			boolean[] cBot = new boolean[nParticles];
			for (int i = 0 ; i < nParticles ; i++) if (conBot.contains(i)) cBot[i] = true;
			boolean[] isPercolating = checkForPercolatingClusters(nParticles, conTop, conBot);
						
			//calc critical pore diameter
			IJ.freeMemory();IJ.freeMemory();			
			
			if (mPSA.calcCriticalPoreDiameter == true) {
								
				IJ.showStatus("Calculating the critical pore diameter ...");
				
				//check which clusters are percolating and write them into a list..
				ArrayList<Integer> percolatingClusters = new ArrayList<Integer>();
				for (int i = 0 ; i < isPercolating.length ; i++) if (isPercolating[i] == true) percolatingClusters.add(i);
				
				//if there is not even one cluster, the critical pore diameters is <= image resolution
				myP.criticalPoreDiameter = -1;
				if (!percolatingClusters.isEmpty()) {  //calculate critical pore diameter										
					
					//creating a list of all pore diameters within a percolating cluster
					ImagePlus thicknessOfPercolatingClusters = new ImagePlus();
					ImageStack thickPercClust = new ImageStack(thickImp.getWidth(), thickImp.getHeight());
					
					ArrayList<Double> poreDiameters = new ArrayList<Double>();		//unique pore diameters
					
					//ignore topmost cross-section if the thickness image was derived from a cuboid or cylindrical ROI
					int startZ = 1;
					if (mPSA.choiceOfRoi.equalsIgnoreCase("Cuboid") | mPSA.choiceOfRoi.equalsIgnoreCase("Cylinder")) startZ = 2;
					
					//loop through horizontal cross-sections
					for (int z = startZ ; z <= thickImp.getNSlices() ; z++) {				
						
						thickImp.setPosition(z);
						ImageProcessor nowIP = thickImp.getProcessor();
						
						//thickImp.updateAndDraw();
						//thickImp.show();
						
						particleLabelTiff.setPosition(z);
						ImageProcessor clustIP = particleLabelTiff.getProcessor();
						
						//particleLabelTiff.updateAndDraw();
						//particleLabelTiff.show();
						
						ImageProcessor outIP = nowIP.duplicate();
						
						//loop through all pixel within cross-section
						for (int x = 0 ; x < thickImp.getWidth() ; x++) {
							for (int y = 0 ; y < thickImp.getHeight() ; y++) {
								
								//check if pixel is contained in a percolating cluster
								int nowCluster = (int)clustIP.getPixelValue(x, y);
								if (percolatingClusters.contains(nowCluster)) {																
									
									double nowPixel = nowIP.getPixelValue(x, y);
									outIP.putPixelValue(x, y, nowPixel); 
																		
									//make an inventory of existing pore diameters..
									boolean isContained = false;
									double myRoundPix = Math.round(nowPixel);
									if (!poreDiameters.isEmpty()) {
										for (int j = 0 ; j < poreDiameters.size() ; j++) if (poreDiameters.get(j) == myRoundPix) isContained = true;										
									}
									if (!isContained) poreDiameters.add(myRoundPix);									
								}
								else outIP.putPixelValue(x, y, 0);   //cut away if not within percolating cluster							
							}							
						}
						
						thickPercClust.addSlice(outIP);
					}
					
					//sort the list and create a new image of just the thicknesses of the percolating pore clusters.. 
					Collections.sort(poreDiameters);
					thicknessOfPercolatingClusters.setStack(thickPercClust);		
					
					//thicknessOfPercolatingClusters.updateAndDraw();
					//thicknessOfPercolatingClusters.show();
					
					IJ.freeMemory();IJ.freeMemory();
								
					//find the smallest of the maximal pore diameters of all cross-sections and use it as starting point
					ProfileStatistics pS = findProfileStatistics(thicknessOfPercolatingClusters, null);				
					double minVertThick = StatUtils.min(pS.maxi);
					
					//set threshold to pore diameters == minVertThick 
					double threshold = minVertThick;
					int thresholdIndex = 0;
					boolean laPerco = false;
					for (int i = 0 ; i < poreDiameters.size() ; i++) if (poreDiameters.get(i) == Math.round(minVertThick)) thresholdIndex = i;
					
					//loop while precolation is false
					while (!laPerco & thresholdIndex >= 0) {
					
						//set all values smaller than the threshold to 0
						ImageStack largerPores = new ImageStack(thickImp.getWidth(), thickImp.getHeight());
						ImagePlus laPores = new ImagePlus();
						for (int z = 1 ; z <= thicknessOfPercolatingClusters.getNSlices() ; z++) {	
							
							thicknessOfPercolatingClusters.setPosition(z);
							ImageProcessor nowIP = thicknessOfPercolatingClusters.getProcessor();
							ImageProcessor laIP = new ByteProcessor(thickImp.getWidth(), thickImp.getHeight());
							
							//loop through all pixel within cross-section
							for (int x = 0 ; x < thickImp.getWidth() ; x++) {
								for (int y = 0 ; y < thickImp.getHeight() ; y++) {							
									double nowPixel = nowIP.getPixelValue(x, y);
									if (nowPixel >= threshold) laIP.putPixel(x, y, 255);
									else laIP.putPixel(x, y, 0);							
								}
							}					
							largerPores.addSlice(laIP);
						}		
						laPores.setStack(largerPores);
						
						//laPores.updateAndDraw();
						//laPores.show();
										
						IJ.freeMemory();IJ.freeMemory();
						
						//find pore clusters in laPores and check if they form a percolating path
						Object[] largestClusters = jJPA.getParticles(laPores, slicesPerChunk, minVol, maxVol, FORE, doExclude);
						int[][] laParticlesLabels = (int[][]) largestClusters[1];
						ImagePlus laParLaTiff = jJPA.displayParticleLabels(laParticlesLabels, laPores);					
						long[] laParticleSizes = jJPA.getParticleSizes(particleLabels);
						final int laNParticles = laParticleSizes.length;
						ArrayList<Integer> laConTop = check4TouchingTheTop(imageNumber, mFC, laParLaTiff, isRealColumn);	
						ArrayList<Integer> laConBot = check4TouchingTheBottom(imageNumber, mFC, laParLaTiff, isRealColumn);
						boolean[] laP = checkForPercolatingClusters(laNParticles, laConTop, laConBot);
						percolatingClusters.clear();
						for (int i = 0 ; i < laP.length ; i++) if (laP[i] == true) percolatingClusters.add(i);
						
						//check if there is a percolating cluster and try again with a smaller threshold in case of not
						if (percolatingClusters.isEmpty()) {
							thresholdIndex--;
							if (thresholdIndex < 0) {
								myP.criticalPoreDiameter = 0;
								break;
							}
							threshold = poreDiameters.get(thresholdIndex);
						}
						else {
							laPerco = true;  // this will break the loop..
							myP.criticalPoreDiameter = threshold;
						}
					}					
				}
			}
			
			//free memory	
			IJ.freeMemory();IJ.freeMemory();
						
			//moments, unit vectors and orientation 
			EigenvalueDecomposition[] eigens = new EigenvalueDecomposition[nParticles];		
			if (mPSA.calcMoments == true | mPSA.calcUnitVectors == true)  eigens = jJPA.getEigens(nowTiffs[0], particleLabels, centroids);	
			
			// init variables for export of results
			int[] id = new int[nParticles - 1];									//ID: unique particle identifier; this number is the label used for the particle in all calculations and output
			double[] volume = new double[nParticles - 1];						//Vol: particle volume
			double[] xCenter = new double[nParticles - 1];						//x Cent: x-coordinate of particle centroid
			double[] yCenter = new double[nParticles - 1];						//y Cent: y-coordinate of particle centroid
			double[] zCenter = new double[nParticles - 1];						//z Cent: z-coordinate of particle centroid
			int[] xMin = new int[nParticles - 1];								//x Min: min x-coordinate of particle 
			int[] yMin = new int[nParticles - 1];								//y Min: min y-coordinate of particle 
			int[] zMin = new int[nParticles - 1];								//z Min: min z-coordinate of particle 
			int[] xMax = new int[nParticles - 1];								//x Max: max x-coordinate of particle 
			int[] yMax = new int[nParticles - 1];								//y Max: max y-coordinate of particle 
			int[] zMax = new int[nParticles - 1];								//z Max: max z-coordinate of particle 
			double[] surfaceArea = new double[nParticles - 1];					//SA: surface area (0 if too small for mesh to be produced; see warning log)		
			boolean[] tTop = new boolean[nParticles - 1];						//cluster is touching top..		
			boolean[] tBot = new boolean[nParticles - 1];						//cluster is touching bot..		
			boolean[] percolating = new boolean[nParticles - 1];				//cluster is percolating or not..			
			double[] momentOfInertiaShortestAxis = new double[nParticles - 1];	//I1: moment of inertia around shortest principal axis
			double[] momentOfInertiamiddleAxis = new double[nParticles - 1];	//I2: moment of inertia around middle principal axis
			double[] momentOfInertiaLongestAxis = new double[nParticles - 1];	//I3: moment of inertia around longest principal axis
			double[] unitVectorInXDirection = new double[nParticles - 1];		//vX: x component of unit vector of longest principal axis (a measure of orientation)
			double[] unitVectorInYDirection = new double[nParticles - 1];		//vY: y component of unit vector of longest principal axis
			double[] unitVectorInZDirection = new double[nParticles - 1];		//vZ: z component of unit vector of longest principal axis
			double[] euler = new double[nParticles - 1];						//Euler (xi): Euler characteristic of the particle
			double[] holes = new double[nParticles - 1];						//Holes (beta1): number of topological holes (handles) in the particle
			double[] cavities = new double[nParticles - 1];						//Cavities (beta2): number of enclosed cavities in the particle
			double[] thickness = new double[nParticles - 1];					//Thickness: mean local thickness of particle
			double[] sdThickness = new double[nParticles - 1];					//SD Thickness: standard deviation of the mean local thickness of particle
			double[] maxThickness = new double[nParticles - 1];					//Max Thickness: maximum local thickness of particle		
				
			//sort results 
			int[] sortedIndices = aa.getIndicesInOrder(volumes);
			
			//write results into a biiiiig table		
			for (int i = 1; i < volumes.length; i++) {
				if (volumes[i] > 0) {			
					
					IJ.showStatus("Sorting out properties of the individual pore-clusters ...");
					
					id[i-1] = sortedIndices[i];				
					volume[i-1] = volumes[sortedIndices[i]];
					xCenter[i-1] = centroids[sortedIndices[i]][0];				
					yCenter[i-1] = centroids[sortedIndices[i]][1];
					zCenter[i-1] = centroids[sortedIndices[i]][2];
				
					xMin[i-1] = (int)Math.floor(limits[sortedIndices[i]][0]);				
					yMin[i-1] = (int)Math.floor(limits[sortedIndices[i]][2]);
					zMin[i-1] = (int)Math.floor(limits[sortedIndices[i]][4]);
					xMax[i-1] = (int)Math.floor(limits[sortedIndices[i]][1]);				
					yMax[i-1] = (int)Math.floor(limits[sortedIndices[i]][3]);
					zMax[i-1] = (int)Math.floor(limits[sortedIndices[i]][5]);
					
					//if (mPSA.calcSurface == true) surfaceArea[i-1] = surfaceAreas[sortedIndecies[i]];		
					
					tTop[i-1] = cTop[sortedIndices[i]];
					tBot[i-1] = cBot[sortedIndices[i]];
					percolating[i-1] = isPercolating[sortedIndices[i]];
									
					EigenvalueDecomposition E = eigens[sortedIndices[i]];
										
					if (mPSA.calcUnitVectors == true) {
						unitVectorInXDirection[i-1] = E.getV().get(0, 0);				
						unitVectorInYDirection[i-1] = E.getV().get(1, 0);
						unitVectorInZDirection[i-1] = E.getV().get(2, 0);
					}
					
					if (mPSA.calcEuler == true) {
						euler[i-1] = eulerCharacters[sortedIndices[i]][0];				
						holes[i-1] = eulerCharacters[sortedIndices[i]][1];
						cavities[i-1] = eulerCharacters[sortedIndices[i]][2];
					}
					
					if (mPSA.calcThickness == true) {
						thickness[i-1] = thick[sortedIndices[i]][0];				
						sdThickness[i-1] = thick[sortedIndices[i]][1];
						maxThickness[i-1] = thick[sortedIndices[i]][2];
					}				
				}
			}
	
			//assign results to output structure		
			myP.id = id;
			myP.volume = volume;
			myP.xCenter = xCenter;
			myP.yCenter = yCenter;		
			myP.zCenter = zCenter;
			myP.minX = xMin;
			myP.minY = yMin;
			myP.minZ = zMin;
			myP.maxX = xMax;
			myP.maxY = yMax;
			myP.maxZ = zMax;
			myP.surfaceArea = surfaceArea;	
			myP.touchesTop = tTop;
			myP.touchesBot = tBot;
			myP.isPercolating = percolating;
			myP.momentOfInertiaLongestAxis = momentOfInertiaLongestAxis;
			myP.momentOfInertiamiddleAxis = momentOfInertiamiddleAxis;
			myP.momentOfInertiaShortestAxis = momentOfInertiaShortestAxis;
			myP.unitVectorInXDirection = unitVectorInXDirection;
			myP.unitVectorInYDirection = unitVectorInYDirection;
			myP.unitVectorInZDirection = unitVectorInZDirection;
			myP.euler = euler;
			myP.holes = holes;		
			myP.cavities = cavities;
			myP.thickness = thickness;
			myP.sdThickness = sdThickness;
			myP.maxThickness = maxThickness;
			
			//calculate connected correlation length and chi
			if (mPSA.calcChi == true) {						
				myP = calculateChiOfPercolation(particleLabelTiff, mFC, myP);
			}
			else {
				double[] dummy = new double[myP.id.length];
				for (int i = 0 ; i < myP.id.length ; i++) dummy[i] = -1;				
				myP.connectedCorrelationLength = dummy;
				myP.chiOfPercolation = dummy;				
			}
			
			IJ.freeMemory();IJ.freeMemory();
			
			//calculate global connection probability
			double sumOfSquaredClusterSizes = 0;
			double rF = 1000000;
			for (int i = 0 ; i < myP.id.length ; i++) {
				double squaredClusterSize = myP.volume[i] / rF * myP.volume[i] / rF;
				sumOfSquaredClusterSizes += squaredClusterSize;
			}
			myP.globVolume = macroPoreVolume(nowTiffs[0]);
			myP.globalConnectionProbability = sumOfSquaredClusterSizes / (myP.globVolume / rF * myP.globVolume / rF);			
			IJ.freeMemory();IJ.freeMemory();
			
			//save desired images
			if (mPSA.plotLabels == true) {
				myOutFolder = "ClusterLabels";
				myOutPath = mFC.myPreOutFolder + "\\" + myOutFolder;
				new File(myOutPath).mkdir();				
				String clusterLabelImageName = nowImageName + "_ClusterLables.tif";
				jIO.tiffSaver(myOutPath, clusterLabelImageName, particleLabelTiff);
				if (!mPSA.plotPercolation) particleLabelTiff.killStack();				
				IJ.freeMemory();IJ.freeMemory();				
			}
			
			if (mPSA.plotThickness == true) {			
				myOutFolder = "PoreThick";
				myOutPath = mFC.myPreOutFolder + "\\" + myOutFolder;
				new File(myOutPath).mkdir();
				String nowDir = mFC.myPreOutFolder + "\\PoreThick";		
				String thicknessImageName = nowImageName + "_Thickness.tif";
				jIO.tiffSaver(nowDir, thicknessImageName, thickImp);	
				thickImp.killStack();
				IJ.freeMemory();IJ.freeMemory();
			}
			
			if (mPSA.plotVolume == true) {
				myOutFolder = "Volume";
				myOutPath = mFC.myPreOutFolder + "\\" + myOutFolder;
				new File(myOutPath).mkdir();
				String nowDir = mFC.myPreOutFolder + "\\Volume";		
				String volumeImageName = nowImageName + "_Volume.tif";
				ImagePlus volImp = jJPA.displayParticleValues(particleLabelTiff, particleLabels, volumes, "Pore-Cluster Volumes");				
				jIO.tiffSaver(nowDir, volumeImageName, volImp);
				volImp.killStack();
				IJ.freeMemory();IJ.freeMemory();
			}
			
			if (mPSA.plotPercolation == true) {
				myOutFolder = "PercolatingVolume";
				myOutPath = mFC.myPreOutFolder + "\\" + myOutFolder;
				new File(myOutPath).mkdir();
				String nowDir = mFC.myPreOutFolder + "\\" + myOutFolder;		
				String volumeImageName = nowImageName + "_PercVol.tif";				
				ArrayList<Integer> percolatingClusters = new ArrayList<Integer>();
				for (int i = 0 ; i < percolating.length ; i++) if (percolating[i]) percolatingClusters.add(id[i]);	
				ImagePlus percPorosityImp = extractPercolatingPorosity(particleLabelTiff, percolatingClusters);
				jIO.tiffSaver(nowDir, volumeImageName, percPorosityImp);
				percPorosityImp.killStack();
				particleLabelTiff.killStack();
				IJ.freeMemory();IJ.freeMemory();
			}
			
		}		//closes if-tag checking whether the 
		
		IJ.freeMemory();IJ.freeMemory();
		
		//calculate fractal dimension
		if (mPSA.calcFractal == true) {
			double[] myFractalResults = calculateFractalDimension(nowTiffs[0], mFC);		
			myP.fractalDimension = myFractalResults[0];	
			myP.fractalDimensionFitR2 = myFractalResults[1];	
		}
		
		IJ.freeMemory();IJ.freeMemory();
		
		//calculate anisotropy
		if (mPSA.globAnisotropy == true) {
			JAnisotropy myA = new JAnisotropy();
			myP.globAnisotropy = myA.autoCalculate(nowTiffs[0]);
		}
		
		IJ.freeMemory();IJ.freeMemory();
		
		//calculate macroporosity
		if (myP.globVolume == 0) myP.globVolume = macroPoreVolume(nowTiffs[0]);		
		
		double[] bulkSoilVolume = {0, 0, 0, 0};
		if (mPSA.choiceOfRoi.equalsIgnoreCase("RealSample") & mPSA.includeSurfaceTopography) {
			int[] myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(mFC.myTiffs[imageNumber], jIO.listInnerCircleFiles(mFC.myInnerCircleFolder, ""), mFC.mySurfaceFileNames);
			ImagePlus mySurfaceFile = jIO.openTiff3D(mFC.mySurfaceFolder + "\\" + mFC.mySurfaceFileNames[myGandS[1]]);
			bulkSoilVolume = calculateTheBulkSoilVolume(mFC.myInnerCircleFiles[myGandS[0]], mySurfaceFile, mPSA);
		}
		else {
			bulkSoilVolume[0] = (int)Math.round(mPSA.areaOfInterest) * nowTiffs[0].getNSlices();
		}
		
		myP.soilBulkVolume = bulkSoilVolume[0];
		myP.volumeRemovedAboveTopSurface = bulkSoilVolume[1];
		myP.volumeRemovedBelowBottomSurface = bulkSoilVolume[2];		
		myP.volumeRemovedFromWalls = bulkSoilVolume[3];	
		
		myP.macroPorosity = myP.globVolume / myP.soilBulkVolume;	
				
		IJ.freeMemory();IJ.freeMemory();
		
		//write Results
		jIO.writePoreClusterProperties2File(nowImageName, outPath, myP);		
		IJ.freeMemory();IJ.freeMemory();
		
	}

	public double tMaxLaplace(ImageStack nowStack, ImageStack lapStack, int myMode) {
		
		HistogramStuff hist = new HistogramStuff();
		ImageManipulator jIM = new ImageManipulator();
		
		ImagePlus lapTiff = new ImagePlus();
		lapTiff.setStack(lapStack);
		
		ImagePlus origTiff = new ImagePlus();
		origTiff.setStack(nowStack);
	
		//get histogram of the gradient image
		int[] myHistGradient = hist.sampleHistogram(lapTiff);				
		myHistGradient[0] = 0;
		int myCutoffGradient = hist.findTheKnee(myHistGradient);
		
		//segment gradient image..
		ImagePlus myGradientMask = lapTiff.duplicate();
		myGradientMask = jIM.binarizeGradientMask(myGradientMask, myCutoffGradient);
		//myGradientMask.show();
		
		//sample upper threshold for this segment
		double weightSum = 0;
		double weightAndSampleSum = 0;
		for (int i = 0 ; i < lapStack.getSize() ; i++) {
			
			IJ.showStatus("Calculating threshold from Laplace mask slice " + (i+1) + "/" + lapTiff.getNSlices());
			
			//set all image to the correct position
			origTiff.setPosition(i+1); //remember that this image has one slice more on top and bottom, respectively..
			lapTiff.setPosition(i+1);
			myGradientMask.setPosition(i+1);
			
			//get each images IPs
			ImageProcessor origIP = origTiff.getProcessor();
			ImageProcessor gradIP = lapTiff.getProcessor();
			ImageProcessor maskIP = myGradientMask.getProcessor();
	
			//because imageCalculator does not work... do it by hand.. again..
			for (int x = 0 ; x < origIP.getWidth(); x++) {
				for (int y = 0 ; y < origIP.getHeight(); y++) {
					if (maskIP.get(x, y) > 0 && origIP.get(x, y) > 0 && origIP.get(x, y) < myMode) {
						double nowWeight = gradIP.getPixelValue(x, y);
						double nowSample = origIP.getPixelValue(x, y);
						weightSum = weightSum + nowWeight;
						weightAndSampleSum = weightAndSampleSum + nowWeight * nowSample;
					}							
				}
			}
		}		
		
		double tMax = weightAndSampleSum / weightSum;
		
		return tMax;
	}

	public double tMaxSobel(ImageStack nowStack, ImageStack sobStack, int myMode) {
	
		//init units
		HistogramStuff hist = new HistogramStuff();
		ImageManipulator jIM = new ImageManipulator();
		
		//init variables
		ImagePlus gradTiff = new ImagePlus();
		gradTiff.setStack(sobStack);
		
		ImagePlus torigTiff = new ImagePlus();
		torigTiff.setStack(nowStack);
		
		//get histogram of the gradient image
		int[] myHistGradient = hist.sampleHistogram(gradTiff);				
		myHistGradient[0] = 0;
		int myCutoffGradient = hist.findTheKnee(myHistGradient);
		
		//segment gradient image..
		ImagePlus myGradientMask = gradTiff.duplicate();
		myGradientMask = jIM.binarizeGradientMask(myGradientMask, myCutoffGradient);
		//myGradientMask.show();
		
		//sample upper threshold for this segment
		double weightSum = 0;
		double weightAndSampleSum = 0;
		for (int i = 0 ; i < sobStack.getSize() ; i++) {
			
			IJ.showStatus("Calculating threshold from Sobel mask slice " + (i+1) + "/" + gradTiff.getNSlices());
			
			//set all image to the correct position
			torigTiff.setPosition(i+1); //remember that this image has one slice more on top and bottom, respectively..
			gradTiff.setPosition(i+1);
			myGradientMask.setPosition(i+1);
			
			//get each images IPs
			ImageProcessor origIP = torigTiff.getProcessor();
			ImageProcessor gradIP = gradTiff.getProcessor();
			ImageProcessor maskIP = myGradientMask.getProcessor();		
	
			//because imageCalculator does not work... do it by hand.. again..
			for (int x = 0 ; x < origIP.getWidth(); x++) {
				for (int y = 0 ; y < origIP.getHeight(); y++) {
					if (maskIP.get(x, y) > 0 && origIP.get(x, y) > 0 && origIP.get(x, y) < myMode) {
						double nowWeight = gradIP.getPixelValue(x, y);
						double nowSample = origIP.getPixelValue(x, y);
						weightSum = weightSum + nowWeight;
						weightAndSampleSum = weightAndSampleSum + nowWeight * nowSample;
					}							
				}
			}
		}		
		
		double tMax = weightAndSampleSum / weightSum;
		
		return tMax;
	}

	public double macroPoreVolume(ImagePlus nowTiff) {
		
		double macroPoreVolume = 0;
		
		for (int i = 0 ; i < nowTiff.getNSlices() ; i++) {
			
			nowTiff.setPosition(i+1);
			ImageProcessor myIP = nowTiff.getProcessor();
			int[] nowHist = myIP.getHistogram();
			
			macroPoreVolume += nowHist[255];
						
		}
		
		return macroPoreVolume;
	}

	public ImagePlus findTheLargestCluster(ImagePlus nowTiff, int indexOfLargestCluster) {
		
		int i;
		ImageStack lStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
		ImagePlus lTiff = new ImagePlus();
		
		//nowTiff.draw();
		//nowTiff.show();
		
		for (i = 0 ; i < nowTiff.getNSlices() ; i++) {
			
			nowTiff.setPosition(i + 1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			IJ.showStatus("Isolating largest cluster in slice #" + (i + 1) + "/" + (nowTiff.getNSlices()+1));				
			
			ImageProcessor modIP = new ByteProcessor(nowTiff.getWidth(), nowTiff.getHeight());
			
			for (int x = 0 ; x < nowIP.getWidth() ; x++) {
				for (int y = 0 ; y < nowIP.getHeight() ; y++) {
					int nowPix = Math.round(nowIP.getPixelValue(x, y));					
					if (nowPix == indexOfLargestCluster) modIP.putPixelValue(x, y, 255);
					else modIP.putPixelValue(x, y, 0);					
				}			
			}								
			lStack.addSlice(modIP);
		}
		lTiff.setStack(lStack);
		
		return lTiff;
	}
	
	public class PoreClusterProperties {
	
		//integral properties
		public double soilBulkVolume;
		public double volumeRemovedAboveTopSurface;
		public double volumeRemovedBelowBottomSurface;
		public double volumeRemovedFromWalls;		
		public double globVolume;
		public double macroPorosity;
		public double globSurface;
		public double globThickness;
		public double criticalPoreDiameter;
		public double chi;		
		public double fractalDimension;					//fractal dimension.. (box counting!)
		public double fractalDimensionFitR2;
		public double globAnisotropy;
		public double globalConnectionProbability; 
					
		//cluster-wise properties
		public boolean containsAParticleAnalysis;
		public int[] id;								//ID: unique particle identifier; this number is the label used for the particle in all calculations and output
		public double[] volume;						//Vol: particle volume
		public double[] xCenter;						//x Cent: x-coordinate of particle centroid
		public double[] yCenter;						//y Cent: y-coordinate of particle centroid
		public double[] zCenter;						//z Cent: z-coordinate of particle centroid
		public double[] surfaceArea;					//SA: surface area (0 if too small for mesh to be produced; see warning log)
		public double[] enclosedVolume;				//Encl. Vol: Volume enclosed by surface mesh (0 if too small for mesh to be produced; see warning log)
		public double[] momentOfInertiaShortestAxis;	//I1: moment of inertia around shortest principal axis
		public double[] momentOfInertiamiddleAxis;	//I2: moment of inertia around middle principal axis
		public double[] momentOfInertiaLongestAxis;	//I3: moment of inertia around longest principal axis
		public double[] unitVectorInXDirection;		//vX: x component of unit vector of longest principal axis (a measure of orientation)
		public double[] unitVectorInYDirection;		//vY: y component of unit vector of longest principal axis
		public double[] unitVectorInZDirection;		//vZ: z component of unit vector of longest principal axis
		public double[] euler;						//Euler (xi): Euler characteristic of the particle
		public double[] holes;						//Holes (beta1): number of topological holes (handles) in the particle
		public double[] cavities;					//Cavities (beta2): number of enclosed cavities in the particle
		public double[] thickness;					//Thickness: mean local thickness of particle
		public double[] sdThickness;					//SD Thickness: standard deviation of the mean local thickness of particle
		public double[] maxThickness;				//Max Thickness: maximum local thickness of particle		
		public double[] majorRadius;					//Major radius: length of best-fit ellipsoid's long radius
		public double[] intRadius;					//Int. radius: length of best-fit ellipsoid's intermediate radius
		public double[] minorRadius;					//Minor radius: length of best-fit ellipsoid's short radius		
					
		public double[] chiOfPercolation;
		public double[] connectedCorrelationLength;
		
		public boolean[] isPercolating;
		public boolean[] touchesTop;
		public boolean[] touchesBot;
		
		public int[] minZ;
		public int[] maxZ;
		public int[] minX; 
		public int[] maxX;
		public int[] minY;
		public int[] maxY;
	}
}

