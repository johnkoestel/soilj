package soilj.tools;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.Vector;

import org.scijava.vecmath.Point3f;

import sc.fiji.localThickness.*;
import soilj.copiedTools.JAnisotropy;
import soilj.copiedTools.JFractalCount;
import soilj.copiedTools.JParticleCounter;
import soilj.tools.InputOutput;
import soilj.tools.ObjectDetector;
import Jama.EigenvalueDecomposition.*;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import marchingcubes.MCTriangulator;
import Utilities.Counter3D;
import Utilities.Object3D;

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

	public int[] calculateTheBulkSoilVolume(String myPreOutPath, String nowGaugePath, ImagePlus soilSurface) {
		
		InputOutput jIO = new InputOutput();
		RoiHandler rH  = new RoiHandler();
		
		int[] bulkVolume = new int[3];
		
		//create output folders		
		String myOutFolder = "Stats";
		String myOutPath = myPreOutPath + "\\" + myOutFolder;
		new File(myOutPath).mkdir();
		
		//read gauge file
		ObjectDetector.ColumnCoordinates jCO = jIO.readGaugeFile(nowGaugePath);	
		PolygonRoi[] pRoi = rH.makeMeAPolygonRoiStack("inner", "tight", jCO, 0);
		
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
		int volumeAboveColumn = 0;
		for (int x = 0 ; x < topIP.getWidth() ; x++) {
			for (int y = 0 ; y < topIP.getHeight() ; y++) {
				int surPix = topIP.getPixel(x, y);
				volumeAboveColumn += surPix;
			}
		}
		
		// .. and below..
		int volumeBelowColumn = 0;
		for (int x = 0 ; x < botIP.getWidth() ; x++) {
			for (int y = 0 ; y < botIP.getHeight() ; y++) {
				int surPix = botIP.getPixel(x, y);
				volumeBelowColumn += surPix;
			}
		}
		
		//calculate volume of complete column
		int columnVolume = 0;
		for (int i = 0 ; i < jCO.heightOfColumn ; i++) {
			columnVolume += Math.round(Math.PI * jCO.innerMajorRadius[i] * jCO.innerMinorRadius[i]);						
		}
		
		bulkVolume[0] = columnVolume - volumeAboveColumn - volumeBelowColumn;
		bulkVolume[1] = volumeAboveColumn;
		bulkVolume[2] = volumeBelowColumn;
		
		return bulkVolume;
		
	}

	public ArrayList<Integer> check4TouchingTheBottom(int imageNumber, InputOutput.MyFolderCollection mFC, ImagePlus nowTiff, PoreClusterProperties myP, boolean isRealColumn) {
	
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

	public ArrayList<Integer> check4TouchingTheTop(int imageNumber, InputOutput.MyFolderCollection mFC, ImagePlus nowTiff, PoreClusterProperties myP, boolean isRealColumn) {
	
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

	public PoreClusterProperties checkForPercolatingClusters(int imageNumber, InputOutput.MyFolderCollection mFC, ImagePlus nowTiff, PoreClusterProperties myP, boolean isRealColumn) {
			
			int i;		
			
			//init list for collecting connected custers
			ArrayList<Integer> conTop = check4TouchingTheTop(imageNumber, mFC, nowTiff, myP, isRealColumn);	
			ArrayList<Integer> conBot = check4TouchingTheBottom(imageNumber, mFC, nowTiff, myP, isRealColumn);	
							
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
			int numOfCluster = myP.id.length;
			boolean[] isPercolating = new boolean[numOfCluster];
			boolean[] touchesTop = new boolean[numOfCluster];
			boolean[] touchesBot = new boolean[numOfCluster];
			
			//init touching information
			for (i = 0 ; i < isPercolating.length ; i++) {
				isPercolating[i] = false;
				touchesTop[i] = false;
				touchesBot[i] = false;
			}
			
			//find the location of the cluster id in the pore cluster properties..
			int[] idPosition = new int[myP.id.length+1];		
			for (i = 0 ; i < myP.id.length ; i++) idPosition[myP.id[i]] = i;
			
			//set the touches that are true
			for (int temp : conTop) touchesTop[idPosition[temp]] = true;
			for (int temp : conBot) touchesBot[idPosition[temp]] = true;
			for (int temp : conPerk) isPercolating[idPosition[temp]] = true;
			
			//write property to PoreClusterProperties
			myP.touchesBot = touchesBot;
			myP.touchesTop = touchesTop;
			myP.isPercolating = isPercolating;
			
	/*		float interval = (System.currentTimeMillis() - startTime) / 1000;		
			DecimalFormat df = new DecimalFormat();
			df.setMaximumFractionDigits(2);
			IJ.error("Checking for percolating clusters took " + df.format(interval) + " seconds ..");*/
			
			return myP;		
		}

	public double[] calculateFractalDimension(ImagePlus nowTiff, InputOutput.MyFolderCollection mFC) {
				 
		JFractalCount myFBC = new JFractalCount();		
		double[] myFractalResults = myFBC.count(nowTiff);
		
		return myFractalResults;
		
	}

	public void tailoredPoreSpaceAnalyses(int imageNumber, InputOutput.MyFolderCollection mFC, ImagePlus[] nowTiffs, MenuWaiter.PoreSpaceAnalyserOptions mPSA) {
				
		InputOutput jIO = new InputOutput();		
		ImageManipulator jIM = new ImageManipulator();
		PoreClusterProperties myP = new PoreClusterProperties();
		AndAllTheRest aa = new AndAllTheRest();
		
		//init some very basic variables..		
		String nowImageName = mFC.myTiffs[imageNumber].substring(0, mFC.myTiffs[imageNumber].length() - 4);
		
		//create output folders and save statistics		
		IJ.showStatus("Saving the pore-cluster properties ...");
		String myOutFolder = "Stats";
		String myOutPath = mFC.myPreOutFolder + "\\" + myOutFolder;
		new File(myOutPath).mkdir();		
		String outPath = mFC.myPreOutFolder + "\\Stats\\" + nowImageName + "_BJResults.txt";
				
		//input parameters for particle analyzer		
		int threshold = 128; //input parameter 4 particle analyzer
		int minVol = 0; //input parameter 4 particle analyzer  ..364 corresponds to 100,000,000 cubic micrometers or 0.1 cubic millimeter
		int maxVol = nowTiffs[0].getWidth() * nowTiffs[0].getHeight() * nowTiffs[0].getNSlices(); //input parameter 4 particle analyzer
		boolean excludeObjectsTouchingTheEdges = false;		
		boolean redirect2Output = true;
		
		//check is it is a real column with non-plane top and bottom surface
		boolean isRealColumn = false;
		if (nowTiffs[1]!=null) isRealColumn = true;
		
		//check if the clusters need to be identified
		myP.containsAParticleAnalysis = false;
		if (mPSA.performParticleAnalyses == true) {
		
			myP.containsAParticleAnalysis = true;
			
			// get the pore clusters	
			IJ.showStatus("Identifying connected pore-clusters ...");
			IJ.freeMemory();IJ.freeMemory();
			Counter3D jOC = new Counter3D(nowTiffs[0], 128, minVol, maxVol, excludeObjectsTouchingTheEdges, redirect2Output);
			jOC.getObjects();			
			
			//init output variables	
			int numOfObjects = jOC.getObjMapAsArray().length;	
			int[] id = jOC.getObjMapAsArray();								//ID: unique particle identifier; this number is the label used for the object in all calculations and output
			float[] volume = new float[id.length];							//Vol: object volume
			float[][] centersOfMass = new float[id.length][3];				//coordinates of CenterOfMass			
			float[][] centroids = new float[id.length][3];					//coordinates of centroid
			int[][] topLeftCornerOfBoundingCube = new int[id.length][3];	//coordinates of top left corner of bounding cube around object
			int[][] bottomRightCornerOfBoundingCube = new int[id.length][3];//coordinates of bottom right corner of bounding cube around object
			int[] numberOfSurfaceVoxels = new int[id.length];
			float[] meanDistance2Surface = new float[id.length];
			float[] medianDistance2Surface= new float[id.length];
			float[] standardDeviationOfDistance2Surface = new float[id.length];
						
			//loop through each pore cluster and get its properties
			IJ.freeMemory();IJ.freeMemory();
			for (int i = 0 ; i < numOfObjects ; i++) {
				
				Object3D nowObject = jOC.getObject(i);  
				
				volume[i] = nowObject.size;
				
				float[] centerOfMass = nowObject.c_mass;
				for (int j = 0 ; j < 3 ; j++) centersOfMass[i][j] = centerOfMass[j];
				
				float[] centroid = nowObject.centroid;
				for (int j = 0 ; j < 3 ; j++) centroids[i][j] = centroid[j];
				
				int[] topLeftcorner = nowObject.bound_cube_TL;
				for (int j = 0 ; j < 3 ; j++) topLeftCornerOfBoundingCube[i][j] = topLeftcorner[j];
				
				int[] bottomRightcorner = nowObject.bound_cube_BR;
				for (int j = 0 ; j < 3 ; j++) bottomRightCornerOfBoundingCube[i][j] = bottomRightcorner[j];				
		
				numberOfSurfaceVoxels[i] = nowObject.surf_size;
				
				/*MCTriangulator myTriangulator = new MCTriangulator();
				int resamplingF = 6;
				boolean[] channelChoice = {true, false, false};
				Vector<int[]> surfaceVoxels = nowObject.surf_voxelsCoord;
				List<Point3f> myTriangles = myTriangulator.getTriangles(nowTiffs[0], threshold, channelChoice, resamplingF);*/
								
				meanDistance2Surface[i] = nowObject.mean_dist2surf;
				medianDistance2Surface[i]= nowObject.median_dist2surf;
				standardDeviationOfDistance2Surface[i] = nowObject.SD_dist2surf;
				
				
			}
			
		/*	//calculate thicknesses
			double[][] thick = new double[nParticles][2];		
			Local_Thickness_Parallel th = new Local_Thickness_Parallel();
			ImagePlus thickImp = null;
			if (mPSA.calcThickness == true) {			
				
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
				thick = jJPA.getMeanStdDev(thickImp, particleLabels, particleSizes, 0);			
				for (int i = 1; i < nParticles; i++) max = Math.max(max, thick[i][2]);
			}
			
			//moments, unit vectors and orientation 
			EigenvalueDecomposition[] eigens = new EigenvalueDecomposition[nParticles];		
			if (mPSA.calcMoments == true | mPSA.calcUnitVectors == true)  eigens = jJPA.getEigens(nowTiffs[0], particleLabels, centroids);	
			
			// init variables for export of results
		
				
			//sort results 
			int[] sortedIndecies = aa.getIndicesInOrder(volumes);
			
			//write results into a biiiiig table		
			for (int i = 1; i < volumes.length; i++) {
				if (volumes[i] > 0) {			
					
					IJ.showStatus("Sorting out properties of the individual pore-clusters ...");
					
					id[i-1] = sortedIndecies[i];				
					volume[i-1] = volumes[sortedIndecies[i]];
					xCenter[i-1] = centroids[sortedIndecies[i]][0];				
					yCenter[i-1] = centroids[sortedIndecies[i]][1];
					zCenter[i-1] = centroids[sortedIndecies[i]][2];
				
					xMin[i-1] = (int)Math.floor(limits[sortedIndecies[i]][0]);				
					yMin[i-1] = (int)Math.floor(limits[sortedIndecies[i]][2]);
					zMin[i-1] = (int)Math.floor(limits[sortedIndecies[i]][4]);
					xMax[i-1] = (int)Math.floor(limits[sortedIndecies[i]][1]);				
					yMax[i-1] = (int)Math.floor(limits[sortedIndecies[i]][3]);
					zMax[i-1] = (int)Math.floor(limits[sortedIndecies[i]][5]);
					
					if (mPSA.calcSurface == true) surfaceArea[i-1] = surfaceAreas[sortedIndecies[i]];				
									
					EigenvalueDecomposition E = eigens[sortedIndecies[i]];
					
					if (mPSA.calcMoments == true) { 
						momentOfInertiaShortestAxis[i-1] = E.getD().get(2, 2);
						momentOfInertiamiddleAxis[i-1] = E.getD().get(1, 1);
						momentOfInertiaLongestAxis[i-1] = E.getD().get(0, 0);
					}
					
					if (mPSA.calcUnitVectors == true) {
						unitVectorInXDirection[i-1] = E.getV().get(0, 0);				
						unitVectorInYDirection[i-1] = E.getV().get(1, 0);
						unitVectorInZDirection[i-1] = E.getV().get(2, 0);
					}
					
					if (mPSA.calcEuler == true) {
						euler[i-1] = eulerCharacters[sortedIndecies[i]][0];				
						holes[i-1] = eulerCharacters[sortedIndecies[i]][1];
						cavities[i-1] = eulerCharacters[sortedIndecies[i]][2];
					}
					
					if (mPSA.calcThickness == true) {
						thickness[i-1] = thick[sortedIndecies[i]][0];				
						sdThickness[i-1] = thick[sortedIndecies[i]][1];
						maxThickness[i-1] = thick[sortedIndecies[i]][2];
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
				
			//calculate percolating clusters ... needs the other parameters.. and is needed for following parameters
			ImagePlus particleLabelTiff = jJPA.displayParticleLabels(particleLabels, nowTiffs[0]);
			myP = checkForPercolatingClusters(imageNumber, mFC, particleLabelTiff, myP, isRealColumn);		
			
			IJ.freeMemory();IJ.freeMemory();
			
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
			myP.globalConnectionProbability = sumOfSquaredClusterSizes / (myP.globVolume / rF * myP.globVolume / rF);
			
			//save desired images
			if (mPSA.plotLabels == true) {
				myOutFolder = "ClusterLabels";
				myOutPath = mFC.myPreOutFolder + "\\" + myOutFolder;
				new File(myOutPath).mkdir();
				ImagePlus outTiff = particleLabelTiff;
				String clusterLabelImageName = nowImageName + "_ClusterLables.tif";
				jIO.tiffSaver(myOutPath, clusterLabelImageName, outTiff);
				IJ.freeMemory();IJ.freeMemory();
			}
			
			if (mPSA.plotThickness == true) {			
				myOutFolder = "PoreThick";
				myOutPath = mFC.myPreOutFolder + "\\" + myOutFolder;
				new File(myOutPath).mkdir();
				String nowDir = mFC.myPreOutFolder + "\\PoreThick";		
				String thicknessImageName = nowImageName + "_Thickness.tif";
				jIO.tiffSaver(nowDir, thicknessImageName, thickImp);	
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
				IJ.freeMemory();IJ.freeMemory();
			}*/
			
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
		myP.globVolume = macroPoreVolume(nowTiffs[0]);
		
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
		public double globVolume;
		public double globSurface;
		public double globThickness;
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

