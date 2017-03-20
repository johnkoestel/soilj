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

import ij.plugin.PlugIn;

/** 
 * REVAnalyses is a SoilJ class containing subroutines for investigations of representative
 * elementary volumes (REVs) for image properties.
 * 
 * @author John Koestel
 *
 */

public class REVAnalyses implements PlugIn {

	public void run(String arg) {
		//ok, this is not needed..
	}
	
	public class REVAnalysesPack {
		
		public String typeOfROI;
		
		public int[] x1;
		public int[] x2;
		public int[] y1;
		public int[] y2;
		public int[] z1;
		public int[] z2;
		
		public int[] radius;
		
		public String[] roiName;		
		
		public int numberOfROIs;
	}
	
	public REVAnalysesPack createROIs4REVAnalyses(MenuWaiter.REVAnalyzerOptions mRA) {
		
		REVAnalysesPack rAP = new REVAnalysesPack();
		
		int numberOfROIs = 1;
		
		rAP.typeOfROI = mRA.choiceOfRoi;		
		
		if (mRA.choiceOfMethod.equalsIgnoreCase("Sub-ROI by division")) {
			numberOfROIs = 1;
			for (int i = 1 ; i <= mRA.divNumber ; i++) {
				numberOfROIs += Math.pow(8, i);
			}
		}
			
		if (mRA.choiceOfMethod.equalsIgnoreCase("Sub-ROI by shrinkage")) {
			numberOfROIs = mRA.stepNumber + 1;
		}
		
		if (mRA.choiceOfMethod.equalsIgnoreCase("Sub-ROI by division and shrinkage")) {
			numberOfROIs = 1;
			for (int i = 1 ; i <= mRA.divNumber ; i++) {
				numberOfROIs += Math.pow(8, i);
			}
			int divROIs = numberOfROIs;
			int shrinkROIs = mRA.stepNumber * (int)Math.round(Math.pow(8, mRA.divNumber));
			numberOfROIs = divROIs + shrinkROIs;
		}
		
		int[] x1 = new int[numberOfROIs];
		int[] x2 = new int[numberOfROIs];
		int[] y1 = new int[numberOfROIs];
		int[] y2 = new int[numberOfROIs];
		int[] z1 = new int[numberOfROIs];
		int[] z2 = new int[numberOfROIs];
		
		String[] roiName = new String[numberOfROIs];
		
		if (mRA.choiceOfMethod.equalsIgnoreCase("Sub-ROI by division")) {
			
			x1[0] = mRA.cubeX1;x2[0] = mRA.cubeX2;
			y1[0] = mRA.cubeY1;y2[0] = mRA.cubeY2;
			z1[0] = mRA.cubeZ1;z2[0] = mRA.cubeZ2;
			
			roiName[0] = "X" + mRA.cubeX1 + "Y" + mRA.cubeY1 + "Z" + mRA.cubeZ1 + "EX" + mRA.edgeX + "EY" + mRA.edgeY + "EZ" + mRA.edgeZ + ".tif"; 
			
			int nowNumOfROIs = 1;
			
			for (int i = 1 ; i <= mRA.divNumber ; i++) {
				
				int nowEdgeX = mRA.edgeX / (int)Math.round(Math.pow(2, i));
				int nowEdgeY = mRA.edgeY / (int)Math.round(Math.pow(2, i));
				int nowEdgeZ = mRA.edgeZ / (int)Math.round(Math.pow(2, i));
				
				for (int x = 0 ; x < (int)Math.round(Math.pow(2, i)) ; x++) {
					for (int y = 0 ; y < (int)Math.round(Math.pow(2, i)) ; y++) {
						for (int z = 0 ; z < (int)Math.round(Math.pow(2, i)) ; z++) {
							x1[nowNumOfROIs] = (int)Math.round(mRA.cubeX1 + x*nowEdgeX); 
							x2[nowNumOfROIs] = (int)Math.round(mRA.cubeX1 + (x+1)*nowEdgeX);
							y1[nowNumOfROIs] = (int)Math.round(mRA.cubeY1 + y*nowEdgeY); 
							y2[nowNumOfROIs] = (int)Math.round(mRA.cubeY1 + (y+1)*nowEdgeY); 
							z1[nowNumOfROIs] = (int)Math.round(mRA.cubeZ1 + z*nowEdgeZ); 
							z2[nowNumOfROIs] = (int)Math.round(mRA.cubeZ1 + (z+1)*nowEdgeZ);
							
							roiName[nowNumOfROIs] = "X" + x1[nowNumOfROIs] + "Y" + y1[nowNumOfROIs] + "Z" + z1[nowNumOfROIs] + "EX" + nowEdgeX + "EY" +nowEdgeY + "EZ" + nowEdgeZ + ".tif";
							
							nowNumOfROIs++;		
						}
					}
				}
				
			}
		}
		
		if (mRA.choiceOfMethod.equalsIgnoreCase("Sub-ROI by shrinkage")) {
			
			x1[0] = mRA.cubeX1;x2[0] = mRA.cubeX2;
			y1[0] = mRA.cubeY1;y2[0] = mRA.cubeY2;
			z1[0] = mRA.cubeZ1;z2[0] = mRA.cubeZ2;
			
			roiName[0] = "X" + mRA.cubeX1 + "Y" + mRA.cubeY1 + "Z" + mRA.cubeZ1 + "EX" + mRA.edgeX + "EY" + mRA.edgeY + "EZ" + mRA.edgeZ + ".tif"; 
			
			int nowNumOfROIs = 1;
			
			for (int i = 1 ; i <= mRA.stepNumber ; i++) {
				
				int nowEdgeX = mRA.edgeX - i * mRA.stepLength;				
				int nowEdgeY = mRA.edgeY - i * mRA.stepLength;
				int nowEdgeZ = mRA.edgeZ - i * mRA.stepLength;
				
				x1[nowNumOfROIs] = (int)Math.round(mRA.cubeX1 + i * (mRA.stepLength / 2)); 
				x2[nowNumOfROIs] = (int)Math.round(x1[nowNumOfROIs] + nowEdgeX);
				y1[nowNumOfROIs] = (int)Math.round(mRA.cubeY1 + i * (mRA.stepLength / 2)); 
				y2[nowNumOfROIs] = (int)Math.round(y1[nowNumOfROIs] + nowEdgeY); 
				z1[nowNumOfROIs] = (int)Math.round(mRA.cubeZ1 + i * (mRA.stepLength / 2)); 
				z2[nowNumOfROIs] = (int)Math.round(z1[nowNumOfROIs] + nowEdgeZ);
				
				roiName[nowNumOfROIs] = "X" + x1[nowNumOfROIs] + "Y" + y1[nowNumOfROIs] + "Z" + z1[nowNumOfROIs] + "EX" + nowEdgeX + "EY" +nowEdgeY + "EZ" + nowEdgeZ + ".tif";
				
				nowNumOfROIs++;		
				
			}
		}
		
		if (mRA.choiceOfMethod.equalsIgnoreCase("Sub-ROI by division and shrinkage")) {
			
			x1[0] = mRA.cubeX1;x2[0] = mRA.cubeX2;
			y1[0] = mRA.cubeY1;y2[0] = mRA.cubeY2;
			z1[0] = mRA.cubeZ1;z2[0] = mRA.cubeZ2;
			
			roiName[0] = "X" + mRA.cubeX1 + "Y" + mRA.cubeY1 + "Z" + mRA.cubeZ1 + "EX" + mRA.edgeX + "EY" + mRA.edgeY + "EZ" + mRA.edgeZ + ".tif"; 
			
			int nowNumOfROIs = 1;
			
			int finalEdgeX = 0; 
			int finalEdgeY = 0;
			int finalEdgeZ = 0;
			
			for (int i = 1 ; i <= mRA.divNumber ; i++) {
				
				int nowEdgeX = mRA.edgeX / (int)Math.round(Math.pow(2, i));
				int nowEdgeY = mRA.edgeY / (int)Math.round(Math.pow(2, i));
				int nowEdgeZ = mRA.edgeZ / (int)Math.round(Math.pow(2, i));
				
				for (int x = 0 ; x < (int)Math.round(Math.pow(2, i)) ; x++) {
					for (int y = 0 ; y < (int)Math.round(Math.pow(2, i)) ; y++) {
						for (int z = 0 ; z < (int)Math.round(Math.pow(2, i)) ; z++) {
							x1[nowNumOfROIs] = (int)Math.round(mRA.cubeX1 + x*nowEdgeX); 
							x2[nowNumOfROIs] = (int)Math.round(mRA.cubeX1 + (x+1)*nowEdgeX);
							y1[nowNumOfROIs] = (int)Math.round(mRA.cubeY1 + y*nowEdgeY); 
							y2[nowNumOfROIs] = (int)Math.round(mRA.cubeY1 + (y+1)*nowEdgeY); 
							z1[nowNumOfROIs] = (int)Math.round(mRA.cubeZ1 + z*nowEdgeZ); 
							z2[nowNumOfROIs] = (int)Math.round(mRA.cubeZ1 + (z+1)*nowEdgeZ);
							
							roiName[nowNumOfROIs] = "X" + x1[nowNumOfROIs] + "Y" + y1[nowNumOfROIs] + "Z" + z1[nowNumOfROIs] + "EX" + nowEdgeX + "EY" +nowEdgeY + "EZ" + nowEdgeZ + ".tif";
							
							nowNumOfROIs++;		
						}
					}
				}
				
				finalEdgeX = nowEdgeX; 
				finalEdgeY = nowEdgeY;
				finalEdgeZ = nowEdgeZ;				
			}
				
			int nextNumOfROIs = nowNumOfROIs;
			int numOfSmallestDivROIs = (int)Math.round(Math.pow(8, mRA.divNumber));
			for (int k = 1 ; k <= mRA.stepNumber ; k++) {
				
				int nowEdgeX = finalEdgeX - k * mRA.stepLength;
				int nowEdgeY = finalEdgeY - k * mRA.stepLength;
				int nowEdgeZ = finalEdgeZ - k * mRA.stepLength;
				
				for (int i = numOfSmallestDivROIs ; i > 0 ; i--) {
					
					x1[nextNumOfROIs] = (int)Math.round(x1[nowNumOfROIs - i] + (k+1) * mRA.stepLength / 2); 
					x2[nextNumOfROIs] = (int)Math.round(x1[nowNumOfROIs - i] + nowEdgeX);
					y1[nextNumOfROIs] = (int)Math.round(y1[nowNumOfROIs - i] + (k+1) * mRA.stepLength / 2);
					y2[nextNumOfROIs] = (int)Math.round(y1[nowNumOfROIs - i] + nowEdgeY); 
					z1[nextNumOfROIs] = (int)Math.round(z1[nowNumOfROIs - i] + (k+1) * mRA.stepLength / 2);
					z2[nextNumOfROIs] = (int)Math.round(z1[nowNumOfROIs - i] + nowEdgeZ); 
				
					roiName[nextNumOfROIs] = "X" + x1[nextNumOfROIs] + "Y" + y1[nextNumOfROIs] + "Z" + z1[nextNumOfROIs] + 
							"EX" + nowEdgeX + "EY" +nowEdgeY + "EZ" + nowEdgeZ + ".tif";
										
					nextNumOfROIs++;
					
				}
			}
				
		}
		
		rAP.x1 = x1;
		rAP.x2 = x2;
		rAP.y1 = y1;
		rAP.y2 = y2;
		rAP.z1 = z1;
		rAP.z2 = z2;
		
		rAP.roiName = roiName;
		
		rAP.numberOfROIs = roiName.length;
		
		return rAP;
	}
	
	public MenuWaiter.PoreSpaceAnalyzerOptions castRA2PSA(MenuWaiter.REVAnalyzerOptions mRA) {
		
		MenuWaiter menu = new MenuWaiter();
		MenuWaiter.PoreSpaceAnalyzerOptions mPSA = menu.new PoreSpaceAnalyzerOptions();
		
		//set cut-out ROI option to true
		mPSA.cutCanvas = true;
		
		//link the rest of the variables
		mPSA.globVolume = mRA.globVolume;
		mPSA.globThickness = mRA.globThickness;
		mPSA.calcFractal = mRA.calcFractal;		
		mPSA.globAnisotropy = mRA.globAnisotropy;
		
		mPSA.performParticleAnalyses = mRA.performParticleAnalyses;
		
		mPSA.calcVolume = mRA.calcVolume;		
		mPSA.calcEuler = mRA.calcEuler;
		mPSA.calcThickness = mRA.calcThickness;
		mPSA.calcCriticalPoreDiameter = mRA.calcCriticalPoreDiameter;		
		mPSA.calcAnisotropy = mRA.calcAnisotropy;
		mPSA.calcPercolation = mRA.calcPercolation;		
		
		mPSA.plotLabels = mRA.plotLabels;
		mPSA.plotVolume = mRA.plotVolume;		
		mPSA.plotThickness = mRA.plotThickness;	
		mPSA.plotPercolation = mRA.plotPercolation;
		
		mPSA.choiceOfRoi = mRA.choiceOfRoi;
		
		return mPSA;		
		
	}
}

