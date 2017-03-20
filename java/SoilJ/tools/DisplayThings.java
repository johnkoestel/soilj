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

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;

import org.apache.commons.math3.stat.StatUtils;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.FreehandRoi;
import ij.gui.GenericDialog;
import ij.gui.Line;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.PointRoi;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.plugin.ContrastEnhancer;
import ij.plugin.PlugIn;
import ij.process.AutoThresholder;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

/** 
 * DisplayThings is a SoilJ class to display things (no, really!). It is mostly used for debugging purposes when developing SoilJ.
 * 
 * @author John Koestel
 *
 */

public class DisplayThings implements PlugIn {

	
	public void run(String arg) {
				//ok, this is not needed..
	}

	public void showMeMyRoi(String name, ImageProcessor myIP, PolygonRoi pRoi, int colorCode) {
	
		ImagePlus nowTiff = new ImagePlus(name, myIP);
		
		if (pRoi != null) {
			Overlay myO = new Overlay(pRoi);		
			myO.setStrokeColor(Color.YELLOW); // if 1
			if (colorCode == 0) myO.setStrokeColor(Color.WHITE);		
			if (colorCode == 2) myO.setStrokeColor(Color.BLUE);
			if (colorCode == 3) myO.setStrokeColor(Color.RED);		
			nowTiff.setOverlay(myO);
		} else nowTiff.setHideOverlay(true);			
		
		nowTiff.setTitle(name);
		nowTiff.updateAndDraw();
		nowTiff.show();
		
		IJ.wait(2500);
		
		//nowTiff.hide();	
	}
	
	public void displayColumnOutlinesByZ(ImagePlus nowTiff, ObjectDetector.ColCoords3D colCoords, InputOutput.MyFolderCollection mFC) {
		
		TailoredMaths math = new TailoredMaths();
		InputOutput jIO = new InputOutput();
		
		double[] myAngels = {0, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4};
		double[] allAngels = {0, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4, Math.PI, 5 * Math.PI / 4, 3 * Math.PI / 2, 7 * Math.PI / 4};
		int radius = (nowTiff.getWidth() + nowTiff.getHeight()) / 4;
		double midX = nowTiff.getWidth() / 2;
		double midY = nowTiff.getHeight() / 2;
		
		float[] x = new float[8];
		float[] y = new float[8];
		float[] xo = new float[8];
		float[] yo = new float[8];
		
		//init output images
		ImageProcessor[] outIP = new ImageProcessor[4];
		ImageProcessor standardIP = new ShortProcessor(2 * radius - 1, nowTiff.getNSlices());
		for (int i = 0 ; i < outIP.length ; i++) outIP[i] = standardIP.duplicate();
		
		float[][][] inner = new float[colCoords.bottomOfColumn - colCoords.topOfColumn][4][2];
		float[][][] outer = new float[colCoords.bottomOfColumn - colCoords.topOfColumn][4][2];
		
		for (int z = 0 ; z < nowTiff.getNSlices() ; z++) {
			
			nowTiff.setPosition(z + 1);
			ImageProcessor nowIP = nowTiff.getProcessor();
			
			IJ.showStatus("Preparing visualization of column outlines for image-slice " + z + " / " + nowTiff.getNSlices());
					
			//if z within the detected column
			if (z >= colCoords.topOfColumn & z < colCoords.bottomOfColumn) {
			
				double[][] xy = new double[8][8];
				
				//inner edge
				int i = z - colCoords.topOfColumn;
				double majRad = colCoords.innerMajorRadius[i];
				double minRad = colCoords.innerMinorRadius[i];
				
				xy = math.getXYOfEllipseFromAngle(allAngels, colCoords.ixmid[i], colCoords.iymid[i], majRad, minRad, colCoords.itheta[i]);
				
				for (int j = 0 ; j < allAngels.length ; j++) {
					x[j] = (float)xy[j][0];			
					y[j] = (float)xy[j][1];
				}
				
				//outer edge
				i = z - colCoords.topOfColumn;
				majRad = colCoords.outerMajorRadius[i];
				minRad = colCoords.outerMinorRadius[i];
				
				xy = math.getXYOfEllipseFromAngle(allAngels, colCoords.xmid[i], colCoords.ymid[i], majRad, minRad, colCoords.theta[i]);
				
				for (int j = 0 ; j < allAngels.length ; j++) {
					xo[j] = (float)xy[j][0];
					yo[j] = (float)xy[j][1];
				}				
			}
			
			for (int angle = 0 ; angle < myAngels.length ; angle++) {
				
				double alpha = myAngels[angle];
				float[][] nowPX = new float[2 * radius - 1][4];
				float[][] nowPY = new float[2 * radius - 1][4];
				
				//get profile				
				int cc = 0;				
				double[] minDist = {radius, radius};
				double[] edgeCoord = new double[2];
				double[] monDist = {radius, radius};
				double[] odgeCoord = new double[2];
				
				for (int r = -radius + 1; r < radius ; r++) {
					
					nowPX[cc][angle] = (float)(midX + r * Math.sin(alpha));
					nowPY[cc][angle] = (float)(midY + r * Math.cos(alpha));					
					
					//if z within the detected column
					if (z >= colCoords.topOfColumn & z < colCoords.bottomOfColumn) {
					
						//check whether one of the inner edges in crossed by the profile
						double distance1 = Math.sqrt(Math.pow(nowPX[cc][angle] - x[angle],2) + Math.pow(nowPY[cc][angle] - y[angle],2));
						double distance2 = Math.sqrt(Math.pow(nowPX[cc][angle] - x[angle + 4],2) + Math.pow(nowPY[cc][angle] - y[angle + 4],2));
						
						if (distance1 < minDist[0]) {
							minDist[0] = distance1;
							edgeCoord[0] = cc;
						}
						if (distance2 < minDist[1]) {
							minDist[1] = distance2;
							edgeCoord[1] = cc;
						}
						
						//check whether one of the outer edges in crossed by the profile
						double oistance1 = Math.sqrt(Math.pow(nowPX[cc][angle] - xo[angle],2) + Math.pow(nowPY[cc][angle] - yo[angle],2));
						double oistance2 = Math.sqrt(Math.pow(nowPX[cc][angle] - xo[angle + 4],2) + Math.pow(nowPY[cc][angle] - yo[angle + 4],2));
						
						if (oistance1 < monDist[0]) {
							monDist[0] = oistance1;
							odgeCoord[0] = cc;
						}
						if (oistance2 < monDist[1]) {
							monDist[1] = oistance2;
							odgeCoord[1] = cc;
						}
					}
					
					cc++;				
					
				}
				
				//transfer detected edges to inner and outer
				if (z >= colCoords.topOfColumn & z < colCoords.bottomOfColumn) {
					inner[z - colCoords.topOfColumn][angle][0] = (float)edgeCoord[0];				
					inner[z - colCoords.topOfColumn][angle][1] = (float)edgeCoord[1];
					outer[z - colCoords.topOfColumn][angle][0] = (float)odgeCoord[0];
					outer[z - colCoords.topOfColumn][angle][1] = (float)odgeCoord[1];
				}
				
				//remember profile				
				for (int i = 0 ; i < nowPX.length ; i++) {
					int nowPixelValue = nowIP.getPixel((int)nowPX[i][angle], (int)nowPY[i][angle]);
					outIP[angle].putPixel(i, z, nowPixelValue);
				}
			}
		}
				
		ImageStack outStack = new ImageStack(outIP[0].getWidth(), outIP[0].getHeight());
		ImagePlus outTiff = new ImagePlus();
		
		for (int j = 0 ; j < myAngels.length ; j++) {
			
			//construct inner overlay
			float[] Z = new float[2 * (colCoords.bottomOfColumn - colCoords.topOfColumn)];
			for (int i = 0 ; i < colCoords.bottomOfColumn - colCoords.topOfColumn ; i++) Z[i] = i + colCoords.topOfColumn;
			for (int i = 0 ; i < colCoords.bottomOfColumn - colCoords.topOfColumn ; i++) Z[i + colCoords.bottomOfColumn - colCoords.topOfColumn] = colCoords.bottomOfColumn - i - 1;
			float[] inin = new float[Z.length];
			for (int i = 0 ; i < Z.length / 2 ; i++) inin[i] = inner[i][j][0];
			for (int i = Z.length / 2 ; i < Z.length ; i++) {
				int ii = i - Z.length / 2;
				int endNumber = colCoords.bottomOfColumn - colCoords.topOfColumn - 1;
				inin[i] = inner[endNumber - ii][j][1];			
			}
			
			//construct outer Overlays
			float[] ZO = new float[colCoords.bottomOfColumn - colCoords.topOfColumn];
			float[] l1 = new float[ZO.length];
			float[] l2 = new float[ZO.length];
			for (int i = 0 ; i < ZO.length ; i++) {
				ZO[i] = i + colCoords.topOfColumn;
				l1[i] = outer[i][j][0];
				l2[i] = outer[i][j][1];
			}
			
			ContrastEnhancer myCE = new ContrastEnhancer();
			myCE.stretchHistogram(outIP[j], 0.5);
			
			ImageProcessor zIP = outIP[j].duplicate();
			ImageProcessor rgbIP = zIP.convertToRGB();
			rgbIP.setColor(Color.RED);
			
			PolygonRoi pRoi = new PolygonRoi(inin, Z, Roi.POLYGON);		
			PolygonRoi oRoi1 = new PolygonRoi(l1, ZO, Roi.FREELINE);
			PolygonRoi oRoi2 = new PolygonRoi(l2, ZO, Roi.FREELINE);
			
			rgbIP.setLineWidth(2);
			pRoi.drawPixels(rgbIP);
			oRoi1.drawPixels(rgbIP);
			oRoi2.drawPixels(rgbIP);
			
			outStack.addSlice(rgbIP);
		}
		
		outTiff.setStack(outStack);
				
		ContrastEnhancer myCE = new ContrastEnhancer();
		myCE.stretchHistogram(outTiff, 0.5);
		
		//save the results
		String myOutFolderName = mFC.myOutFolder + "\\ReviewFoundOutlines";
		new File(myOutFolderName).mkdir();
		jIO.tiffSaver(myOutFolderName, mFC.fileName, outTiff);
		
	}
	
	
	public ImagePlus jDisplayParticleLabels(int[][] particleLabels, ImagePlus imp) {  //this object was written by Doube.. but private in BoneJ.. therefore the copy and past.
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getImageStackSize();
		final int wh = w * h;
		ImageStack stack = new ImageStack(w, h);
		double max = 0;
		for (int z = 0; z < d; z++) {
			float[] slicePixels = new float[wh];
			for (int i = 0; i < wh; i++) {
				slicePixels[i] = (float) particleLabels[z][i];
				max = Math.max(max, slicePixels[i]);
			}
			stack.addSlice(imp.getImageStack().getSliceLabel(z + 1),
					slicePixels);
		}
		ImagePlus impParticles = new ImagePlus(imp.getShortTitle() + "_parts",
				stack);
		impParticles.setCalibration(imp.getCalibration());
		impParticles.getProcessor().setMinAndMax(0, max);
		return impParticles;
	}
	
	public void displayIP(ImageProcessor myIP, String label) {		
		
		ImagePlus test = new ImagePlus(label, myIP);		
		test.draw();test.show();
		
	}

	public void plotDepthProfiles(ObjectDetector.RadialModes myRM, int[] donuts2Plot) {
		
		double[] z = new double[myRM.maskingThreshold.length];
		double[] myData1 = new double[myRM.maskingThreshold.length];
		double[] myData2 = new double[myRM.maskingThreshold.length];
				
		for (int i = 0 ; i < myRM.maskingThreshold.length ; i++) {
					
			z[i] = myRM.maskingThreshold.length - i - 1;	
			myData1[i] = myRM.maskedRadialMinima[i][donuts2Plot[0]];	
			myData2[i] = myRM.maskedRadialMinima[i][donuts2Plot[1]];			
		}
		
		//find maxima and minima
		double[] mini = new double[2];mini[0]=StatUtils.min(myData1);mini[1]=StatUtils.min(myData2);
		double[] maxi = new double[2];maxi[0]=StatUtils.max(myData1);maxi[1]=StatUtils.max(myData2);
		
		Plot rmo = new Plot(donuts2Plot[0] + " and " + donuts2Plot[1],"brightness", "depth");
		
		rmo.setLimits(StatUtils.min(mini)-50, StatUtils.max(maxi)+50, 0, myRM.maskingThreshold.length);
		rmo.setColor(Color.BLUE);
		rmo.addPoints(myData1, z, Plot.LINE);
		rmo.setColor(Color.RED);
		rmo.addPoints(myData2, z, Plot.LINE);		
		rmo.draw();rmo.show();
		
	}
	
	public void plotVerticalProfile(double[] thing2Plot, String title, String xLabel) {
		
		double[] z = new double[thing2Plot.length];
		
		for (int i = 0 ; i < thing2Plot.length ; i++) {
					
			z[i] = thing2Plot.length - i - 1;	
			
		}
		
		//find maxima and minima
		double mini = StatUtils.min(thing2Plot);
		double maxi = StatUtils.max(thing2Plot);
		double std = Math.sqrt(StatUtils.variance(thing2Plot));
		
		Plot rmo = new Plot(title, xLabel, "depth");
		
		rmo.setLimits(mini - 0.1 * std, maxi + 0.1 * std, 0, z.length);
		rmo.setColor(Color.BLUE);
		rmo.addPoints(thing2Plot, z, Plot.LINE);
		rmo.draw();rmo.show();
		
	}

	public void plotRadialProfiles(ObjectDetector.RadialModes myRM, int[] depth2Plot) {
		
		double[] myData1 = new double[myRM.radius.length];
		double[] myData2 = new double[myRM.radius.length];
		
		for (int i = 0 ; i < myRM.radius.length ; i++) {	
			myData1[i] = myRM.maskedRadialMinima[depth2Plot[0]][i];	 
			myData2[i] = myRM.maskedRadialMinima[depth2Plot[1]][i];
		}
		
		//find maxima and minima
		double[] mini = new double[2];mini[0]=StatUtils.min(myData1);mini[1]=StatUtils.min(myData2);
		double[] maxi = new double[2];maxi[0]=StatUtils.max(myData1);maxi[1]=StatUtils.max(myData2);
		
		Plot rmo = new Plot(depth2Plot[0] + " and " + depth2Plot[1], "radial distance", "brightness");	
		
		rmo.setLimits(0, StatUtils.max(myRM.radius), StatUtils.min(mini)-50, StatUtils.max(maxi)+50);
		rmo.setColor(Color.BLUE);
		rmo.addPoints(myRM.radius, myData1, Plot.CIRCLE);
		rmo.setColor(Color.RED);
		rmo.addPoints(myRM.radius, myData2, Plot.CIRCLE);		
		rmo.draw();rmo.show();		
	}
	
	public boolean plotSpecialWallFinderRadii(double[] x, double[] y, double[] ygrad, double[] x2) {
		

		RollerCaster rC = new RollerCaster();
		
		///////////////////////////////////////////////////
		// plot grey values
		////////////////////////////////////////////////////
		
		
		//find maxima and minima
		double mini=StatUtils.min(y);				;
		double maxi=StatUtils.max(y);
		double minmax = maxi-mini;
		
		double[] ylim = {mini-0.05*minmax, maxi+0.05*minmax}; 
		double[] outerWallPerimeter = {x2[0], x2[0]};
		double[] innerWallPerimeter = {x2[1], x2[1]};
		
		
		Plot rgv = new Plot("median filtered grey values", "radial coordinate (vx)", "grey value");	
		
		rgv.setLimits(StatUtils.min(x), StatUtils.max(x), ylim[0], ylim[1]);
		rgv.setColor(Color.BLUE);
		rgv.addPoints(rC.castDouble2Float(x), rC.castDouble2Float(y), Plot.LINE);
		rgv.setColor(Color.RED);
		rgv.addPoints(rC.castDouble2Float(outerWallPerimeter), rC.castDouble2Float(ylim), Plot.LINE);
		rgv.addPoints(rC.castDouble2Float(innerWallPerimeter), rC.castDouble2Float(ylim), Plot.LINE);
		rgv.draw();
		PlotWindow win1 = rgv.show();	

		//find position for second plot
		int newPosX = win1.getX();
		int newPosY = win1.getY() + win1.getHeight();
		
		///////////////////////////////////////////////////
		//plot Gradient
		////////////////////////////////////////////////////
		
		//find maxima and minima
		mini=StatUtils.min(ygrad);				;
		maxi=StatUtils.max(ygrad);
		minmax = maxi-mini;
		
		ylim[0] = mini-0.05*minmax;
		ylim[1] = maxi+0.05*minmax; 
		
		Plot rmo = new Plot("radial gradients of median filtered grey values", "radial coordinate (vx)", "gradient");	
		
		rmo.setLimits(StatUtils.min(x), StatUtils.max(x), mini-3, maxi+3);
		rmo.setColor(Color.GREEN);
		rmo.addPoints(rC.castDouble2Float(x), rC.castDouble2Float(ygrad), Plot.LINE);
		rmo.setColor(Color.RED);
		rmo.addPoints(rC.castDouble2Float(outerWallPerimeter), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.addPoints(rC.castDouble2Float(innerWallPerimeter), rC.castDouble2Float(ylim), Plot.LINE);		
		rmo.draw();
		
		PlotWindow win2 = rmo.show();	
		win2.setLocation(newPosX, newPosY);
				
		GenericDialog gd = new GenericDialog("");
		gd.showDialog();
	    if (gd.wasCanceled()) return false;
	    else new WaitForUserDialog("press any key");
	    gd.removeAll();
	    	    
	    win1.close();
	    win2.close();
	    rgv.dispose();
	    rmo.dispose();
	    
	    IJ.freeMemory();IJ.freeMemory();
	    
	    return true;
		
	}
	
	
	public boolean plotHistogramWithThresholds(double[] x, double[] y, double[] threshes, String title, String xLabel, String yLabel) {
		
		RollerCaster rC = new RollerCaster();		
		
		//find maxima and minima
		double mini=StatUtils.min(y);
		double maxi=StatUtils.max(y);
		
		Plot rmo = new Plot(title, xLabel, yLabel);	
		
		rmo.setLimits(StatUtils.min(x), StatUtils.max(x), mini-3, maxi+3);
		rmo.setColor(Color.BLACK);
		rmo.addPoints(rC.castDouble2Float(x), rC.castDouble2Float(y), Plot.LINE);
				
		//add thresholds
		double[] ylim = {mini, maxi};
		double[] meth1 = {threshes[0], threshes[0]};
		double[] meth2 = {threshes[1], threshes[1]};
		double[] meth3 = {threshes[2], threshes[2]};
		double[] meth4 = {threshes[3], threshes[3]};
		double[] meth5 = {threshes[4], threshes[4]};
		double[] meth6 = {threshes[5], threshes[5]};
		double[] meth7 = {threshes[6], threshes[6]};
		double[] meth8 = {threshes[7], threshes[7]};
		double[] meth9 = {threshes[8], threshes[8]};
		
		rmo.setColor(Color.RED); //default (also some isodata)
		rmo.addPoints(rC.castDouble2Float(meth1), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.setColor(Color.BLUE); // Otsu			
		rmo.addPoints(rC.castDouble2Float(meth2), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.setColor(Color.GREEN); //Huang
		rmo.addPoints(rC.castDouble2Float(meth3), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.setColor(Color.CYAN);  //max entropy			
		rmo.addPoints(rC.castDouble2Float(meth4), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.setColor(Color.MAGENTA);  //minimum
		rmo.addPoints(rC.castDouble2Float(meth5), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.setColor(Color.ORANGE);   //minError
		rmo.addPoints(rC.castDouble2Float(meth6), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.setColor(Color.YELLOW);   //Renyi Entropy
		rmo.addPoints(rC.castDouble2Float(meth7), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.setColor(Color.PINK);   ///triangle
		rmo.addPoints(rC.castDouble2Float(meth8), rC.castDouble2Float(ylim), Plot.LINE);
		rmo.setColor(Color.GRAY);   //IsoData
		rmo.addPoints(rC.castDouble2Float(meth9), rC.castDouble2Float(ylim), Plot.LINE);		
	
		//show histogram and thresholds..
		rmo.draw();
		
		PlotWindow win1 = rmo.show();
		GenericDialog gd = new GenericDialog("");
		gd.showDialog();
	    if (gd.wasCanceled()) return false;
	    else new WaitForUserDialog("press any key");
	    gd.removeAll();
	    	    
	    win1.close();	    
	    rmo.dispose();
	    
	    IJ.freeMemory();IJ.freeMemory();
	    
	    IJ.freeMemory();IJ.freeMemory();
	    
	    return true;
		
	}
	
	public boolean plotXYXY(double[] x, double[] y, double[] x2, double[] y2,String title, String xLabel, String yLabel) {
		
		RollerCaster rC = new RollerCaster();	
		
		//find maxima and minima
		//find maxima and minima
		double[] mini = new double[2];mini[0]=StatUtils.min(y);mini[1]=StatUtils.min(y2);
		double[] maxi = new double[2];maxi[0]=StatUtils.max(y);maxi[1]=StatUtils.max(y2);
		
		Plot rmo = new Plot(title, xLabel, yLabel);		
		
		rmo.setLimits(StatUtils.min(x), StatUtils.max(x), StatUtils.min(mini) - 3, StatUtils.max(maxi) + 3);
		rmo.setColor(Color.BLUE);
		rmo.addPoints(rC.castDouble2Float(x), rC.castDouble2Float(y), Plot.LINE);
		rmo.setColor(Color.RED);
		rmo.addPoints(rC.castDouble2Float(x2), rC.castDouble2Float(y2), Plot.LINE);
		rmo.draw();
		PlotWindow pW = rmo.show();		 
		
		GenericDialog gd = new GenericDialog("");
		gd.showDialog();
	    if (gd.wasCanceled()) return false;
	    else new WaitForUserDialog("press any key");
	    gd.removeAll();
	    
	    pW.removeAll();
	    pW.dispose();
	    rmo.dispose();
	    
	    IJ.freeMemory();IJ.freeMemory();
	    
	    return true;
		
	}

	public void plotRadialProfilesWithSpecialFunctionFit(ObjectDetector.RadialModes myRM, FitStuff.FittingResults myResults, int depth2Plot) {
		
		FitStuff fittings = new FitStuff();
		
		FitStuff.HyperbolicFunction2Params mLF = fittings.new HyperbolicFunction2Params(); 
		double[] myData = new double[myRM.radius.length];
		double[] myParams = new double[myResults.numberOfParams];
		double[] realRadius = new double[(int)StatUtils.max(myRM.radius)];
		double[] myFit = new double[realRadius.length];
		
		mLF.setup(StatUtils.max(myRM.radius), myRM.maskedRadialMinima[depth2Plot][myRM.radius.length - 1]);
		
		for (int i = 0 ; i < myRM.radius.length ; i++) myData[i] = myRM.maskedRadialMinima[depth2Plot][i];
		for (int i = 0 ; i < myResults.numberOfParams ; i++) myParams[i] = myResults.params[depth2Plot][i];		
		for (int i = 0 ; i < myFit.length ; i++) {		
			realRadius[i] = i;
			myFit[i] = mLF.value(realRadius[i], myParams);
		}
		
		//find maxima and minima
		double mini = StatUtils.min(myData);
		double maxi = StatUtils.max(myData);
		
		Plot rmo = new Plot("" + depth2Plot, "radial distance", "brightness");	
		
		rmo.setLimits(0, StatUtils.max(myRM.radius), mini - 50, maxi + 50);
		rmo.setColor(Color.BLUE);
		rmo.addPoints(myRM.radius, myData, Plot.CIRCLE);
		rmo.setColor(Color.RED);
		rmo.addPoints(realRadius, myFit, Plot.LINE);		
		rmo.draw();rmo.show();		
	}
}

