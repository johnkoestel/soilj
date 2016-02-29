package SoilJ.tools;

import java.awt.Color;

import org.apache.commons.math3.stat.StatUtils;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.WaitForUserDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class DisplayThings implements PlugIn {

	
	public void run(String arg) {
				//ok, this is not needed..
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
	
	public void plotXYXY(double[] x, double[] y, double[] x2, double[] y2,String title, String xLabel, String yLabel) {
		
		RollerCaster rC = new RollerCaster();
			
		//find maxima and minima
		double mini=StatUtils.min(y);
		double maxi=StatUtils.max(y);
		
		Plot rmo = new Plot(title, xLabel, yLabel);	
		
		rmo.setLimits(StatUtils.min(x), StatUtils.max(x), mini-3, maxi+3);
		rmo.setColor(Color.BLUE);
		rmo.addPoints(rC.castDouble2Float(x), rC.castDouble2Float(y), Plot.CIRCLE);
		rmo.setColor(Color.RED);
		rmo.addPoints(rC.castDouble2Float(x2), rC.castDouble2Float(y2), Plot.CIRCLE);
		rmo.draw();rmo.show();		
		
		GenericDialog gd = new GenericDialog("");
		gd.showDialog();
	    if (gd.wasCanceled());
	    else new WaitForUserDialog("press any key");
	    gd.removeAll();
		
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
