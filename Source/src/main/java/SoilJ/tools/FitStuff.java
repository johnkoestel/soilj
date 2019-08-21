package SoilJ.tools;

import java.awt.Color;
import java.awt.Font;

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

import java.util.ArrayList;
import java.util.Random;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.measure.CurveFitter;
import ij.plugin.ContrastEnhancer;
import ij.plugin.PlugIn;
import ij.process.EllipseFitter;
import ij.process.ImageProcessor;

/** 
 * FitStuff is a SoilJ class containing subroutines for function fitting.
 * 
 * @author John Koestel
 *
 */

public class FitStuff implements PlugIn {

	
	public void run(String arg) {
				//ok, this is not needed..
	}

	public LinearFunction fitLinearFunction(double[] x, double[] y) {
		
		LinearFunction myLF = new LinearFunction();
		
		SimpleRegression regression = new SimpleRegression();
		for (int i = 0 ; i < y.length ; i++) regression.addData(x[i], y[i]);
		
		myLF.intercept = regression.getIntercept();
		myLF.slope = regression.getSlope();
		myLF.R2 = regression.getRSquare();
	
		return myLF;
		
	}

	public double evalLinearFunction(LinearFunction myLF, double x) {
		
		double data = myLF.slope * x + myLF.intercept;
		
		return data;
		
	}	


	public GoodnessOfFit calculateR2OfEllipsoidFit(double[] angleAtThisAngle, float[] xOD, float[] yOD, EllipseFitter jEF) {
		
		int j;
		GoodnessOfFit gOF = new GoodnessOfFit();
		
		//calculate R2
		float[] xfit = new float[angleAtThisAngle.length];
		float[] yfit = new float[angleAtThisAngle.length];
		float[] xhelp = new float[angleAtThisAngle.length];
		float[] yhelp = new float[angleAtThisAngle.length];		
		
		double SSX = 0; double SSY = 0;	
		double[] srx = new double[angleAtThisAngle.length];
		double[] sry = new double[angleAtThisAngle.length];
		double[] d2m = new double[angleAtThisAngle.length];  //distance to midpoint;
		
		double MX = 0; double MY = 0;			
		for (j = 0 ; j < xOD.length ; j++) MX = MX + xOD[j];
		for (j = 0 ; j < yOD.length ; j++) MY = MY + yOD[j];
		MX = MX / xOD.length;MY = MY / yOD.length;
		for (j = 0 ; j < xOD.length ; j++) SSX = SSX + Math.pow(xOD[j] - MX, 2);
		for (j = 0 ; j < yOD.length ; j++) SSY = SSY + Math.pow(yOD[j] - MY, 2);
		for (j = 0 ; j < angleAtThisAngle.length ; j++) {			
			double theta = jEF.theta;
			double alpha = angleAtThisAngle[j] - theta + Math.PI/2;
			double a = jEF.major / 2;
			double b = jEF.minor / 2;			
			xfit[j] = (float)(jEF.xCenter - a * Math.cos(alpha) * Math.cos(theta) + b * Math.sin(alpha) * Math.sin(theta));
			yfit[j] = (float)(jEF.yCenter + a * Math.cos(alpha) * Math.sin(theta) + b * Math.sin(alpha) * Math.cos(theta));
			
			srx[j] = Math.pow(xOD[j] - xfit[j], 2);
			sry[j] = Math.pow(yOD[j] - yfit[j], 2);
			
			xhelp[j] = xOD[j];
			yhelp[j] = yOD[j];		
			
			d2m[j] = Math.sqrt((jEF.xCenter - xOD[j]) * (jEF.xCenter - xOD[j]) + (jEF.yCenter - yOD[j]) * (jEF.yCenter - yOD[j]));
	
		}
		
		double R2 = (2 - StatUtils.sum(srx) / SSX - StatUtils.sum(sry) / SSY) / 2;
		
		//plot for visual inspection
		//Plot myPlot = new Plot("verification", "X","Y");
		//myPlot.setLimits(0, 1500, 0, 1500);
		//myPlot.setColor(Color.BLUE);
		//myPlot.addPoints(xhelp, yhelp, Plot.CIRCLE);
		//myPlot.setColor(Color.RED);
		//myPlot.addPoints(xfit, yfit, Plot.CROSS);
		//myPlot.draw();		
		//PlotWindow nowWindow = myPlot.show();
		//nowWindow.close();
		
		gOF.R2 = R2;
		gOF.srx = srx;
		gOF.sry = sry;
		gOF.d2m = d2m;
		gOF.median2d2m = StatUtils.percentile(d2m, 50);
		
		return gOF;
	}
	
	public FittingResult fitWithCurveFitter(String equation, double[] initialGuess, double[] xData, double[] yData) {
		
		FittingResult mFR = new FittingResult();	
		
		//do the fit		
		CurveFitter mCF = new CurveFitter(xData, yData);
		mCF.doCustomFit(equation, initialGuess, false);
		
		mFR.params = mCF.getParams();
		mFR.R2 = mCF.getFitGoodness();
		
		return mFR;
		
	}
	
	public FittingResults fitTripleExponential(ObjectDetector.RadialFacts myRF) {//define equation to be fitted
		
		FittingResults mFR = new FittingResults();
		
		int numOfLayers = myRF.imageHeight;
		double[] min = myRF.min;
		double[] max = myRF.max;
		double[][] yData = myRF.radialProfile;
		double[] x = new double[myRF.standardRadius];
		for (int i = 0 ; i < myRF.standardRadius ; i++) x[i] = i;
		
		for (int i = 0 ; i < numOfLayers ; i++) {

			Random rand = new Random();		
			double dy = max[i] - min[i];
			String equation1 = "d + e / exp(a * x - a)"; 
			String equation2 = "d + e / exp(b * x - b)"; 
			String equation3 = "d + e / exp(c * x - c)"; 
			String equation = "(" + equation1 + " + " + equation2 + " + " + equation3 + ") / 3";
			
			//make an initial guesses
			double a = 0.1 * rand.nextDouble(); 
			double b = 0.1 * rand.nextDouble();	
			double c = 0.1 * rand.nextDouble();	
			double d = min[i];
			double e = dy;
			double[] initialGuess = {a, b, c, d, e};			
			
			double[] y = new double[x.length];
			for (int j = 0 ; j < y.length ; j++) y[j] = yData[i][j];
			
			FittingResult nowFR = fitWithCurveFitter(equation, initialGuess, x, y);
						
			for (int j = 0 ; j < nowFR.params.length ; j++)	mFR.params[i][j] = nowFR.params[j];
			mFR.R2[i] = nowFR.R2;
			
		}
		
		return mFR;
		
	}

	/*
	 * public FittingResults fitFractionalFunction2Params(ObjectDetector.RadialModes
	 * myModes, int maxEval) {
	 * 
	 * FittingResults myFit = new FittingResults(); int numberOfStacksInTiff =
	 * myModes.maskingThreshold.length; double[][] mySpecialFunParams = new
	 * double[numberOfStacksInTiff][2]; double[] R2 = new
	 * double[numberOfStacksInTiff];
	 * 
	 * for (int i = 0 ; i < numberOfStacksInTiff ; i++) {
	 * 
	 * IJ.
	 * showStatus("Fitting correction function to soil air-phase illumination data ... "
	 * + (i + 1) + "/" + numberOfStacksInTiff);
	 * 
	 * double[] nowIll = new double[myModes.radius.length]; for (int j = 0 ; j <
	 * nowIll.length ; j++) nowIll[j] = myModes.maskedRadialMinima[i][j];
	 * 
	 * //init optimization LevenbergMarquardtOptimizer optimizer = new
	 * LevenbergMarquardtOptimizer(); CurveFitter<FractionalFunction2Params> fitter
	 * = new CurveFitter<FractionalFunction2Params>(optimizer);
	 * FractionalFunction2Params myFF = new FractionalFunction2Params();
	 * myFF.setup(nowIll[0], nowIll[nowIll.length - 1]);
	 * 
	 * //place the data for (int j = 0 ; j < nowIll.length ; j++)
	 * fitter.addObservedPoint(myModes.radius[j], nowIll[j]);
	 * 
	 * //fit it int cc = 0; double[] fittedLog = new double[2]; boolean hasConverged
	 * = false; while (hasConverged == false & cc < 10) { try {
	 * 
	 * Random rand = new Random();
	 * 
	 * //make an initial guess double a = 1 + rand.nextDouble(); double b = 1 +
	 * rand.nextDouble();
	 * 
	 * double[] initialGuess = {a, b};
	 * 
	 * fittedLog = fitter.fit(maxEval, myFF, initialGuess);
	 * 
	 * hasConverged = true;
	 * 
	 * } catch (ConvergenceException e) { } catch (NotStrictlyPositiveException e) {
	 * } catch (TooManyEvaluationsException e) { } finally { cc++; }
	 * 
	 * }
	 * 
	 * //create a vector for each radial coordinate //double[] corrFun = new
	 * double[(int)max(myModes.radius)]; for (int j = 0 ; j < 2 ; j++)
	 * mySpecialFunParams[i][j] = fittedLog[j];
	 * 
	 * //calculate R2 if (fittedLog[0] == 0) R2[i] = 0; else { double[]
	 * squaredResiduals = new double[nowIll.length]; double[] squaredTotals = new
	 * double[nowIll.length]; double meanOfData = StatUtils.mean(nowIll); for (int j
	 * = 0 ; j < nowIll.length ; j++) { double residual = nowIll[j] -
	 * myFF.value(myModes.radius[j], fittedLog); double total = nowIll[j] -
	 * meanOfData; squaredResiduals[j] = residual * residual; squaredTotals[j] =
	 * total * total; } R2[i] = 1 - (StatUtils.sum(squaredResiduals) /
	 * StatUtils.sum(squaredTotals));
	 * 
	 * }
	 * 
	 * }
	 * 
	 * myFit.numberOfParams = 2; myFit.params = mySpecialFunParams; myFit.R2 = R2;
	 * 
	 * return myFit;
	 * 
	 * }
	 */

//	public FittingResults fitGeneralizedLogisticFunction(ObjectDetector.RadialModes myModes, int maxEval) {
//				
//		FittingResults myFit = new FittingResults();
//		TailoredMaths m = new TailoredMaths();
//		AndAllTheRest aa = new AndAllTheRest(); 
//		
//		int numberOfStacksInTiff = myModes.maskingThreshold.length;
//		double[][] myLogisticParams = new double[numberOfStacksInTiff][6];		
//		double[] R2 = new double[numberOfStacksInTiff];
//	    
//		for (int i = 0 ; i < numberOfStacksInTiff ; i++) {
//			
//			IJ.showStatus("Fitting correction function to soil matrix illumination data ... " + (i + 1) + "/" + numberOfStacksInTiff);
//			
//			double[] nowIll = new double[myModes.radius.length]; 
//			for (int j = 0 ; j < nowIll.length ; j++) nowIll[j] = myModes.maskedRadialModes[i][j];
//						
//			LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();			
//		    CurveFitter<Logistic.Parametric> fitter = new CurveFitter<Logistic.Parametric>(optimizer);
//		    Parametric myLogL = new Logistic.Parametric();
//		 	    		    
//		    //check if A is on the left of K
//		    double A = m.min(nowIll);
//		    double K = m.max(nowIll);
//		    int ap, kp;
//		    ArrayList<Integer> aPos = aa.findFirstPositionInArray(nowIll, A);
//		    ArrayList<Integer> kPos = aa.findFirstPositionInArray(nowIll, K);		    		    
//		    //if (aPos.size() == 1) ap = aPos.get(0);
//		    //if (kPos.size() == 1) kp = kPos.get(0);
//		    ap = aPos.get(0);
//		    kp = kPos.get(0);
//		    
//		    //place the data 
//		    for (int j = 0 ; j < nowIll.length ; j++) fitter.addObservedPoint(myModes.radius[j], nowIll[j]);
//		   
//		    //fit it
//		    int cc = 0;			
//		    double[] fittedLog = new double[6]; 
//		    boolean hasConverged = false;
//		    while (hasConverged == false & cc < 10) {
//		    	try {
//	
//		    		Random rand = new Random();
//		    				    		
//		    		//make an initial guess
//				    double B = 0.1 * rand.nextDouble(); 
//				    double Q = 10 * rand.nextDouble();
//				    double ny  = rand.nextDouble();
//				    double M = m.max(myModes.radius) / 3 + rand.nextDouble() * m.max(myModes.radius) / 3;
//				    
//				    if (ap > kp) B = -B;				    	
//		    		
//				    double[] initialGuess = {K, M, B, Q, A, ny};
//				    
//		    		fittedLog = fitter.fit(maxEval, myLogL, initialGuess);				    		
//		    		
//		    		hasConverged = true;
//		    		
//		    	} catch (ConvergenceException e) {		    		
//		    	} catch (NotStrictlyPositiveException e) {		    		
//		    	} catch (TooManyEvaluationsException e) {		    		
//		    	} finally {
//		    		cc++;
//		    	}
//		    	
//		    }
//		    
//		    //create a vector for each radial coordinate
//		    //double[] corrFun = new double[(int)max(myModes.radius)];
//		    for (int j = 0 ; j < 6 ; j++) myLogisticParams[i][j] = fittedLog[j];
//		    
//		    //calculate R2
//		    if (fittedLog[0] == 0) R2[i] = 0;
//		    else {
//		    	double[] fittedCurve = new double[nowIll.length];
//		    	double[] squaredResiduals = new double[nowIll.length];	
//		    	double[] squaredTotals = new double[nowIll.length];	
//		    	double meanOfData = StatUtils.mean(nowIll);
//		    	for (int j = 0 ; j < nowIll.length ; j++) {
//		    		fittedCurve[j] = myLogL.value(myModes.radius[j], fittedLog);
//		    		double residual = nowIll[j] - myLogL.value(myModes.radius[j], fittedLog);
//		    		double total = nowIll[j] - meanOfData;
//		    		squaredResiduals[j] = residual * residual;
//		    		squaredTotals[j] = total * total;
//		    	}
//		    	R2[i] = 1 - (StatUtils.sum(squaredResiduals) / StatUtils.sum(squaredTotals));
//		    	
//		    	//check fitted Curve for reasonable values.. if not reasonable, set R2 to 0;
//		    	double[] dR = new double[myModes.radius.length - 1];
//		    	double[] diffFittedCurve = new double[myModes.radius.length - 1];
//		    	double[] gradient = new double[myModes.radius.length - 1];
//		    	for (int j = 0 ; j < dR.length - 1; j++) {
//		    		dR[j] = myModes.radius[j+1] - myModes.radius[j];
//		    		diffFittedCurve[j] = fittedCurve[j+1] - fittedCurve[j];
//		    		gradient[j] = diffFittedCurve[j] / dR[j];		    		
//		    	}
//		    	
//		    	double thresh = 0.25;
//		    	for (int j = 1 ; j < gradient.length - 1 ; j++) {
//		    		if (gradient[j-1] < -thresh & gradient[j] > thresh & gradient[j+1] < -thresh) R2[i] = -1;
//		    		if (gradient[j-1] > thresh & gradient[j] < -thresh & gradient[j+1] > thresh) R2[i] = -1;
//		    	}
//		    	
//		    	//have a look at it...
//		    	//DisplayThings dT = new DisplayThings();
//		    	//dT.plotXYXY(myModes.radius, nowIll, myModes.radius, fittedCurve, "GLF" + i, "radial coord.", "median grey value");
//		    	
//		    }
//		    
//		}		
//		
//		myFit.numberOfParams = 6;
//		myFit.params = myLogisticParams;
//		myFit.R2 = R2;
//		
//		return myFit;
//		
//	}

//	public FittingResults fitHyperbolicFunction2Params(ObjectDetector.RadialModes myModes, int maxEval) {
//		
//		FittingResults myFit = new FittingResults();
//		int numberOfStacksInTiff = myModes.maskingThreshold.length;
//		double[][] mySpecialFunParams = new double[numberOfStacksInTiff][2];
//		double[] R2 = new double[numberOfStacksInTiff];
//		
//		for (int i = 0 ; i < numberOfStacksInTiff ; i++) {
//			
//			IJ.showStatus("Fitting correction function to soil air-phase illumination data ... " + (i + 1) + "/" + numberOfStacksInTiff);
//			
//			double[] nowIll = new double[myModes.radius.length]; 
//			for (int j = 0 ; j < nowIll.length ; j++) nowIll[j] = myModes.maskedRadialMinima[i][j]; 
//			
//			//init optimization
//			LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();			
//		    CurveFitter<HyperbolicFunction2Params> fitter = new CurveFitter<HyperbolicFunction2Params>(optimizer);
//		    HyperbolicFunction2Params myLogF = new HyperbolicFunction2Params();
//		    myLogF.setup(StatUtils.max(myModes.radius), nowIll[nowIll.length - 1]);
//		    
//		    //place the data 
//		    for (int j = 0 ; j < nowIll.length ; j++) fitter.addObservedPoint(myModes.radius[j], nowIll[j]);
//		   
//		    //fit it
//		    int cc = 0;			
//		    double[] fittedLog = new double[2]; 
//		    boolean hasConverged = false;
//		    while (hasConverged == false & cc < 10) {
//		    	try {
//	
//		    		Random rand = new Random();
//		    				    		
//		    		//make an initial guess
//				    double a = nowIll[0];
//				    double b = 0.05 * rand.nextDouble();
//		    		
//				    double[] initialGuess = {a, b};
//				    
//		    		fittedLog = fitter.fit(maxEval, myLogF, initialGuess);			    		
//		    		
//		    		hasConverged = true;
//		    		
//		    	} catch (ConvergenceException e) {		    		
//		    	} catch (NotStrictlyPositiveException e) {		    		
//		    	} catch (TooManyEvaluationsException e) {		    		
//		    	} finally {
//		    		cc++;
//		    	}
//		    	
//		    }
//		    
//		    //create a vector for each radial coordinate
//		    //double[] corrFun = new double[(int)max(myModes.radius)];
//		    for (int j = 0 ; j < 2 ; j++) mySpecialFunParams[i][j] = fittedLog[j];
//		    
//		    //calculate R2
//		    if (fittedLog[0] == 0) R2[i] = 0;
//		    else {
//		    	double[] fittedCurve = new double[nowIll.length];
//		    	double[] squaredResiduals = new double[nowIll.length];	
//		    	double[] squaredTotals = new double[nowIll.length];	
//		    	double meanOfData = StatUtils.mean(nowIll);
//		    	for (int j = 0 ; j < nowIll.length ; j++) {
//		    		fittedCurve[j] = myLogF.value(myModes.radius[j], fittedLog);
//		    		double residual = nowIll[j] - myLogF.value(myModes.radius[j], fittedLog);
//		    		double total = nowIll[j] - meanOfData;
//		    		squaredResiduals[j] = residual * residual;
//		    		squaredTotals[j] = total * total;
//		    	}
//		    	R2[i] = 1 - (StatUtils.sum(squaredResiduals) / StatUtils.sum(squaredTotals));
//		    	
//		    	//check fitted Curve for reasonable values.. if not reasonable, set R2 to 0;
//		    	double[] dR = new double[myModes.radius.length - 1];
//		    	double[] diffFittedCurve = new double[myModes.radius.length - 1];
//		    	double[] gradient = new double[myModes.radius.length - 1];
//		    	for (int j = 0 ; j < dR.length - 1; j++) {
//		    		dR[j] = myModes.radius[j+1] - myModes.radius[j];
//		    		diffFittedCurve[j] = fittedCurve[j+1] - fittedCurve[j];
//		    		gradient[j] = diffFittedCurve[j] / dR[j];		    		
//		    	}
//		    	
//		    	double thresh = 0.25;
//		    	for (int j = 1 ; j < gradient.length - 1 ; j++) {
//		    		if (gradient[j-1] < -thresh & gradient[j] > thresh & gradient[j+1] < -thresh) R2[i] = -1;
//		    		if (gradient[j-1] > thresh & gradient[j] < -thresh & gradient[j+1] > thresh) R2[i] = -1;
//		    	}
//		    	
//		    	//have a look at it...
//		    	//DisplayThings dT = new DisplayThings();
//		    	//dT.plotXYXY(myModes.radius, nowIll, myModes.radius, fittedCurve, "HF" + i, "radial coord.", "median grey value");
//		    	
//		    }
//		    
//		}		
//		
//		myFit.numberOfParams = 2;
//		myFit.params = mySpecialFunParams;
//		myFit.R2 = R2;
//		
//		return myFit;
//		
//	}
	
	public class FittedEllipse {
		
		public double xCenter;
		public double yCenter;
		public double zCenter;
		public double minorRadius;
		public double majorRadius;				
		public double theta;		
		public double R2;
		
		public ArrayList<Double> relDevFromMedianRadiusOfDiscPoints;
		public double medianAbsDevFromMedianRadius;
		public double medianRadius;
		
		public FittedEllipse duplicate() {
			
			FittedEllipse fD = new FittedEllipse();
			
			fD.xCenter = this.xCenter;
			fD.yCenter = this.yCenter;
			fD.zCenter = this.zCenter;
			fD.minorRadius = this.minorRadius;
			fD.majorRadius = this.majorRadius;				
			fD.theta = this.theta;		
			fD.R2 = this.R2;
			
			fD.relDevFromMedianRadiusOfDiscPoints = this.relDevFromMedianRadiusOfDiscPoints;
			fD.medianAbsDevFromMedianRadius = this.medianAbsDevFromMedianRadius;
			fD.medianRadius = this.medianRadius;
			
			return fD;
		}
		
	}
	
	public FittedEllipse doRobustEllipseFit(int i, float[] xD, float[] yD, double[] myAngle, ImageProcessor myIP, MenuWaiter.ColumnFinderMenuReturn jCFS, double para) {
		
		FittedEllipse fE = new FittedEllipse();
		EllipseFitter jEF = new EllipseFitter();
		
		//kick out zero entries
		ArrayList<Integer> badones = new ArrayList<Integer>();
		for (int j = 0 ; j < xD.length ; j++) {
			if (xD[j] == 0 | yD[j] == 0) badones.add(j);
		}
		float[] xC = new float[xD.length - badones.size()];
		float[] yC = new float[xD.length - badones.size()];
		double[] angleC = new double[xD.length - badones.size()];
		int precc = 0;
		for (int j = 0 ; j < xD.length ; j++) {
			if (badones.contains(j)); 
			else {
				xC[precc] = xD[j];
				yC[precc] = yD[j];
				angleC[precc] = myAngle[j];
				precc++;
			}			
		}
		
		//fit with all available values
		PolygonRoi pRoi = new PolygonRoi(xC, yC, Roi.POLYLINE); // create pRoi with outer Wall coordinates
		ImageProcessor copy1 = myIP.duplicate();
		copy1.setRoi(pRoi); 
		jEF.fit(copy1, copy1.getStatistics());
		FitStuff.GoodnessOfFit gOF = calculateR2OfEllipsoidFit(angleC, xC, yC, jEF);		
		double[] relDevFromMedianRadius = new double[xC.length];
		for (int j = 0 ; j < xC.length ; j++) {
			relDevFromMedianRadius[j] = (gOF.d2m[j] - gOF.median2d2m) / gOF.median2d2m;
		}
		
		PolygonRoi sRoi = new PolygonRoi(xC, yC, Roi.POLYLINE); // create pRoi with outer Wall coordinates
		ImageProcessor copy2 = myIP.duplicate();
		copy2.setRoi(sRoi);				
		
		ImagePlus newImg = new ImagePlus("layer " + i, copy2);
		Overlay myO = new Overlay(sRoi);	
		
		if (jCFS.showFit) {
			PointRoi innerFoundEdges = new PointRoi(xC, yC, yC.length);		
			myO.add(innerFoundEdges);  
		
			copy2.setColor(Color.YELLOW);
			Font font = new Font("Verdana", Font.PLAIN, 40);
			TextRoi tRoi = new TextRoi((int)(0.01f * newImg.getWidth()), (int)(0.01f * newImg.getHeight()), String.format("%1.2f\t",(float)para), font);
			tRoi.drawPixels(copy2);
		
			ContrastEnhancer myCE = new ContrastEnhancer();
			myCE.stretchHistogram(copy2, 0.5);
			myO.setStrokeColor(Color.RED);
			newImg.setOverlay(myO);
			newImg.updateAndDraw();
			newImg.show();
		
			IJ.wait(500);
			newImg.hide();
			newImg.flush();		
		}
		
		//assign output variables
		fE.xCenter = jEF.xCenter;
		fE.yCenter = jEF.yCenter;
		fE.zCenter = i;
		fE.minorRadius = jEF.minor / 2;
		fE.majorRadius = jEF.major / 2;				
		fE.theta = jEF.theta;		
		fE.R2 = gOF.R2;		
		
		//try to do the fit again but without outliers
		ArrayList<Double> kickout = new ArrayList<Double>();
		if (gOF.R2 < jCFS.r2Thresh) {			
				
			//kick out values that are above threshold
			double kickThresh = 0.01;
			int maxDev = 0;
			double maxDEV = 0;
			for (int j = 0 ; j < relDevFromMedianRadius.length ; j++) {				
				if (Math.abs(relDevFromMedianRadius[j]) > kickThresh) kickout.add(relDevFromMedianRadius[j]);	
				if (Math.abs(relDevFromMedianRadius[j]) > maxDEV) {
					maxDev = j;
					maxDEV = Math.abs(relDevFromMedianRadius[j]);
				}
			}	
			
			//also add the badones as negative ones.. 

			for (int j = 0 ; j < badones.size() ; j++) {				
				kickout.add((double) -1);				
			}
				 
			//remove kickouts stepwise					
			int cc = 0;
			ArrayList<Integer> kicked = new ArrayList<Integer>();
			kicked.add(maxDev);
			
			//abort fitting if too many points were not found
			if (badones.size() > xC.length / 6) cc = 999;
			
			while (gOF.R2 < jCFS.r2Thresh & cc <= jCFS.maxNumberOfOutliers4PerimeterFits) {
			
				//compile new array without kickouts	
				float[] xOS = new float[xC.length - kicked.size()];
				float[] yOS = new float[xC.length - kicked.size()];
				double[] angleS = new double[xC.length - kicked.size()];
				
				int nowcc = 0;
				for (int j = 0 ; j < xC.length ; j++){
					if (!kicked.contains(j)) {
						xOS[nowcc] = xC[j];
						yOS[nowcc] = yC[j];
						angleS[nowcc] = angleC[j];	
						nowcc++;
					}						
				}
				
				sRoi = new PolygonRoi(xOS, yOS, Roi.POLYLINE); // create pRoi with outer Wall coordinates
				copy2 = myIP.duplicate();
				copy2.setRoi(sRoi);				
				
				if (jCFS.showFit) {
					newImg = new ImagePlus("layer " + i, copy2);
					copy2.setColor(Color.YELLOW);
					myO = new Overlay(sRoi);	
					PointRoi innerFoundEdges = new PointRoi(xOS, yOS, yOS.length);
					myO.add(innerFoundEdges);  
					
					Font font = new Font("Verdana", Font.PLAIN, 40);
					TextRoi tRoi = new TextRoi((int)(0.01f * newImg.getWidth()), (int)(0.01f * newImg.getHeight()), String.format("%1.2f\t",(float)para), font);
					tRoi.drawPixels(copy2);
					
					ContrastEnhancer myCE = new ContrastEnhancer();
					myO.setStrokeColor(Color.RED);
					myCE.stretchHistogram(copy2, 0.5);
					newImg.setOverlay(myO);
					newImg.updateAndDraw();
					newImg.show();
					
					IJ.wait(500);
					newImg.hide();
					newImg.flush();
				}
								

				EllipseFitter jEFS = new EllipseFitter();
				jEFS.fit(copy2, copy2.getStatistics());			
				gOF = calculateR2OfEllipsoidFit(angleS, xOS, yOS, jEFS);
				double[] newRelDevFromMedianRadius = new double[xC.length];
				ArrayList<Double> newAbsDevFromMedianRadius = new ArrayList<Double>();
				int kickRedu = 0;
				for (int j = 0 ; j < xOS.length ; j++) {
					if (!kicked.contains(j)) {
						newRelDevFromMedianRadius[j] = Math.abs(gOF.d2m[j - kickRedu] - gOF.median2d2m) / gOF.median2d2m;		
					}
					else {
						newRelDevFromMedianRadius[j] = 0;
						newAbsDevFromMedianRadius.add(gOF.d2m[j - kickRedu] - gOF.median2d2m);
						kickRedu++;
					}
				}
				maxDev = 0;
				maxDEV = 0;
				for (int j = 0 ; j < newRelDevFromMedianRadius.length ; j++) {				
					if (newRelDevFromMedianRadius[j] > maxDEV) {
						maxDev = j;
						maxDEV = newRelDevFromMedianRadius[j];
					}
				}	
				
				kicked.add(maxDev);
				
				double[] absDeviation = new double[newAbsDevFromMedianRadius.size()];
				for (int j = 0 ; j < absDeviation.length ; j++) absDeviation[j] = newAbsDevFromMedianRadius.get(j);
				
				fE.xCenter = jEFS.xCenter;
				fE.yCenter = jEFS.yCenter;
				fE.zCenter = i;
				fE.minorRadius = jEFS.minor / 2;
				fE.majorRadius = jEFS.major / 2;				
				fE.theta = jEFS.theta;		
				fE.R2 = gOF.R2;
				fE.relDevFromMedianRadiusOfDiscPoints = kickout;
				
				fE.medianAbsDevFromMedianRadius = StatUtils.percentile(absDeviation, 50);
				fE.medianRadius =  gOF.median2d2m;
				
				cc++;
			}						
		}
		
		return fE;
		
	}

	public class FractionalFunction2Params implements ParametricUnivariateFunction {
	
		private double top;
		private double bottom;

		public void setup(double top, double bottom) {			
			this.top = top;
			this.top = bottom;
		}
	
		public double value(double x, double[] params) {
		
			double a = params[0];				
			double b = params[1];			
		
			double f = top - (top - bottom) * a * Math.pow(x, b);
		
			return f;			
		}
	
	
		public double[] gradient(double x, double[] params) {
		
			double a = params[0];			
			double b = params[1];			
		
			double[] gradient = new double[2];
		
			gradient[0] = - (top - bottom) * Math.pow(x, b);		    
			gradient[1] = - (top - bottom) * a * b * Math.pow(x, b - 1);	
		
			return gradient;		
		}
	
	}

	public class HyperbolicFunction2Params implements ParametricUnivariateFunction {
	
		private double rmax;
		private double bottom;

		public void setup(double rmax, double bottom) {
			this.rmax = rmax;
			this.bottom = bottom;
		}
	
	public double value(double x, double[] params) {
		
			double a = params[0];				
			double b = params[1];			
		
			double f = a - (a - bottom) / (1 + b * (rmax - x));
		
			return f;			
		}		
	
	public double[] gradient(double x, double[] params) {
		
			double a = params[0];			
			double b = params[1];			
		
			double[] gradient = new double[2];
		
			gradient[0] = 1 - 1 / (1 + b * (rmax - x));		    
			gradient[1] = a * (rmax - x) / Math.pow((1 + b * (rmax - x)), 2) - bottom * (rmax - x) / Math.pow((1 + b * (rmax - x)), 2);	
		
			return gradient;			
		}
	
	}

	public class GoodnessOfFit {
	
		double[] srx;	
		double[] sry;
	
		double[] d2m; // distance to midpoint
		double median2d2m;
	
		double R2;
	
	}

	public class LinearFunction {

		double slope;
		double intercept;
		double R2;
	
	}

	public class FittingResult {
		
		double[] params;
		double R2;
	
	}
	
	public class FittingResults {
	
		int numberOfParams;
		double[][] params;
		double[] R2;
	
	}
}