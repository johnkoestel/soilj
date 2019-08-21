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
import ij.plugin.PlugIn;

import SoilJ.tools.DisplayThings;
import SoilJ.tools.MorphologyAnalyzer;
import SoilJ.tools.MorphologyAnalyzer.ProfileStatistics;
import SoilJ.tools.RollerCaster;

/** 
 * PlotVerticalProfile is a SoilJ plugin that extracts vertical profiles of the greyscale
 * statistics within horizontal slices of 3-D images. 
 * 
 * @author John Koestel
 *
 */

public class CalcT1ContributionFrom2ComplexMRIImages_ extends ImagePlus implements PlugIn  {

	public void run(String arg) {

		String pathSep = "/";
		
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = JointThresholdDetection_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		
		
		
	}
	
}