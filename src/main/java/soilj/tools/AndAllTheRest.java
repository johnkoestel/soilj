package SoilJ.tools;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class AndAllTheRest implements PlugIn {

	
	public void run(String arg) {
				//ok, this is not needed..
	}

	public ArrayList<Integer> findFirstPositionInArray(double[] myArray, double value) {
		
		ArrayList<Integer> myPosition = new ArrayList<Integer>();
		for (int i = 0 ; i < myArray.length ; i++) if (myArray[i] == value) myPosition.add(i);
		
		return myPosition;
				
	}
	
	public boolean isContainedIn(int j, int[] array) {
	
		boolean isItTrue = false;
		
		for (int i = 0 ; i < array.length ; i++) {
			if (j == array[i]) {
				isItTrue = true;
			}
		}
		
		return isItTrue;
	}

	public static int[] getIndicesInOrder(double[] array) { //stolen from a JAVA student in the WWW.. god bless him!!
	    Map<Integer, Double> map = new HashMap<Integer, Double>(array.length);
	    for (int i = 0; i < array.length; i++)
	        map.put(i, array[i]);
	
	    List<Entry<Integer, Double>> l = new ArrayList<Entry<Integer, Double>>(map.entrySet());
	
	    Collections.sort(l, new Comparator<Entry<?, Double>>() {
	            @Override
	            public int compare(Entry<?, Double> e1, Entry<?, Double> e2) {
	                return e2.getValue().compareTo(e1.getValue());
	            }
	        });
	
	    int[] result = new int[array.length];
	    for (int i = 0; i < result.length; i++)
	        result[i] = l.get(i).getKey();
	
	    return result;
	}
	
}