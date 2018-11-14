package test;

import java.util.ArrayList;

import base.Array;
import base.Fit;
import base.LBFGS;
import base.Minimizer;
import base.Model;
import base.sLBFGS;

public class GenericTestFile {
	
	public static void main(String[] args) {
		SVMTest svm = new SVMTest(0.001);
		Minimizer min = new Minimizer();
		Fit fit = null;
		
		// Minimize with LBFGS
		double tStart = System.nanoTime();
		min = new LBFGS(svm, 200, 1E-5, 1000, true, true);
		double tEnd = System.nanoTime();
		
		try {
			fit = min.doMinimize(null, null);
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println(fit.likelihood);
	
		min = new sLBFGS(svm, 10, 100, 10, 30, 10, 1000, 0.01, 0, true);
		
		try {
			fit = min.doMinimize(null, null);
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("LBFGS Time consumed: "+(tEnd-tStart)/1E9);
		System.out.println(fit.likelihood);
	}
}
