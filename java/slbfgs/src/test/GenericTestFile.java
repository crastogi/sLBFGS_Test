package test;

import base.Fit;
import base.LBFGS;
import base.Minimizer;
import base.sLBFGS;

public class GenericTestFile {
	
	public static void main(String[] args) {
		double tStart = 0, tEnd = 0;
		Fit fit = null;
		SVMTest svm = new SVMTest(0.001);
		Minimizer min = new Minimizer();
		
		// Minimize with LBFGS
		min = new LBFGS(svm, 200, 1E-5, 1000, true, true);
		try {
			tStart = System.nanoTime();
			fit = min.doMinimize(null, null);
			tEnd = System.nanoTime();
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println(fit.likelihood);
		System.out.println("LBFGS Time consumed: "+(tEnd-tStart)/1E9);
	
		// Now minimize with sLBFGS
		min = new sLBFGS(svm, 20, 200, 10, 30, 10, 500, 0.1, 1E-5, 0, true);
		try {
			tStart = System.nanoTime();
			fit = min.doMinimize(new double[svm.nDim], null);
			tEnd = System.nanoTime();
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println(fit.likelihood);
		System.out.println("sLBFGS Time consumed: "+(tEnd-tStart)/1E9);
	}
}
