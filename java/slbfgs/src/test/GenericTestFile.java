package test;

import base.Fit;
import base.LBFGS;
import base.Minimizer;

public class GenericTestFile {
	
	public static void main(String[] args) {
		SVMTest svm = new SVMTest(0.001);
		Minimizer min = new Minimizer();
		Fit fit = null;
		
		// Minimize with LBFGS
		min = new LBFGS(svm, 200, 1E-5, 1000, true, true);
		
		try {
			fit = min.doMinimize(null, null);
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println(fit.likelihood);
	}
}
