package test;

import base.*;

public class SVM_Sweep_Test {
	
	public static void main(String[] args) {
		int nSamples = 100;
		int memoryDepth = 50;
		boolean useReducedSpace = false;
		SVM svm = new SVM(0.001);
		int gradientBatch = 20;
		int hessianPeriod = 5;
		int stochIters = 10;
		int hessianBatch = 10*gradientBatch;
		double stepSize = .2;
		double epsilon = 1E-5;
		int maxEpochs = 1000;
		int svrgSubBatch = 1;
		double[] x_min, dists, dLoops, epochs;
		Fit fit = null;
		Minimizer min = new Minimizer();
		
		// Minimize with LBFGS to find true minimum
		min = new LBFGS(svm, 200, 1E-6, 1000, true, true);
		try {
			fit = min.doMinimize(new double[svm.nDim], null);
		}	catch (Exception e) {
			e.printStackTrace();
		}
		x_min = fit.finalPosition;
		
		// Test
		svm.evaluatedDataPoints = 0;
		min = new sLBFGS_kSVRG(svm, 4, 10, 250, 1000, 125, 5, .025, 1E-5, 0, true);
		//min = new kSVRG(svm, 4, true, false, 10, 1000, .5, 1E-5, false, true);
		//min = new SGD(svm, false, 20, 1000, true, 0, .01, 1E-5, false, true, true);
		min.setXStar(x_min);
		try {
			fit = min.doMinimize(new double[svm.nDim], null);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		if (true)
			return;
	
		// Now minimize with sLBFGS and compute averages
		epochs = new double[nSamples];
		dLoops = new double[nSamples];
		dists = new double[nSamples];
		for (int currSample = 0; currSample<nSamples; currSample++) {
			min = new sLBFGS_Test(svm, gradientBatch, hessianBatch, memoryDepth, maxEpochs, 
					hessianPeriod, stochIters, stepSize, epsilon, 0, false);
			if (useReducedSpace) {
				((sLBFGS_Test) min).useReducedSpace = true;
			}
			if (svrgSubBatch > 1) {
				((sLBFGS_Test) min).svrgSubBatch = svrgSubBatch;
			}
			try {
				fit = min.doMinimize(new double[svm.nDim], null);
			} catch (Exception e) {
				e.printStackTrace();
				continue;
			}
			dists[currSample] = Array.dist(fit.finalPosition, x_min);
			dLoops[currSample]= fit.dataLoops;
			epochs[currSample]= fit.fitSteps;
		}
		System.out.println("Distance to optimum: "+Array.mean(dists)+" ± "+Math.sqrt(Array.var(dists)));
		System.out.println("Data Loops: "+Array.mean(dLoops)+" ± "+Math.sqrt(Array.var(dLoops)));
		System.out.println("Epochs: "+Array.mean(epochs)+" ± "+Math.sqrt(Array.var(epochs)));
		Array.print(dists);
		Array.print(dLoops);
		Array.print(epochs);
	}
}
