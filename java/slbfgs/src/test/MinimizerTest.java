package test;

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Random;

import base.Array;
import base.Fit;
import base.LBFGS;
import base.Minimizer;
import base.sLBFGS;

public class MinimizerTest {
	public static int testLoops = 10;				//Number of random starts to test
	public static double xlim = 1, ylim = 1;		//Bounds of random starts
	public static int maxMemoryDepth = 10, maxIterations = 100;
	public static double convergence = 1E-5;
	public static double stochStepSize = 0.01;
	public static int dimensionality = 2;
	
	public static void main(String[] args) {
		double tStart, tEnd;
		//Initialize random generator
		double[] seed, fitTime, iterations, dataLoops, successes, correctPos = null;
		Random generator = new Random();
		
		//Create null output stream
		PrintStream originalOut	= System.out;
		PrintStream originalErr = System.err;
		PrintStream nullStream  = new PrintStream(new NullOutputStream());  
		
		//Loop over function types
		for (int fType = 1; fType<=0; fType++) {
			//Initialize test functions
			MinimizerTestFunctions testFunc = new MinimizerTestFunctions(fType, dimensionality);
			Fit fit = null;
			Minimizer minimize = new Minimizer();
			
			//Define final position
			correctPos = new double[dimensionality];
			switch (fType) {
				case 1:		for (int i=0; i<dimensionality; i++) {
								correctPos[i] = 1;
							}
							break;
				case 2:		for (int i=0; i<dimensionality; i++) {
								correctPos[i] = 0;
							}
							break;
				case 3:		for (int i=0; i<dimensionality; i+=2) {
								correctPos[i] = 3;
								correctPos[i+1] = 0.5;
							}
							break;
				case 4:		for (int i=0; i<dimensionality; i++) {
								correctPos[i] = 0;
							}
							break;
				case 5:		for (int i=0; i<dimensionality; i+=2) {
								correctPos[i] = 0;
								correctPos[i+1] = -1;
							}
			}
			
			//Loop over multiple random starts
			System.setOut(nullStream);
			System.setErr(nullStream);
			fitTime = new double[2];
			iterations = new double[2];
			successes = new double[2];
			dataLoops = new double[2];
			for (int currLoop = 0; currLoop<testLoops; currLoop++) {
				//create new seed
				seed = new double[dimensionality];
				for (int i=0; i<dimensionality; i++) {
					seed[i] = generator.nextDouble()*xlim*2-xlim;
				}
				
				// Run LBFGS with MCSearch
				minimize = new LBFGS(testFunc, maxMemoryDepth, convergence, maxIterations, true,false);
				try {
					tStart = System.nanoTime();
					fit = minimize.doMinimize(seed, null);
					tEnd = System.nanoTime();
					if (Array.dist(fit.finalPosition, correctPos)<2*convergence) {
						successes[0]++;
						fitTime[0]		+= (tEnd-tStart)/1E9;
						iterations[0]	+= fit.fitSteps;
						dataLoops[0]	+= fit.dataLoops;
					}
				} catch (Exception e) {
					//Do nothing
				}
				
				// Run sLBFGS
				minimize = new sLBFGS(testFunc, 1, 1, maxMemoryDepth, maxIterations, 10, 500, stochStepSize, convergence, 0, true);
				try {
					tStart = System.nanoTime();
					fit = minimize.doMinimize(seed, null);
					tEnd = System.nanoTime();
					if (Array.dist(fit.finalPosition, correctPos)<2*convergence) {
						successes[1]++;
						fitTime[1]		+= (tEnd-tStart)/1E9;
						iterations[1]	+= fit.fitSteps;
						dataLoops[1]	+= fit.dataLoops;
					}
				} catch (Exception e) {
					//Do nothing
				}
			}
			System.setOut(originalOut);
			System.setErr(originalErr);
			System.out.println("Function Type:\t"+testFunc.fName);
			System.out.println("\t\tLBFGS\t\tsLBFGS");
			System.out.println("Successes:     \t"+successes[0]+"\t\t"+successes[1]);
			System.out.println("Fit Time:      \t"+fitTime[0]+"\t"+fitTime[1]);
			System.out.println("Iterations:    \t"+iterations[0]+"\t\t"+iterations[1]);
			System.out.println("Data Loops:    \t"+dataLoops[0]+"\t\t"+dataLoops[1]);
			System.out.println("\n");
		}
		
		// Run SVM test. Begin with a convergent L-BFGS run
		SVM svm = new SVM(0.001);
		Minimizer min = new Minimizer();
		Fit lbfgsFit = null, slbfgsFit = null;
		fitTime = new double[2];
		iterations = new double[2];
		successes = new double[2];
		dataLoops = new double[2];
		
		// Minimize with LBFGS
		System.setOut(nullStream);
		System.setErr(nullStream);
		for (int currLoop = 0; currLoop<testLoops; currLoop++) {
			// Here, memory depth = 200
			try {
				min = new LBFGS(svm, 200, convergence, maxIterations, false, false);
				tStart = System.nanoTime();
				lbfgsFit = min.doMinimize(null, null);
				tEnd = System.nanoTime();
				fitTime[0] += (tEnd-tStart)/1E9;
				iterations[0] += lbfgsFit.fitSteps;
				successes[0]++;
				dataLoops[0] += lbfgsFit.dataLoops;
			} catch (Exception e) {
				// Do nothing
			}
			
			try {
				min = new sLBFGS(svm, 20, 200, maxMemoryDepth, maxIterations, 10, 500, stochStepSize, convergence, 0, false);
				tStart = System.nanoTime();
				slbfgsFit = min.doMinimize(null, null);
				tEnd = System.nanoTime();
				fitTime[1] += (tEnd-tStart)/1E9;
				iterations[1] += slbfgsFit.fitSteps;
				successes[1]++;
				dataLoops[1] += slbfgsFit.dataLoops;
			} catch (Exception e) {
				// Do nothing
			}
		}
		System.setOut(originalOut);
		System.setErr(originalErr);
		System.out.println("Function Type:\t"+svm.fName);
		System.out.println("\t\tLBFGS\t\tsLBFGS");
		System.out.println("Successes:     \t"+successes[0]+"\t\t"+successes[1]);
		System.out.println("Fit Time:      \t"+fitTime[0]+"\t"+fitTime[1]);
		System.out.println("Iterations:    \t"+iterations[0]+"\t\t"+iterations[1]);
		System.out.println("Data Loops:    \t"+dataLoops[0]+"\t\t"+dataLoops[1]);
		System.out.println("\n");
	}
	
	private static class NullOutputStream extends OutputStream {
	    @Override
	    public void write(int b){
	         return;
	    }
	    @Override
	    public void write(byte[] b){
	         return;
	    }
	    @Override
	    public void write(byte[] b, int off, int len){
	         return;
	    }
	    public NullOutputStream(){
	    }
	}
}
