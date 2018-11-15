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
	public static int testLoops = 1000;			//Number of random starts to test
	public static double xlim = 1, ylim = 1;		//Bounds of random starts
	public static int maxMemoryDepth = 10, maxIterations = 30;
	public static double convergence = 1E-5;
	public static double d0 = 1;
	public static double stochStepSize = 0.01;
	public static int dimensionality = 2;
	
	public static void main(String[] args) {
		//Initialize random generator
		double[] seed, fitTime, functionCalls, iterations, successes, correctPos = null;
		Random generator = new Random();
		
		//Create null output stream
		PrintStream original	= System.out;
		PrintStream nullStream  = new PrintStream(new NullOutputStream());  
		
		//Loop over function types
		for (int fType = 1; fType<=5; fType++) {
			//Initialize test functions (2D case only)
			MinimizerTestFunctions testFunc = new MinimizerTestFunctions(fType, dimensionality);
			Fit fit = null;
			Minimizer minimize = new Minimizer();
			
			//Define final position
			switch (fType) {
				case 1:		correctPos = new double[]{1,1};
							break;
				case 2:		correctPos = new double[]{0,0};
							break;
				case 3:		correctPos = new double[]{3,0.5};
							break;
				case 4:		correctPos = new double[]{0,0};
							break;
				case 5:		correctPos = new double[]{0,-1};
			}
			
			//Loop over multiple random starts
			System.setOut(nullStream);
			System.setErr(nullStream);
			fitTime = new double[5];
			functionCalls = new double[5];
			iterations = new double[5];
			successes = new double[5];
			for (int currLoop = 0; currLoop<testLoops; currLoop++) {
				//create new seed
				seed = new double[dimensionality];
				for (int i=0; i<dimensionality; i++) {
					seed[i] = generator.nextDouble()*xlim*2-xlim;
				}
				
				minimize = new LBFGS(testFunc, maxMemoryDepth, convergence, maxIterations, false, false);
				try {
					fit = minimize.doMinimize(seed, null);
					if (Array.dist(fit.finalPosition, correctPos)<2*convergence) {
						successes[2]++;
						fitTime[2] += fit.fitTime;
						functionCalls[2] += fit.functionCalls;
						iterations[2] += fit.fitSteps;
					}
				} catch (Exception e) {
					//Do nothing
				}
				
				minimize = new LBFGS(testFunc, maxMemoryDepth, convergence, maxIterations, true,false);
				try {
					fit = minimize.doMinimize(seed, null);
					if (Array.dist(fit.finalPosition, correctPos)<2*convergence) {
						successes[3]++;
						fitTime[3] += fit.fitTime;
						functionCalls[3] += fit.functionCalls;
						iterations[3] += fit.fitSteps;
					}
				} catch (Exception e) {
					//Do nothing
				}
				
				minimize = new sLBFGS(testFunc, 1, 1, maxMemoryDepth, maxIterations, 10, 500, stochStepSize, convergence, 0, true);
				try {
					fit = minimize.doMinimize(seed, null);
					if (Array.dist(fit.finalPosition, correctPos)<2*convergence) {
						successes[4]++;
						fitTime[4] += fit.fitTime;
						functionCalls[4] += fit.functionCalls;
						iterations[4] += fit.fitSteps;
					}
				} catch (Exception e) {
					//Do nothing
				}
			}
			System.setOut(original);
			System.setOut(original);
			System.out.println("Function Type:          "+fType);
			System.out.println("PS\tPS+Random Axis\tLBFGS\tLBFGS+MCSearch\tsLBFGS");
			System.out.print("Successes:     \t");
			Array.print(successes);
			System.out.print("Fit Time:      \t");
			Array.print(fitTime);
			System.out.print("Function Calls:\t");
			Array.print(functionCalls);
			System.out.print("Iterations:    \t");
			Array.print(iterations);
			System.out.println("\n");
		}
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
