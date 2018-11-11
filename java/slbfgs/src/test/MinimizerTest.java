package test;

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Random;

import base.Array;
import minimizers.*;

public class MinimizerTest {
	public static int testLoops = 10000;			//Number of random starts to test
	public static double xlim = 1, ylim = 1;		//Bounds of random starts
	public static int maxMemoryDepth = 7, maxIterations = 1000;
	public static double convergence = 1E-5;
	public static double d0 = 1;
	
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
			MinimizerTestFunctions testFunc = new MinimizerTestFunctions(fType, 2);
			MinimizerFit fit = null;
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
			fitTime = new double[4];
			functionCalls = new double[4];
			iterations = new double[4];
			successes = new double[4];
			for (int currLoop = 0; currLoop<testLoops; currLoop++) {
				//create new seed
				seed = new double[2];
				seed[0] = generator.nextDouble()*xlim*2-xlim;
				seed[1] = generator.nextDouble()*ylim*2-ylim;
				
				minimize = new LBFGS(testFunc, maxMemoryDepth, convergence, maxIterations, false, false);
				try {
					fit = (MinimizerFit) minimize.minimize(seed, null)[0];
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
					fit = (MinimizerFit) minimize.minimize(seed, null)[0];
					if (Array.dist(fit.finalPosition, correctPos)<2*convergence) {
						successes[3]++;
						fitTime[3] += fit.fitTime;
						functionCalls[3] += fit.functionCalls;
						iterations[3] += fit.fitSteps;
					}
				} catch (Exception e) {
					//Do nothing
				}
			}
			System.setOut(original);
			System.setOut(original);
			System.out.println("Function Type:          "+fType);
			System.out.println("PS\tPS+Random Axis\tLBFGS\tLBFGS+MCSearch");
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
