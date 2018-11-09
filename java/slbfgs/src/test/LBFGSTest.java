package test;

import java.util.Random;

import base.Array;

public class LBFGSTest {
	public static void main(String[] args) {
		int nLoops	= 100;
		double oldCalls= 0;
		double oldGNorm= 0;
		double oldFVal = 0;
		int oldLoops= 0;
		double newCalls= 0;
		double newIters= 0;
		double newGNorm= 0;
		double newFVal = 0;
		int newLoops= 0;
		double simCalls= 0;
		double simIters= 0;
		double simGNorm= 0;
		double simFVal = 0;
		int simLoops= 0;
		int nDim	= 2;
		int fType	= 1;
		int memDepth= 5;
		double[] x0	= new double[nDim];
		for (int j=0; j<nDim; j+=2) {				//Initialize the starting vector X.
			x0 [j]	= -1.2;							//Creates a vector with alternating values between -1.2 and 1.
			x0 [j+1]= 1;
		}
		
		LBFGSFunctions f = new LBFGSFunctions(fType, nDim);	
		LBFGSMinimizers testMinimizer	= new LBFGSMinimizers(f, memDepth, 1e-10, .0001, .9, 1e-20, 1e20, 1e-5, 400, 20, 1);
		LBFGSMinimizers testMinimizerMC	= new LBFGSMinimizers(f, memDepth, 1e-10, .0001, .9, 1e-20, 1e20, 1e-5, 400, 20, 2);
		Random generator = new Random();

		for (int loop=0; loop<nLoops; loop++) {
			for (int i=0; i<nDim; i++) {
				x0[i] = generator.nextDouble()*.5-0;
			}
			
			boolean provideHk0	= false;
			int iterations		= 0;
			double functionValue= 0;
			double epsilon		= 1.0e-5;
			double xtol			= 1.0e-16;
			int[] errorFlag		= new int[1];
			int[] iprint		= new int[2];
			double[] x			= Array.clone(x0);
			double[] g			= new double[nDim];		//x, g, errorFlag and diag are used to pass objects by reference between the different functions
			double[] hk0		= new double[nDim];		//so they can all alter the same value/memory location rather than being passed by value (see errorFlag)

			iprint [0]	= 0;
			iprint [1]	= 0;
			errorFlag[0]=0;


			do {										//Loop until convergence (errorFlag=0) or if iterations exceed 200 
				LBFGSFunctions.CompactOutput output = f.evaluate(x);
				functionValue = output.value;
				g = output.gradient;
				iterations += 1;
				
				try {
					LBFGSOldMinimizer.lbfgs (nDim, memDepth, x, functionValue, g, provideHk0, hk0, iprint, epsilon, xtol, errorFlag);
				} catch (LBFGSOldMinimizer.ExceptionWithIflag e) {
					System.err.println( "Sdrive: lbfgs failed.\n"+e );
				}
			} while (errorFlag[0]!=0 && iterations<=400 );
			if (iterations <= 400) {
				oldCalls += iterations;
				oldLoops++;
				oldFVal  += functionValue;
				oldGNorm += Array.norm(g);
			}
			
			try {
				testMinimizerMC.minimize(x0);
				newCalls += testMinimizerMC.nFunctionEvals;
				newIters += testMinimizerMC.reportIters;
				newFVal	 += testMinimizerMC.reportFVal;
				newGNorm += testMinimizerMC.reportGNorm;
				newLoops++;
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			try {
				testMinimizer.minimize(x0);
				simCalls += testMinimizer.nFunctionEvals;
				simIters += testMinimizer.reportIters;
				simFVal  += testMinimizer.reportFVal;
				simGNorm += testMinimizer.reportGNorm;
				simLoops++;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		System.out.println(oldCalls/oldLoops+"\t"+oldLoops+"\t"+oldCalls/oldLoops+"\t"+oldFVal/oldLoops+"\t"+oldGNorm/oldLoops);
		System.out.println(newCalls/newLoops+"\t"+newLoops+"\t"+newIters/newLoops+"\t"+newFVal/newLoops+"\t"+newGNorm/newLoops);
		System.out.println(simCalls/simLoops+"\t"+simLoops+"\t"+simIters/simLoops+"\t"+simFVal/simLoops+"\t"+simGNorm/simLoops);
		System.out.println(testMinimizerMC.counter);
		System.out.println(testMinimizerMC.counter2);
	}
}
