package minimizers;

import java.util.ArrayList;
import java.util.Formatter;

import base.*;

public class Minimizer {
	protected double epsilon, randomSeedScale = .05;
	protected Model model;
	
	/* This is the function where the actual minimizer algorithm must be
	 * defined. The public functions below are wrapper functions that
	 * automatically provide the user with added functionality (see below).
	 */
	protected Fit doMinimize(double[] seed, String trajectoryFile) 
			throws Exception {
		return null;
	}
	
	//Special case of minimize to remove negative eigenvalues
	public Fit[] minimize(double[] seed, String trajectoryFile) 
			throws Exception {
		int retryCount = 0;
		double originalSeedScale = randomSeedScale;
		Fit fitOutput;
		
		/* Try reseeding the fit if no seed is provided, as sometimes the random
		 * start position can make it hard for LBFGS to converge. This will not 
		 * trigger for Pattern Search. This refitting step will be attempted 
		 * maxRetries number of times. Additionally, the random seed magnitude
		 * will drop by a factor of 10 every failure. 
		 */
		if (seed==null) {
			while (true) {
				try {
					fitOutput = doMinimize(seed, trajectoryFile);
					break;
				} catch (ArithmeticException e) {
					if (retryCount>model.maxFitRetries) {
						throw new Exception("Random reseed failed: General "
								+ "Minimizer Failure.", e);
					}
					System.out.println("Minimizer failed. Attempting refit with"
							+ " new random seed.");
					retryCount++;
					randomSeedScale /= 10;
					continue;
				} catch (Exception e) {
					throw e;
				}
			}
			randomSeedScale = originalSeedScale;
		} else {
			fitOutput = doMinimize(seed, trajectoryFile);
		}

		return new Fit[]{fitOutput};
	}
	
	public Fit[] minimize(double[] seed, String trajectoryFile, ArrayList<Object[]> datasets) throws Exception {
		Fit output = null;
		
		System.out.println("Fitting full dataset.");
		output = minimize(seed, trajectoryFile)[0];

		return new Fit[]{output};
	}
	
	/* gradientEval uncompresses the input (desymmetrizes) and compresses
	 * the gradient (symmetrizes) so that functions work in the compressed
	 * (reduced) space.
	 */
	protected Model.CompactGradientOutput gradientEval(double[] input, 
			boolean normalize) throws Exception {
		double[] tempGradient;
		Model.CompactGradientOutput out;
		
		model.setParams(input);
		out 				= model.gradientEval();
		//Handle the case where the gradient has not been defined
		if (out==null)	throw new UnsupportedOperationException("The function "
				+ "gradientEval() has not been defined!");
		tempGradient		= out.gradientVector;
		if (normalize) {
			out.functionValue	*= model.likelihoodNormalizer();
			tempGradient		= Array.scalarMultiply(tempGradient, model.likelihoodNormalizer());
		}
		out.gradientVector	= tempGradient;
		return out;	
	}
	
	protected void printStep(int iterations, int calls, double likelihood, 
			double distance, double ... params) {
		Formatter fmt = new Formatter();
		
		System.out.printf("   %7d      %7d   ", iterations, calls);
		fmt.format("%18.18s   %18.18s", String.format("%10.15f", likelihood), 
				String.format("%10.15f", distance));
		System.out.print(fmt);
		fmt = new Formatter();
		for (int i=0; i<params.length; i++) {
			fmt.format("   %18.18s", String.format("%10.15f", params[i]));
		}
		System.out.print(fmt+"\n");
		fmt.close();
	}
}
