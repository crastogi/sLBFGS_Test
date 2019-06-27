package base;

public class SGD extends Minimizer{
	private boolean isVerbose, useSVRG, useADAM, useAvg = true;
	private int d, N, bs, M, eIters;
	private double eta, epsilon, gradientNormBound = 100;
	private Fit fitOutput;
	
	// NOTE: ADAM parameters below
	private double b1 = .9, b2 = .999, e = 1E-8;
	
	// SGD object constructor; load basic minimization parameters. To be used 
	// for all subsequent minimizations using this object.
	public SGD(Model model, boolean isStochastic, int gradientBatch, int maxEpoch, 
			boolean singleEpochLoop, int epochIterations, double stepsize, 
			double epsilon, boolean useSVRG, boolean useADAM, boolean isVerbose) {
		this.model		= model;
		d				= model.getNFeatures();	// Dimensionality
		N				= model.N;				// Number of data points 
		if (isStochastic) {
			bs = 1;
		} else {
			if (gradientBatch<1) {
				throw new IllegalArgumentException("Gradient batchsize must be >1!");
			}
			bs = gradientBatch;
		}
		// The maximum number of epochs before termination
		if (maxEpoch<1) throw new IllegalArgumentException("Max epochs must be >1!");
		M				= maxEpoch;
		// Set the number of iterations per epoch
		if (singleEpochLoop) {
			// Make it so that the dataset is looped once per epoch
			eIters = (int) Math.floor(((double) N)/((double) bs));
 		} else {
 			if (epochIterations<1) {
 				throw new IllegalArgumentException("Epoch iterations must be > 1!");
 			}
 			eIters = epochIterations;
 		}
		// SGD stepsize
		if (stepsize<=0) {
			throw new IllegalArgumentException("Step size must be positive!");
		}
		eta = stepsize;
		// The accuracy with which the solution needs to be found
		if (epsilon<0) {
			throw new IllegalArgumentException("Convergence criteria must be positive!");
		}
		this.epsilon	= epsilon;
		this.useSVRG	= useSVRG;			// Should SVRG be used
		this.useADAM	= useADAM;			// Should ADAM be used
		this.isVerbose	= isVerbose;		// Verbosity
	}
	
	public Fit doMinimize(double[] seed, String trajectoryFile) throws Exception {
		int totIters = 0;
		double currFVal, a_t, v_tSum, gNorm, divisor;
		double[] x0, x_t, g_t, g_k, w_k, m_t, v_t, x_t_prev;
		double[] avg_g = null, avg_x = null, mu_k = null;
		Model.CompactGradientOutput out;
		
		/*Check to see if seed is of correct length and symmetrize seed, if need
		 * be. If no seed is provided, use a random start point to break 
		 * symmetry. */
		if (seed!=null) {
			x0 = Array.clone(seed);
		} else {
			x0 = Array.scalarMultiply(Array.randomDouble(d), .05);
		}
		fitOutput = model.generateFit(seed);

		// Initialize Parameters, even if they are not needed
		x_t		= Array.clone(x0);
		x_t_prev= Array.clone(x0);
		g_t 	= new double[d];
		w_k 	= Array.clone(x0);
		m_t 	= new double[d];
		v_t 	= new double[d];
		
		// Print reporting line
		if (isVerbose) {
			if (useSVRG) {
				System.out.print("    Epochs   Data Loops           "+
						"Likelihood       Distance Moved        Gradient Norm");
			} else {
				System.out.println("    Epochs   Data Loops       "+
						"Distance Moved        Gradient Norm");
			}
			if (xStar!=null) {
				System.out.println("       Distance to x*");
			} else {
				System.out.println("");
			}
		}
		
		// Loop over epochs
		for (int m=0; m<=M; m++) {
			// Perform SVRG update, if being used
			if (useSVRG) {
				out 		= evaluate(w_k);
				currFVal	= out.functionValue;
				mu_k 		= Array.clone(out.gradientVector);
				if (isVerbose) {
					if (xStar!=null) {
						printStep(m, model.evaluatedDataPoints/N, currFVal,
								Array.norm(Array.subtract(w_k, x_t_prev)), 
								Array.norm(mu_k), Array.norm(Array.subtract(w_k, xStar)));
					} else {
						printStep(m, model.evaluatedDataPoints/N, currFVal,
								Array.norm(Array.subtract(w_k, x_t_prev)), 
						Array.norm(mu_k));
					}
				}
				x_t_prev = w_k;
				// Can use this SVRG update step to declare convergence
				if (Array.norm(mu_k)/Math.max(1, Array.norm(w_k)) < epsilon) {
					fitOutput.recordFit(m, -1, -1, currFVal, model);
					System.out.println("Convergence criteria met.");
					return fitOutput;
				}
			} else if (m>0) {		// Need to ensure one epoch has passed
				// Use averaged gradients and function values to test convergence
				avg_x = Array.scalarMultiply(avg_x, 1.0/((double) eIters));
				avg_g = Array.scalarMultiply(avg_g, 1.0/((double) eIters));
				if (isVerbose) {
					if (xStar!=null) {
						printStep(m, model.evaluatedDataPoints/N, 0,
								Array.norm(Array.subtract(avg_x, x_t_prev)), 
								Array.norm(avg_g), Array.norm(Array.subtract(x_t, xStar)));
					} else {
						printStep(m, model.evaluatedDataPoints/N, 0,
								Array.norm(Array.subtract(avg_x, x_t_prev)), 
								Array.norm(avg_g));
					}
				}
				x_t_prev = avg_x;
				if (Array.norm(avg_g)/Math.max(1, Array.norm(avg_x)) < epsilon) {
					// Ensure convergence using a full batch evaluation
					out 	= evaluate(x_t);
					if (Array.norm(out.gradientVector)/Math.max(1, Array.norm(x_t)) < epsilon) {
						fitOutput.recordFit(m, -1, -1, out.functionValue, model);
						System.out.println("Convergence criteria met.");
						return fitOutput;
					}
				}
			}
			
			// Allow to check convergence until last epoch, otherwise terminate
			if (m==M) {
				break;
			}
			
			// Init vectors used to store average position and gradient values
			avg_x = new double[d];
			avg_g = new double[d];
			
			// Loop over epoch iterates
			for (int t=0; t<eIters; t++) {
				// Update average stores
				for (int currIdx=0; currIdx<d; currIdx++) {
					avg_x[currIdx] += x_t[currIdx];
					avg_g[currIdx] += g_t[currIdx];
				}
				totIters++;						// For ADAM update
				model.sampleBatch(bs);			// Set batchsize
				// Compute grad @x_t
				g_t = stochasticEvaluate(x_t).gradientVector;
				if (useSVRG) {					// SVRG update
					g_k = stochasticEvaluate(w_k).gradientVector;
					// g_t = g_t - g_k + mu_k
					for (int currIdx=0; currIdx<d; currIdx++) {
						g_t[currIdx] += -g_k[currIdx] + mu_k[currIdx];
					}
				}
				// Condition gradient to stabilize optimization
				gNorm = Array.norm(g_t);
				divisor = 1;
				while (gNorm/divisor > gradientNormBound) {
					divisor *= 10;
				}
				
				// Adjust step if using ADAM
				if (useADAM) {
					// Unwrapped array functions
					v_tSum = 0;
					for (int currIdx=0; currIdx<d; currIdx++) {
						// m_t = b1*m_t + (1-b1)*g_t;
						m_t[currIdx] = b1*m_t[currIdx] + (1-b1)*g_t[currIdx];
						// v_t = b2*v_t + (1-b2)*g_t.^2
						v_t[currIdx] = b2*v_t[currIdx] + (1-b2)*g_t[currIdx]*g_t[currIdx];
						v_tSum += v_t[currIdx];
					}
					// Use alternate expansion: alpha_t = alpha*sqrt(1-b2^t)/[(1-b1^t)*(sqrt(sum(vt))+e)]
					a_t = eta*Math.sqrt(1-Math.pow(b2, totIters))/((1-Math.pow(b1, totIters))*(Math.sqrt(v_tSum)+e));
					// Step update: x_t = x_t - a_t*m_t
					x_t = Array.addScalarMultiply(x_t, -a_t/divisor, m_t);
				} else {
					// Use typical update step
					x_t = Array.addScalarMultiply(x_t, -eta/divisor, g_t);
				}				
			}
			// Update w_k if SVRG is being used
			if (useSVRG) {
				if (useAvg) {
					w_k = Array.scalarMultiply(avg_x, 1.0/((double) eIters));
				} else {
					w_k = x_t;
				}
			}
		}
		throw new Exception("SGD Failure: maximum epochs exceeded without convergence!");
	}
}