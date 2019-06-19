package base;

import java.util.ArrayList;
import java.util.Random;
import Jama.*;
import test.SVM;

public class sLBFGS_Test extends Minimizer{
	private boolean isVerbose, randomSelect = false;
	private int d, N, b, bH, M, m, L, currDepth, maxEpoch;
	private double eta, delta, fdHVPStepSize = 5E-2, gradientNormBound = 100;
	private double[] rho;
	private double[][] s, y, U, UT;
	private Fit fitOutput;
	private Model.CompactGradientOutput fOut;
	
	public boolean useReducedSpace = false;
	public int svrgSubBatch = 1;
	
	//LBFGS object constructor; load basic minimization parameters. To be used for all subsequent minimizations using this object.
	public sLBFGS_Test(Model model, int gradientBatch, int hessianBatch, int memorySize, 
			int maxEpoch, int hessianPeriod, int epochIterations, double stepsize,
			double epsilon, double delta, boolean isVerbose) {
		this.model		= model;
		d				= model.getNFeatures();
		N				= model.N;				// Number of data points 
		if (gradientBatch<1) throw new IllegalArgumentException("Gradient batch size must be greater than 1!");
		b				= gradientBatch;		// Batch size for stochastic gradient updates
		if (hessianBatch<1) throw new IllegalArgumentException("Hessian batch size must be greater than 1!");
		bH				= hessianBatch;			// Batch size for stochastic Hessian updates
		if (memorySize<=0)	throw new IllegalArgumentException("Memory size must be positive!");
		M				= memorySize;			// Maximum memory depth
		if (maxEpoch<1) throw new IllegalArgumentException("The maximum number of epochs must be positive!");
		this.maxEpoch	= maxEpoch;				// The maximum number of epochs before termination
		if (hessianPeriod<2) throw new IllegalArgumentException("The hessian update period must be larger than 1!");
		L				= hessianPeriod;		// Number of iterations before the inverse hessian is updated
		if (epochIterations<=hessianPeriod) throw new IllegalArgumentException("The number of iterations per epoch must be greater than the hessian period!");
		m				= epochIterations;		// Number of iterations to run per epoch
		if (stepsize<=0) throw new IllegalArgumentException("Eta (step size) must be greater than 0!");
		eta				= stepsize;				// The fixed step size value
		if (epsilon<0) throw new IllegalArgumentException("Epsilon cannot be negative!");
		this.epsilon	= epsilon;				// The accuracy with which the solution needs to be found
		if (delta<0) throw new IllegalArgumentException("Delta (inverse hessian regularization parameter) cannot be negative!");
		this.delta		= delta;				// An optional inverse hessian regularization parameter
		this.isVerbose	= isVerbose;
	}
	
	public Fit doMinimize(double[] seed, String trajectoryFile) throws Exception {
		// Construct transformation matrices if we are to use reduced space
		if (useReducedSpace) {
			computeHessian(seed);
		} else {
			U = null;
			UT= null;
		}		
		int r = 0;									// Number of currently computed Hessian correction pairs
		currDepth = 0;
		double tStart, egNorm, divisor;
		// Full gradient, variance reduced gradient, and components of variance reduced gradient
		double[] mu_k = new double[d], v_t = new double[d], grad_f_xt = new double[d], grad_f_wk = new double[d];
		// Effective gradient
		double[] effGrad = new double[d];
		// Positions in current iterations
		double[] w_k = new double[d], w_k_prev = new double[d], x_t = new double[d];
		// Average of path traveled in the current and previous inverse hessian updates
		double[] u_r = new double[d], u_r_prev = new double[d];
		// Components of two-loop update
		double[] s_r = new double[d], y_r= new double[d];
		rho = new double[M];
		s = new double[M][d];
		y = new double[M][d];
		// Stores the history of all previous steps in the current epoch
		ArrayList<double[]> x_t_hist = new ArrayList<double[]>();
		Random generator = new Random();
		
		// Deal with a potential seed and initialize
		if (seed!=null) {
			try {
				w_k = Array.clone(compress(seed));
			} catch (Exception e) {
				throw new IllegalArgumentException("Improper seed!");
			}
		} else {
			w_k = Array.scalarMultiply(Array.randomDouble(d), randomSeedScale);
		}
		fitOutput = model.generateFit(seed);
		w_k_prev = Array.clone(w_k);
			
		// Init the number of loops over data points
		model.evaluatedDataPoints = 0;
		tStart	= System.nanoTime();
		// Loop over epochs 
		for (int k=0; k<=maxEpoch; k++) {
			fitOutput.addStep(uncompress(w_k));
			// Compute full gradient for variance reduction
			if (svrgSubBatch<2) {
				fOut = evaluate(uncompress(w_k));
			} else {
				//model.setSuperBatch(N/svrgSubBatch);
				model.sampleBatch(N/svrgSubBatch);
				fOut = stochasticEvaluate(uncompress(w_k));
			}
			// Print some information about the current epoch
			if (isVerbose) {
				if (k==0) {
					System.out.println("Starting Function Value: "+fOut.functionValue);
					System.out.println("    Epochs   Data Loops           Likelihood       Distance Moved"+
							"        Gradient Norm");
				} else {
					printStep(k, model.evaluatedDataPoints/N, fOut.functionValue, Array.norm(Array.subtract(w_k, w_k_prev)), 
							Array.norm(compress(fOut.gradientVector)));					
				}
			}

			// Check for convergence, final epoch
			if (Double.isNaN(fOut.functionValue) || Double.isInfinite(fOut.functionValue)) {
				throw new Exception("sLBFGS Failure: NaN encountered! Try reducing the step size...");
			}
			if (Array.norm(compress(fOut.gradientVector))/Math.max(1, Array.norm(w_k)) < epsilon) {
				fitOutput.recordFit(k, -1, (System.nanoTime()-tStart)/1E9, fOut.functionValue, model);
				System.out.println("Convergence criteria met.");
				return fitOutput;
			}
			if (k==maxEpoch) {
				break;			// This allows convergence to be tested on the FINAL epoch iteration without computation
			}

			mu_k = Array.clone(compress(fOut.gradientVector));		// Assign variance reduced gradient
			x_t = Array.clone(w_k);							// Set x_t to current value of w_k
						
			// Perform m stochastic iterations before a full gradient computation takes place
			for (int t=1; t<=m; t++) {
				// Compute the current stochastic gradient estimate; begin by sampling a minibatch
				model.sampleBatch(b);
				// Next, compute the reduced variance gradient
				grad_f_xt = compress(stochasticEvaluate(uncompress(x_t)).gradientVector);
				grad_f_wk = compress(stochasticEvaluate(uncompress(w_k)).gradientVector);
				v_t = Array.subtract(Array.add(grad_f_xt, mu_k), grad_f_wk);
				
				// Update u_r with current position
				u_r = Array.add(u_r, x_t);			
				// Compute next iteration step; condition the gradient so as not to produce rapid oscilations (extreme function values)
				x_t_hist.add(Array.clone(x_t));		// Need to store the history of iteration steps				
				if (r < 1) {						// Until a single hessian correction has taken place, H_0 = I
					effGrad = v_t;
				} else {							// Compute the two-loop recursion product
					effGrad = twoLoopRecursion(v_t);
				}			
				// Bound the effective gradient update step
				egNorm = Array.norm(effGrad);
				divisor = 1;
				while (egNorm/divisor > gradientNormBound) {
					divisor *= 10;
				}
				x_t = Array.addScalarMultiply(x_t, -eta, Array.scalarMultiply(effGrad, 1/divisor));
				
				// Check to see if L iterations have passed (triggers hessian update)
				if (t % L == 0) {
					// Increment the number of hessian correction pairs
					r++;
					// Finish computing u_r
					u_r = Array.scalarMultiply(u_r, 1.0/((double) L));
					
					// Use HVP to compute hessian updates. Begin by sampling a minibatch
					model.sampleBatch(bH);
					// Compute s_r update
					s_r = Array.subtract(u_r, u_r_prev);
					// Compute y_r estimate using HVP
					y_r = compress(Array.subtract(stochasticEvaluate(uncompress(Array.addScalarMultiply(u_r, fdHVPStepSize, s_r))).gradientVector, 
							stochasticEvaluate(uncompress(Array.addScalarMultiply(u_r, -fdHVPStepSize, s_r))).gradientVector));
					y_r = Array.scalarMultiply(y_r, 1.0/(2*fdHVPStepSize));
					// Store latest values of s_r and y_r
					add(s_r, y_r);
					
					// Resetting u_r for next evaluation
			        u_r_prev = Array.clone(u_r);
			        u_r = new double[d];
				}
			}
			// Choose either the last position vector x_t or a random one from the previous epoch as the starting point for the next
			if (randomSelect) {
				w_k = Array.clone(x_t_hist.get(0));
				for (int currPos=1; currPos<x_t_hist.size(); currPos++) {
					w_k = Array.add(w_k, x_t_hist.get(currPos));
				}
				w_k = Array.scalarMultiply(w_k, 1.0/((double) x_t_hist.size()));
//				w_k = Array.clone(x_t_hist.get(generator.nextInt(m)));
			} else {
				w_k = Array.clone(x_t);
			}
			x_t_hist = new ArrayList<double[]>();
		}
		throw new Exception("sLBFGS Failure: maximum epochs exceeded without convergence!");
	}
	
	private double[] twoLoopRecursion(double[] v_t) {
		double alpha, beta, gamma;
		double[] q = new double[d], r = new double[d];
		double[] alphas = new double[currDepth];
		
		// Begin by cloning the input gradient
		q = Array.clone(v_t);
		
		// The first loop (starts from the latest entry and goes to the earliest)
		for (int i=0; i<currDepth; i++) {
			// Compute and store alpha_i = rho_u*s_i*q
			alpha = rho[i]*Array.dotProduct(s[i], q);
			alphas[i] = alpha;
			// Update q: q = q - alpha_i*y_i
			q = Array.addScalarMultiply(q, -alpha, y[i]);
		}		
		// Start computing R. To do so, begin by computing gamma_k = s_k*y_k/(y_k*y_k)
		gamma = Array.dotProduct(s[currDepth-1], y[currDepth-1])/
				Array.dotProduct(y[currDepth-1], y[currDepth-1]);
		// r = gamma_k*q/(1 + delta*gamma_k); the denominator includes the pseudo-hessian 
		// regularization parameter delta NOTE: There is no need to multiply by I here, 
		// as that will anyway produce a dot product 
		r = Array.scalarMultiply(q, gamma/(1.0 + delta*gamma));
		
		// Second loop (goes in reverse, starting from the earliest entry)
		for (int i=currDepth-1; i>=0; i--) {
			// beta = rho_i*y_i*r
			beta = rho[i]*Array.dotProduct(y[i], r);
		    // r = r + s_i*(alpha_i-beta)
			r = Array.addScalarMultiply(r, alphas[i]-beta, s[i]);
		}
		return r;
	}
	
	private void add(double[] inS, double[] inY) {
		// Compute rho in the lbfgs two-loop method: rho_j = 1/s_j^T*y_j
		cycleDown(rho, 1.0/Array.dotProduct(inS, inY));
		// Add values to structure
		cycleDown(s, inS);
		cycleDown(y, inY);
		currDepth = Math.min(currDepth+1, M);
	}
	
	// Cycle rows in matrix downwards and add a new row on top
	private void cycleDown(double[][] input, double[] newRow) {	
		for (int i=input.length-1; i>0; i--) {
			for (int j=0; j<input[0].length; j++) {
				input[i][j] = input[i-1][j];
			}
		}
		for (int i=0; i<input[0].length; i++) {
			input[0][i] = newRow[i];
		}
	}
	
	private void cycleDown(double[] input, double newValue) {	
		for (int i=input.length-1; i>0; i--) {
			input[i] = input[i-1];
		}
		input[0] = newValue;
	}
	
	private void computeHessian(double[] seed) {
		int trueDims = 0;
		double cutoff = 1E-10;
		double[] baseBetas, modBetas, diagS;
		double[][] forwardDifference= new double[d][d];
		double[][] reverseDifference= new double[d][d];
		double[][] fdHessian		= new double[d][d];
		double[][] fullU;
		SingularValueDecomposition svd;
		
		// First, deal with the starting point of the hessian computation
		if (seed!=null) {
			try {
				baseBetas = Array.clone(seed);
			} catch (Exception e) {
				throw new IllegalArgumentException("Improper seed!");
			}
		} else {
			baseBetas = Array.scalarMultiply(Array.randomDouble(d), randomSeedScale);
		}
		
		// Next, construct the fdHessian by computing the gradients
		for (int i=0; i<d; i++) {
			try {
				modBetas			= Array.clone(baseBetas);
				modBetas[i]			+= fdHVPStepSize;
				forwardDifference[i]= evaluate(modBetas).gradientVector;
				modBetas[i]			-= 2*fdHVPStepSize;
				reverseDifference[i]= evaluate(modBetas).gradientVector;
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		// Find symmetrized finite differences (forward/backward FD)
		for (int i=0; i<d; i++) {
			for (int j=0; j<d; j++) {
				fdHessian[i][j] = (forwardDifference[j][i]-reverseDifference[j][i])/(4*fdHVPStepSize) + 
									(forwardDifference[i][j]-reverseDifference[i][j])/(4*fdHVPStepSize);
			}
		}
		
		// Now compute SVD.
		svd = new SingularValueDecomposition(new Matrix(fdHessian));
		
		// Extract column space. First find true dimensionality
		diagS = svd.getSingularValues();
		for (int i=0; i<diagS.length; i++) {
			if (diagS[i]-((SVM) model).lambda>cutoff) {
				trueDims++;
			}
		}
		// Now extract column space
		fullU = svd.getU().getArrayCopy();
		U = new double[d][trueDims];
		for (int col = 0; col<trueDims; col++) {
			for (int i = 0; i<d; i++) {
				U[i][col] = fullU[i][col];
			}
		}
		UT = Array.transpose(U);
		d = trueDims;
	}
	
	private double[] compress(double[] input) {
		if (useReducedSpace) {
			return Array.matrixVectorMultiply(UT, input);
		} else {
			return input;
		}
	}
	
	private double[] uncompress(double[] input) {
		if (useReducedSpace) {
			return Array.matrixVectorMultiply(U, input);
		} else {
			return input;
		}
	}
}