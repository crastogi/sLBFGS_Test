package minimizers;

import base.*;

public class sLBFGS_Clean extends Minimizer{
	private boolean isVerbose;
	private int d, N, b, bH, M, m, L, currDepth, maxEpoch;
	private double eta, delta, fdHVPStepSize = 5E-2, gradientNormBound = 10;
	private double[] rho;
	private double[][] s, y;
	private Fit fitOutput;
	private Model.CompactGradientOutput fOut;
	
	//LBFGS object constructor; load basic minimization parameters. To be used for all subsequent minimizations using this object.
	public sLBFGS_Clean(Model model, int gradientBatch, int hessianBatch, int memorySize, 
			int maxEpoch, int hessianPeriod, int epochIterations, double stepsize,
			double epsilon, double delta, boolean isVerbose) {
		this.model		= model;
		d				= model.getNDimensions();
		N				= model.nCount;			// Number of data points 
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
		int r = 0;							// # of computed Hessian correction pairs
		currDepth = 0;
		double tStart, dirNorm, divisor, boundedDirNorm, dirNormMean = 0;
		double maxDirGrowth = 2.718282;		// Maximum allowed growth rate
		double dirDecayRate = .1;			// EMA decay rate
		double maxNormBound = Math.sqrt(d)*gradientNormBound;
		double[] mu_k;						// Full gradient
		double[] v_t;						// Variance reduced gradient
		double[] grad_f_xt, grad_f_wk;		// Variance reduced gradient components
		double[] dir;						// Effective gradient
		double[] w_k, w_k_prev, x_t;		// Positions in current iterations
		// Average of path traveled in the current and previous inverse hessian updates
		double[] u_r = new double[d], u_r_prev = new double[d];
		double[] s_r, y_r;					// Components of two-loop update
		rho = new double[M];
		s = new double[M][d];
		y = new double[M][d];
		
		// Deal with a potential seed and initialize
		if (seed!=null) {
			try {
				w_k = Array.clone(seed);
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
			fitOutput.addStep(w_k);
			// Compute full gradient for variance reduction
			fOut = gradientEval(w_k, true);
			if (k==0) {
				dirNormMean= Array.norm(fOut.gradientVector);
			}
			
			// Print some information about the current epoch
			if (isVerbose) {
				if (k==0) {
					System.out.println("Starting Function Value: "+fOut.functionValue);
					System.out.println("    Epochs   Data Loops           Likelihood       Distance Moved"+
							"        Gradient Norm");
				} else {
					if (xStar!=null) {
						printStep(k, model.evaluatedDataPoints/N, fOut.functionValue, Array.dist(w_k, w_k_prev), 
								Array.norm(fOut.gradientVector), Array.dist(w_k, xStar));	
					} else {
						printStep(k, model.evaluatedDataPoints/N, fOut.functionValue, Array.dist(w_k, w_k_prev), 
								Array.norm(fOut.gradientVector));	
					}
				}
			}

			// Check for convergence, final epoch
			if (Double.isNaN(fOut.functionValue) || Double.isInfinite(fOut.functionValue)) {
				throw new Exception("sLBFGS Failure: NaN encountered! Try reducing the step size...");
			}
			if (Array.norm(fOut.gradientVector)/Math.max(1, Array.norm(w_k)) < epsilon) {
				fitOutput.recordFit(k, -1, (System.nanoTime()-tStart)/1E9, fOut.functionValue, model);
				System.out.println("Convergence criteria met.");
				model.hessianEval();
				fitOutput.storeHessian(model.getHessian());
				return fitOutput;
			}
			if (k==maxEpoch) {
				break;			// This allows convergence to be tested on the FINAL epoch iteration without computation
			}

			mu_k = Array.clone(fOut.gradientVector);		// Assign variance reduced gradient
			x_t = Array.clone(w_k);							// Set x_t to current value of w_k
			w_k_prev = Array.clone(w_k);					
						
			// Perform m stochastic iterations before a full gradient computation takes place
			for (int t=1; t<=m; t++) {
				// Compute the current stochastic gradient estimate; begin by sampling a minibatch
				model.sampleBatch(b);
				// Next, compute the reduced variance gradient
				grad_f_xt = stochasticEvaluate(x_t).gradientVector;
				grad_f_wk = stochasticEvaluate(w_k).gradientVector;
				v_t = Array.subtract(Array.add(grad_f_xt, mu_k), grad_f_wk);
				
				// Update u_r with current position
				u_r = Array.add(u_r, x_t);			
				// Compute next iteration step; condition the gradient so as not to produce rapid oscilations (extreme function values)				
				if (r < 1) {						// Until a single hessian correction has taken place, H_0 = I
					dir = v_t;
				} else {							// Compute the two-loop recursion product
					dir = twoLoopRecursion(v_t);
				}			
				// Bound the effective gradient update step
				dirNorm = Array.norm(dir);
				
				// See if dirNorm is spiking
				if (dirNorm > maxDirGrowth*dirNormMean && r>1) {
					boundedDirNorm = dirNormMean;
				} else {
					boundedDirNorm = dirNorm;
				}
				// Bound gradient update
				boundedDirNorm = Math.min(boundedDirNorm, maxNormBound);
				// Compute exponential moving average
				dirNormMean = dirDecayRate*boundedDirNorm+(1-dirDecayRate)*dirNormMean;
				divisor = dirNorm/dirNormMean;
				// Compute next step
				x_t = Array.addScalarMultiply(x_t, -eta/divisor, dir);
								
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
					y_r = Array.subtract(stochasticEvaluate(Array.addScalarMultiply(u_r, fdHVPStepSize, s_r)).gradientVector, 
							stochasticEvaluate(Array.addScalarMultiply(u_r, -fdHVPStepSize, s_r)).gradientVector);
					y_r = Array.scalarMultiply(y_r, 1.0/(2*fdHVPStepSize));
					// Store latest values of s_r and y_r
					add(s_r, y_r);
										
					// Resetting u_r for next evaluation
			        u_r_prev = Array.clone(u_r);
			        u_r = new double[d];
				}
			}
			w_k = Array.clone(x_t);
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
}