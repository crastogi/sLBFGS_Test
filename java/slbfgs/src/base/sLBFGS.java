package base;

import java.util.ArrayList;

public class sLBFGS extends Minimizer{
	private boolean isVerbose;
	private int d, N, b, bH, M, m, L, maxEpoch;
	private double eta, delta;
	private Fit fitOutput;
	private Model.CompactGradientOutput fOut;
	
	private double IHTime = 0, TwoLoopTime = 0;

	//TODO: Add convergence criteria cutoff!
	//TODO: Include option to toggle randomized selection of final point!
	//TODO: Test the total number of passes over data points!
	
	//LBFGS object constructor; load basic minimization parameters. To be used for all subsequent minimizations using this object.
	public sLBFGS(Model model, int gradientBatch, int hessianBatch, int memorySize, 
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
		double tEpochStart = 0, tEpochEnd = 0;
		
		int r = 0;									// Number of currently computed Hessian correction pairs
		double tStart, tEnd;
		double[] v_t, v_t_prev;						// Current and previous iteration's variance reduced gradient estimate
		double[] u_r, u_r_prev, s_r;				// Average of path travelled in the current and previous inverse hessian updates
		double[] mu_k;								// Full gradient for variance-reduced Gradient
		double[] w_k, w_k_prev;						// Position in current epoch iteration
		double[] x_t;								// Position in the current iteration
		double[] y_r;
		double[] grad_f_xt, grad_f_wk;				// Components of the variance reduced gradient v_t
		InverseHessian IH = new InverseHessian();	// Initialize a new inverse hessian object
		ArrayList<double[]> x_t_hist;				// Stores the history of all previous steps in the current epoch
		
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
		v_t = new double[d];
		u_r = new double[d];
		u_r_prev = new double[d];
//		x_t_hist = new ArrayList<double[]>();
		
		// Loop over epochs 
		System.out.println("Epochs   Data Loops           Likelihood"
				+ "       Distance Moved");
		tStart	= System.nanoTime();
		for (int k=0; k<maxEpoch; k++) {
			// Compute full gradient for variance reduction
			fOut = evaluate(w_k);
			// Check for convergence
			if (Array.norm(fOut.gradientVector)/Math.max(1, Array.norm(w_k)) < epsilon) {
				System.out.println("converged");
				fitOutput.recordFit(0, 0, tStart, fOut.functionValue, model);
				System.out.println("Total time consumed: "+(System.nanoTime()-tStart)/1E9);
				System.out.println("IH Time consumed: "+IHTime);
				System.out.println("Two Loop Time consumed: "+TwoLoopTime);
				System.out.println("Batch Sampling Time consumed: "+model.SampleTime);
				System.out.println("Stochastic Gradient Time consumed: "+model.StochasticTime);
				return fitOutput;
			}
			
			//TODO: Fix printing out of update steps
			// Print some information about the current epoch
			if (isVerbose) {
				printStep(k, model.evaluatedDataPoints/N, fOut.functionValue, Array.norm(Array.subtract(w_k, w_k_prev)), (tEpochEnd-tEpochStart)/1E9);
			}

			tEpochStart = System.nanoTime();
			
			// Assign variance reduced gradient
			mu_k = Array.clone(fOut.gradientVector);
			// Set x_t to current value of w_k
			x_t = Array.clone(w_k);
			
			// Perform m stochastic iterations before a full gradient computation takes place
			for (int t=1; t<=m; t++) {
				// Compute the current stochastic gradient estimate.
				v_t_prev = Array.clone(v_t);
				// Begin by sampling a minibatch
				model.sampleBatch(b);
				// Next, compute the reduced variance gradient
				grad_f_xt = stochasticEvaluate(x_t).gradientVector;
				grad_f_wk = stochasticEvaluate(w_k).gradientVector;
				v_t = Array.subtract(Array.add(grad_f_xt, mu_k), grad_f_wk);
				// Update u_r with current position
				u_r = Array.add(u_r, x_t);
				
				// Compute next iteration step
//				x_t_hist.add(Array.clone(x_t));		// Need to store the history of iteration steps
				if (r < 1) {						// Until a single hessian correction has taken place, H_0 = I
					x_t = Array.addScalarMultiply(x_t, -eta, v_t);
				} else {							// Compute the two-loop recursion product
					x_t = Array.addScalarMultiply(x_t, -eta, twoLoopRecursion(IH, v_t));
				}
				
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
					//TODO: turn fd movement step here into a tunable parameter
					y_r = Array.subtract(stochasticEvaluate(Array.addScalarMultiply(u_r, 1e-2, s_r)).gradientVector, 
							stochasticEvaluate(Array.addScalarMultiply(u_r, -1e-2, s_r)).gradientVector);
					y_r = Array.scalarMultiply(y_r, 1.0/(2*1e-3));
					// Store latest values of s_r and y_r
					IH.add(s_r, y_r);
					
					// Resetting u_r for next evaluation
			        u_r_prev = Array.clone(u_r);
			        u_r = new double[d];
				}
			}
			//TODO: provide option for selecting between this and random selection
			w_k = Array.clone(x_t);
//			x_t_hist = new ArrayList<double[]>();
			tEpochEnd = System.nanoTime();
		}
		//TODO: Find way to complain about not converging if all epochs are exhausted!
		
		// Compute full gradient for variance reduction
		fOut = evaluate(w_k);
		//TODO: Fix printing out of update steps
		// Print some information about the current epoch
		if (isVerbose) {
			printStep(maxEpoch, model.evaluatedDataPoints/N, fOut.functionValue, Array.norm(Array.subtract(w_k, w_k_prev)));
		}
		fitOutput.recordFit(0, 0, tStart, fOut.functionValue, model);
		System.out.println("Total time consumed: "+(System.nanoTime()-tStart)/1E9);
		System.out.println("IH Time consumed: "+IHTime);
		System.out.println("Two Loop Time consumed: "+TwoLoopTime);
		System.out.println("Batch Sampling Time consumed: "+model.SampleTime);
		System.out.println("Stochastic Gradient Time consumed: "+model.StochasticTime);
		return fitOutput;
	}
	
	private double[] twoLoopRecursion(InverseHessian IH, double[] v_t) {
		double tStart = System.nanoTime();
		
		int currMemDepth = IH.rho.size();
		double alpha, beta, gamma;
		double[] q, r;
		double[] alphas = new double[currMemDepth];
		
		// Begin by cloning the input gradient
		q = Array.clone(v_t);
		
		// The first loop (starts from the latest entry and goes to the first)
		for (int i=currMemDepth-1; i>=0; i--) {
			// Compute and store alpha_i = rho_u*s_i*q
			alpha = IH.rho.get(i)*Array.dotProduct(IH.s.get(i), q);
			alphas[i] = alpha;
			// Update q: q = q - alpha_i*y_i
			q = Array.addScalarMultiply(q, -alpha, IH.y.get(i));
		}
		
		// Start computing R. To do so, begin by computing gamma_k = s_k*y_k/(y_k*y_k)
		gamma = Array.dotProduct(IH.s.get(currMemDepth-1), IH.y.get(currMemDepth-1))/
				Array.dotProduct(IH.y.get(currMemDepth-1), IH.y.get(currMemDepth-1));
		// r = gamma_k*q/(1 + delta*gamma_k); the denominator includes the pseudo-hessian 
		// regularization parameter delta NOTE: There is no need to multiply by I here, 
		// as that will anyway produce a dot product 
		r = Array.scalarMultiply(q, gamma/(1.0 + delta*gamma));
		
		// Second loop (goes in reverse, starting from the first entry)
		for (int i=0; i<currMemDepth; i++) {
			// beta = rho_i*y_i*r
			beta = IH.rho.get(i)*Array.dotProduct(IH.y.get(i), r);
		    // r = r + s_i*(alpha_i-beta)
			r = Array.addScalarMultiply(r, alphas[i]-beta, IH.s.get(i));
		}
		
		TwoLoopTime += (System.nanoTime()-tStart)/1E9;
		
		return r;
	}
	
	private class InverseHessian {
		public ArrayList<Double> rho = new ArrayList<Double>();
		public ArrayList<double[]> s = new ArrayList<double[]>();
		public ArrayList<double[]> y = new ArrayList<double[]>();

		public InverseHessian() {
			// Do nothing
		}
		
		public void add(double[] inS, double[] inY) {
			double tStart = System.nanoTime();
			
			// Compute rho in the lbfgs two-loop method: rho_j = 1/s_j^T*y_j
			double currRho = 1.0/Array.dotProduct(inS, inY);
			
			// Check to see if memory depth has been reached.
			if (rho.size()==M) {
				// Memory depth reached. remove oldest (first) values
				rho.remove(0);
				s.remove(0);
				y.remove(0);
			}
			// Add values to structure
			rho.add(currRho);
			s.add(inS);
			y.add(inY);
			IHTime += (System.nanoTime()-tStart)/1E9;
		}
	}
}