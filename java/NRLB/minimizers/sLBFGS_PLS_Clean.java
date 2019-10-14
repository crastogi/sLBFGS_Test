package minimizers;

import base.*;

public class sLBFGS_PLS_Clean extends Minimizer{
	private boolean isVerbose, isBracketed;
	private int d, N, b, bH, M, m, L, currDepth, maxEpoch;
	private int maxMCSIterations	= 10;
	private int maxPLSIterations	= 10;
	private double eta, delta, fdHVPStepSize = 5E-2, dirNormBound = 10;
	private double fFull, effEtaMean;
    private double etaEMARate = 0.1;			// Eta EMA decay rate
    private double etaGrowthRate = 1.005;		// effEta growth rate
    private double maxEtaBound = 100;			// Maximum growth of eta
	private double dirEMARate = 0.1;			// Direction EMA decay rate
	private double maxDirGrowth = 2.718282;		// Maximum allowed dir growth rate
	// Parameters and constants for line search
	private double alphaMin, alphaMax, alphaL, fL, gL, fT, gT, alphaU, fU, gU;
	private double stepAlphaMin = 1e-20;
	private double stepAlphaMax = 1e20;
	private double uTol			= 1e-10;
	private double c1 			= .0001;
	private double c2 			= .999;
	private double WolfeThreshold= .1;
	private double[] rho, w_k, mu_k;
	private double[][] s, y;
	private Fit fitOutput;
	private Model.CompactGradientOutput fOut;
	private GP_Clean gp;						// Gaussian process used for evaluations
	
	//LBFGS object constructor; load basic minimization parameters. To be used for all subsequent minimizations using this object.
	public sLBFGS_PLS_Clean(Model model, int gradientBatch, int hessianBatch, int memorySize, 
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
		double maxDirBound = Math.sqrt(d)*dirNormBound;
		double[] dir;						// Effective gradient
		double[] w_k_prev, x_t;				// Positions in current iterations
		// Average of path traveled in the current and previous inverse hessian updates
		double[] u_r = new double[d], u_r_prev = new double[d];
		double[] s_r, y_r;                  // Components of two-loop update
		rho = new double[M];
		s = new double[M][d];
		y = new double[M][d];
		GP_Clean.CompactOutput plsOut = null;
		
		// Deal with a potential seed and initialize
		if (seed!=null) {
			try {
				w_k = Array.clone(seed);
			} catch (Exception e) {
				throw new IllegalArgumentException("Improper seed!");
			}
		} else {
			w_k		= Array.scalarMultiply(Array.randomDouble(d), randomSeedScale);
		}
		fitOutput	= model.generateFit(seed);
		w_k_prev	= Array.clone(w_k);
			
		// Init the number of loops over data points
		model.evaluatedDataPoints = 0;
		tStart	= System.nanoTime();
		// Loop over epochs 
		for (int k=0; k<=maxEpoch; k++) {
			fitOutput.addStep(w_k);
			// Compute full gradient for variance reduction
			fOut = gradientEval(w_k, true);
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

			fFull= fOut.functionValue;						// Assign variance reduced fVal
			mu_k = Array.clone(fOut.gradientVector);		// Assign variance reduced gradient
			x_t  = Array.clone(w_k);						// Set x_t to current value of w_k
			w_k_prev = Array.clone(w_k);					
			
			// Init update steps at the first iteration
            if (k==0) {
            	effEtaMean = eta;
                dirNormMean = Array.norm(mu_k);
                u_r = Array.add(u_r, x_t);
                plsOut = probLineSearch(x_t, fFull, mu_k, 0, Array.scalarMultiply(mu_k, 0), eta, Array.scalarMultiply(mu_k, -1.0/dirNormMean));
            }
			
			// Perform m stochastic iterations before a full gradient computation takes place
			for (int t=1; t<=m; t++) {
				// Check to see if L iterations have passed to trigger Hessian
				// update; need to do this NOW before NEXT pls call
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
				
				// Update u_r with current position
				u_r = Array.add(u_r, plsOut.position);			
				// Compute next iteration step; condition the gradient so as not to produce rapid oscilations (extreme function values)
				if (r < 1) {						// Until a single hessian correction has taken place, H_0 = I
					dir = plsOut.gradientVector;
				} else {							// Compute the two-loop recursion product
					dir = twoLoopRecursion(plsOut.gradientVector);
					// Check for NaNs in dir after 2loop recursion
					if (Double.isNaN(Array.sum(dir)) || Double.isInfinite(Array.sum(dir))) {
						throw new Exception("Optimizer Failure: NaN encountered " +
								"during two-loop recursion!");
					}
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
				boundedDirNorm = Math.min(boundedDirNorm, maxDirBound);
				// Compute exponential moving average
				dirNormMean = dirEMARate*boundedDirNorm+(1-dirEMARate)*dirNormMean;
				divisor = dirNorm/dirNormMean;
				// Compute next step using line search
				plsOut = probLineSearch(plsOut, Array.scalarMultiply(dir, -1.0/divisor));				
			}
			w_k = Array.clone(plsOut.position);
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
		// r = gamma_k*q/(1 + delta*gamma_k); the denominator includes the 
		// pseudo-hessian regularization parameter delta NOTE: There is no need 
		// to multiply by I here, as that will anyway produce a dot product 
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
	
	private GP_Clean.CompactOutput probLineSearch(GP_Clean.CompactOutput input, double[] dir) 
			throws Exception{
		return probLineSearch(input.position, input.functionValue, 
				input.gradientVector, input.varF, input.varDF, input.effEta, dir);
	}
		
	private GP_Clean.CompactOutput probLineSearch(double[] x0, double f0, double[] df0, 
			double var_f0, double[] var_df0, double alpha0, double[] dir) throws Exception{
		boolean wolfeCrit;
		int idx;
		double tt, beta, wolfeValue, gpMean, temp;
		
		// Create new GP and update at start
		beta		= Array.norm(df0);		// Compute scaling factor beta
		gp 			= new GP_Clean(c1, c2, f0, beta, alpha0, dir);
		tt			= 0;
	    gp.updateGP(tt, x0, f0, df0, var_f0, var_df0);
		
		// Set new start point and loop until budget is exhausted
		tt  		= 1.0;
		while (true) {
		    // Begin by evaluating function and updating GP at new tt value.
			evaluateFunction(tt);
		    
		    // Next, verify that the GP gradient at the origin is negative. If 
			// it is not, it means the GP believes the objective function is 
			// increasing, and the gradient direction is NOT a descent 
			// direction. Reset
		    if (gp.d1m(0) > 0) {
		        tt = eta/alpha0;
		        return makeOuts(tt, true);
		    }
		    
		    // Next, confirm if the new point satisfies the wolfe conditions
	    	// Loop over all evaluated positions, except the original point; use
	        // Wolfe probabilities and GP means to accept points. Ensure that 
		    // idx=0 is ignored in the following evaluation.
		    idx			= -1;					// Throw error just in case
		    gpMean		= Double.MAX_VALUE;
		    wolfeCrit	= false;
		    for (int currT=1; currT<gp.N; currT++) {
		    	wolfeValue		= gp.probWolfe(currT);
		    	if (wolfeValue > WolfeThreshold) {
		    		// Found a point that matches Wolfe criteria
		    		wolfeCrit	= true;
		    		// Check to see if the new gpMean is lower than the old one
		    		temp		= gp.m(currT);
		    		if (temp < gpMean) {
		    			idx		= currT;		// Set to new idx
		    			gpMean	= temp;
		    		}
		    	}
		    }
	        // If point(s) match Wolfe criteria, return optimal one
	        if (wolfeCrit) {
	            tt = gp.T[idx];
	            return makeOuts(tt, false); 
	        }
		    
		    // If not, check to see if the number of line search steps has been 
	        // exceeded; testing here allows the final point to be checked for 
	        // passing the Wolfe criteria
		    if (gp.N >= maxPLSIterations) {
		    	break;
		    }
		    
		    // None of the currently evaluated points are acceptable. Compute a 
		    // new test point using mcsearch on the gp.
		    try {
		        tt = mcsearch(Array.max(gp.T));
		        // Check to see if tt has already been tested (occurs when an
		        // existing point satisfies the non-probabilistic wolfe conditions)
		        if (gp.contains(tt)) {
		            // Check to see if tt is the same as the largest point tested
		        	if (tt==Array.max(gp.T)) {
		        		tt = Array.max(gp.T)*2;
		        	} else {
		                // New tt is the average of it and a point to the right
		        		idx = Array.matchIdx(gp.ordT, tt);
		                tt = (gp.ordT[idx]+gp.ordT[idx+1])/2;
		            }
		        }
		    } catch (Exception e) {
		    	// An error has taken place in the GP line search; reset
		    	tt = eta/alpha0;
		        return makeOuts(tt, true);
		    }
		}

		// Algorithm reached limit without finding acceptable point; choose 
		// point with lowest mean
        idx		= 1;
        gpMean	= Double.MAX_VALUE;		// Ignore idx = 0
        // Evaluate GP Means
        for (int currT=2; currT<gp.N; currT++) {
    		// Check to see if the new gpMean is lower than the old one
    		temp		= gp.m(currT);
    		if (temp < gpMean) {
    			idx		= currT;		// Set to new idx
    			gpMean	= temp;
    		}
        }
		tt = gp.T[idx];
        return makeOuts(tt, true);
	}
	
	private void evaluateFunction(double tt) throws Exception {
		double f;
		double[] x_t, df = new double[d];
		Model.CompactGradientOutput f_xt, f_wk;
		
		// Loop to ensure NaNs are not encountered during evaluation
		while (true) {
			// Start by setting the minibatch and computing the new position
			model.sampleBatch(b);
			x_t = gp.truePosition(tt);
			// Next, compute the reduced variance function and gradient
			f_xt= stochasticEvaluate(x_t);
			f_wk= stochasticEvaluate(w_k);
			f	= (f_xt.functionValue + fFull) - f_wk.functionValue;
			// Ensure that function value is NOT NaN or Inf
			if (!(Double.isNaN(f) || Double.isInfinite(f))) {
				break;
			}
			// If it is, reduce step size
			tt	= tt/10;
			if (tt < 1E-9) {
				// Unable to rescue optimization; terminate
				throw new Exception("Optimizer Failure: NaN encountered during "
					+ "line search! Try using a smaller step size.");
			}
		}
		// No numerical instability; update GP
		for (int idx=0; idx<d; idx++) {
			df[idx] = (f_xt.gradientVector[idx] + mu_k[idx]) - f_wk.gradientVector[idx];
		}
	    gp.updateGP(tt, x_t, f, df, f_xt.varF, f_xt.varDF);
	}
	
	private GP_Clean.CompactOutput makeOuts(double tt, boolean resetAlpha0) throws Exception {
		double effEta;
		
		// First, ensure that point exists; if not, compute it
    	if (!gp.contains(tt)) {
        	evaluateFunction(tt);
        }
    	// Next, set eta
		if (resetAlpha0) {
			effEta 		= eta;
			effEtaMean	= eta;
		} else {
			// Compute effective eta
			effEta = gp.alpha0*tt;
			// Bound alpha0 to average step size if it goes out of bounds
			if (effEta > 100*effEtaMean || effEta < .01*effEtaMean) { 
				effEta = effEtaMean;
			} else {
				// Update effEtaMean
				effEtaMean = (1-etaEMARate)*effEtaMean*etaGrowthRate + etaEMARate*effEta;
				effEtaMean = Math.min(effEtaMean, eta*maxEtaBound);
				effEta = effEtaMean;
			}
		}

		// Finally, create output, ensuring next step size is smoothed eta
		return gp.makeOuts(tt, effEta);
	}
	
	private double mcsearch(double alphaT) throws Exception {
		boolean useModifiedFunc;
		int nLSEvals;
		double currIntWidth, prevIntWidth, f0, g0, fTTest;
		
		// Check to see if guess alphaT is within bounds
		if (alphaT<stepAlphaMin || alphaT>stepAlphaMax) {
			throw new Exception("GP Line search failure: initial step alpha guess out of bounds!");
		}
		// Initialize variables and see if search direction is a descent 
		// direction; x0 is implicitly set to 0
		f0				= gp.m(0);
		g0				= gp.d1m(0);
		if (g0 > 0) {
			throw new Exception("GP Line search failure: the search direction is not a descent direction!");
		}
		alphaL			= 0;
		fL				= f0;
		gL				= g0;
		alphaU			= 0;
		fU				= f0;
		gU				= g0;
		currIntWidth	= stepAlphaMax-stepAlphaMin;
		prevIntWidth	= 2*currIntWidth;
		useModifiedFunc	= true;
		isBracketed		= false;
		nLSEvals		= 0;
		
		//Create first interval, and ensure that alphaT lies within the acceptable range of alpha values
		alphaT			= Math.min(Math.max(alphaT, stepAlphaMin), stepAlphaMax);
		alphaMin		= alphaL;
		alphaMax		= alphaT + 4*(alphaT-alphaL);
		
		while (nLSEvals<=maxMCSIterations) {
			// Evaluate function
			fT	= gp.m(alphaT);
			gT	= gp.d1m(alphaT);
			nLSEvals++;
			
			// Check for convergence
			fTTest = f0+alphaT*c1*g0;
			if (alphaT==stepAlphaMax && (fT-fTTest<=0 && gT-c1*g0<0)) {
				// \psi(\alpha_t) \leq 0 and \psi^'(\alpha_t) < 0
				throw new ArithmeticException("GP Line search failure: search terminated, step alpha at maximum!");
			} else if (alphaT==stepAlphaMin && (fT-fTTest>0 || gT-c1*g0>=0)) {
				// \psi(\alpha_t) > 0 and \psi^'(\alpha_t) \geq 0
				throw new ArithmeticException("GP Line search failure: search terminated, optimal step alpha is lower than minimum!");
			} if (fT-fTTest<=0 && Math.abs(gT)<=c2*Math.abs(g0)) {
				//converged under strong Wolfe conditions
				return alphaT;
			}
			
			// Check if modified update should be used
			// \psi(\alpha_t) \leq 0 and \phi^'(\alpha_t)>0
			if (useModifiedFunc && fT-fTTest<=0 && gT>=0) {
				useModifiedFunc = false;
			}
			//Generate a new safeguarded alphaT and interval I
			if (useModifiedFunc) {
				alphaT = mcstep(alphaL, fL-alphaL*c1*g0, gL-c1*g0, alphaT, fT-alphaT*c1*g0, gT-c1*g0, 
						alphaU, fU-alphaU*c1*g0, gU-c1*g0);
			} else {
				alphaT = mcstep(alphaL, fL, gL, alphaT, fT, gT, alphaU, fU, gU);
			}

			//Ensure that the interval I decreases sufficiently using a bisection prediction  
			if (isBracketed && (Math.abs(alphaU-alphaL)>=.66*prevIntWidth) ) {
				alphaT = alphaL + (alphaU-alphaL)/2;
			}
			alphaT = Math.min(Math.max(alphaT, stepAlphaMin), stepAlphaMax);
			//Define bounds of new interval
			if (isBracketed) {
				alphaMin = Math.min(alphaL, alphaU);
				alphaMax = Math.max(alphaL, alphaU);
			} else {
				alphaMin = alphaL;
				alphaMax = alphaT+4*(alphaT-alphaL);
			}
			prevIntWidth = currIntWidth;
			currIntWidth = alphaMax-alphaMin;
			//Check to see that the interval is larger than uTol (machine precision)
			if (isBracketed && currIntWidth<=uTol*alphaMax) {
				throw new Exception("GP Line search failure: The relative width of the interval of uncertainty is less than uTol!");
			}
			//If unusual termination occurs, extrapolate from largest point
			if ( (isBracketed && (alphaT<=alphaMin || alphaT>=alphaMax)) || nLSEvals>=maxMCSIterations-1) {
				alphaT = Array.max(gp.T)*2;	
				return alphaT;
			}
		}
		throw new ArithmeticException("GP Line search failure: number of line search function calls exceeded maximum!");
	}
	
	private double mcstep(double alphaLM, double fLM, double gLM, double alphaTM, double fTM, double gTM, 
			double alphaUM, double fUM, double gUM) throws Exception {
		boolean isBound = false;
		double d1, d2, s, p, q, cubic, quadratic, stpf;
		double derivativeSign = Math.signum(gTM)*Math.signum(gLM);
		
		if (gLM*(alphaTM-alphaLM)>=0.0) {
			//Violates Theorem 2.1
			throw new ArithmeticException("GP Line search failure: Interval cannot contain a minimizer!");
		}
		
		//Case 1: Higher function value. Corresponds with case U1
		if(fTM>fLM) {
			isBracketed = true;
			d1 = gLM + gTM + 3*(fLM-fTM)/(alphaTM-alphaLM);
			s = Math.max(Math.abs(d1), Math.max(Math.abs(gLM), Math.abs(gTM)));
			d2 = s*Math.sqrt((d1/s)*(d1/s) - (gLM/s)*(gTM/s));
			if (alphaTM<alphaLM) d2 = -d2;
			p = d1 + d2 - gLM;
			q = gTM + 2*d2 - gLM;
			cubic = alphaLM + p/q*(alphaTM-alphaLM);
			quadratic = alphaLM + (alphaTM-alphaLM)*(gLM/( ( (fLM-fTM)/(alphaTM-alphaLM) ) +gLM) )/2;
			if ( Math.abs(cubic-alphaLM) < Math.abs(quadratic-alphaLM) ) {
				stpf = cubic;
			} else {
				stpf = cubic + (quadratic-cubic)/2;
			}
		} //Case 2: lower function value with derivatives of opposite sign. Corresponds with case U3
		else if(derivativeSign<0) {
			isBracketed = true;
			d1 = gLM + gTM + 3*(fLM-fTM)/(alphaTM-alphaLM);
			s = Math.max(Math.abs(d1), Math.max(Math.abs(gLM), Math.abs(gTM)));
			d2 = s*Math.sqrt((d1/s)*(d1/s) - (gLM/s)*(gTM/s));
			if (alphaTM>alphaLM) d2 = -d2;
			p = d1 + d2 - gTM;
			q = gLM + 2*d2 - gTM;
			cubic = alphaTM + p/q*(alphaLM - alphaTM);
			quadratic = alphaTM + (gTM/(gTM-gLM))*(alphaLM-alphaTM);
			if (Math.abs(cubic-alphaTM) > Math.abs(quadratic-alphaTM)) {
				stpf = cubic;
			} else {
				stpf = quadratic;
			}
		} //Case 3: lower function value, derivatives of the same sign, with decreasing magnitude. Case U2
		else if(Math.abs(gTM)<Math.abs(gLM)) {	
			isBound = true;
			d1 = gLM + gTM + 3*(fLM-fTM)/(alphaTM-alphaLM);
			s = Math.max(Math.abs(d1), Math.max(Math.abs(gLM), Math.abs(gTM)));
			d2 = s*Math.sqrt(Math.max(0, ((d1/s)*(d1/s) - (gLM/s)*(gTM/s)) ));
			if (alphaTM>alphaLM) d2 = -d2;
			p = d1 + d2 - gTM;
			q = gLM + 2*d2 - gTM;
			if (p/q<0 && d2!=0) {		//If cubic tends to infinity in the direction of the step
				cubic = alphaTM + p/q*(alphaLM-alphaTM);
			} 
			else if (alphaTM>alphaLM) {	//If the cubic estimation will be out of bounds, restrict it
				cubic = alphaMax;
			} else {
				cubic = alphaMin;
			}
			quadratic = alphaTM + (gTM/(gTM-gLM))*(alphaLM-alphaTM);
			if(isBracketed) {
				if(Math.abs(alphaTM-cubic) < Math.abs(alphaTM-quadratic)) {
					stpf = cubic;
				} else {
					stpf = quadratic;
				}
			} //Since the interval has not been bracketed, along with case U2
			//the following conditions satisfy safeguarding conditions 2.2+2.3 and 2.4+2.5
			else {
				if (alphaTM>alphaLM) {
					stpf = alphaMax;
				} else {
					stpf = alphaMin;
				}
			}
		} //Case 4: Lower function value, derivatives have same sign without decreasing magnitude. Case U2
		else {
			if (isBracketed) {
				d1 = gUM + gTM + 3*(fTM-fUM)/(alphaUM-alphaTM);
				s = Math.max(Math.abs(d1), Math.max(Math.abs(gUM), Math.abs(gTM)));
				d2 = 2*Math.sqrt((d1/s)*(d1/s) - (gUM/s)*(gTM/s));
				if (alphaTM>alphaUM) d2 = -d2;
				p = d1 + d2 - gTM;
				q = gUM + 2*d2 - gTM;
				cubic = alphaTM + p/q*(alphaUM-alphaTM);
				stpf = cubic;
			} //Since the interval has not been bracketed, along with case U2
			//the following conditions satisfy safeguarding conditions 2.2+2.3 and 2.4+2.5
			else if (alphaTM>alphaLM) {
				stpf = alphaMax;
			} else {
				stpf = alphaMin;
			}
		}
		
		//Update interval of uncertainty (Updating algorithm)
		if (fTM>fLM) {
			alphaU = alphaTM;
			fU = fTM;
			gU = gTM;
		} else {
			if (derivativeSign<0) {
				alphaU = alphaLM;
				fU = fLM;
				gU = gLM;
			}
			alphaL = alphaTM;
			fL = fTM;
			gL = gTM;
		}
		
		//Compute new safeguarded step
		stpf = Math.min(Math.max(stpf, alphaMin), alphaMax);
		if (isBracketed && isBound) {	//Modified version of bounds in case3?
			if (alphaUM>alphaLM) {
				stpf = Math.min(alphaLM + 0.66*(alphaUM-alphaLM), stpf);
			} else {
				stpf = Math.max(alphaLM + 0.66*(alphaUM-alphaLM), stpf);
			}
		}
		return stpf;
	}
}