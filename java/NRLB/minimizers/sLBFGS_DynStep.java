package minimizers;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

import base.*;

public class sLBFGS_DynStep extends Minimizer{
	private boolean isVerbose, averageSteps = false;
	private int d, N, b, bH, M, m, L, currDepth, maxEpoch;
	private double eta, delta, fdHVPStepSize = 5E-2, dirNormBound = 10;
	private double fFull, effEtaMean;
	private double[] rho, w_k, mu_k;
	private double[][] s, y;
	private Fit fitOutput;
	private Model.CompactGradientOutput fOut;
	// Gaussian process used for evaluations
	GP gp;
	
	// Test parameters
    boolean testNegGradStepSize		= false;
    boolean testExceptionStepSize	= false;
    boolean testLSFailureStepSize	= false;
    double negGradStepSize 			= .01;
    double exceptionStepSize		= .01;
    double lsFailureStepSize		= .01;
    double etaGrowthRate			= 1.005;
    double etaDecayRate				= .1;
    
    PrintStream gpOutputFile, original;
    
	
	//LBFGS object constructor; load basic minimization parameters. To be used for all subsequent minimizations using this object.
	public sLBFGS_DynStep(Model model, int gradientBatch, int hessianBatch, int memorySize, 
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
		boolean firstEval = true;
		int r = 0;									// Number of currently computed Hessian correction pairs
		currDepth = 0;
		double tStart, dirNorm, divisor, f_t, effEta;
		// Full gradient, variance reduced gradient, and components of variance reduced gradient
		double[] v_t = new double[d];
		// Effective gradient
		double[] dir = new double[d];
		// Positions in current iterations
		double[] w_k_prev = new double[d], x_t = new double[d];
		// Average of path travelled in the current and previous inverse hessian updates
		double[] u_r = new double[d], u_r_prev = new double[d];
		// Components of two-loop update
		double[] s_r = new double[d], y_r= new double[d];
		rho = new double[M];
		s = new double[M][d];
		y = new double[M][d];
		// Stores the history of all previous steps in the current epoch
		ArrayList<double[]> x_t_hist = new ArrayList<double[]>();
		GP.CompactOutput plsOut = null;
		Model.CompactGradientOutput f_xt, f_wk;
		
		// Init 
		effEta		= eta;
		effEtaMean	= eta;
		w_k			= new double[d];
		mu_k		= new double[d];
		
		//TODO:
		double gNormMean = 0, dirNormMean = 0, desiredNorm = 1;
		double gDecayRate = 1;			// Decay rate for function gradient
		double egDecayRate= .1;			// egNorm decay rate
		PrintStream outputFile	= new PrintStream(new FileOutputStream("/Users/chaitanya/Documents/GitWorkspaces/slbfgs/java/output/RawPath_2.txt"));
		original	= System.out;
		gpOutputFile= new PrintStream(new FileOutputStream("/Users/chaitanya/Documents/GitWorkspaces/slbfgs/java/output/GP_Path.txt"));
		System.setOut(gpOutputFile);
		System.out.println("c1: NA");
		System.out.println("c2: NA");
		System.setOut(original);
		
		// Deal with a potential seed and initialize
		if (seed!=null) {
			try {
				w_k = Array.clone(seed);
			} catch (Exception e) {
				outputFile.close();
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
					//TODO:
					gNormMean = Array.norm(fOut.gradientVector);
					dirNormMean= gNormMean;
					
					System.out.println("Starting Function Value: "+fOut.functionValue);
					System.out.println("    Epochs   Data Loops           Likelihood       Distance Moved"+
							"        Gradient Norm");
				} else {
					if (xStar!=null) {
						printStep(k, model.evaluatedDataPoints/N, fOut.functionValue, Array.dist(w_k, w_k_prev), 
								Array.norm(fOut.gradientVector), Array.dist(w_k, xStar));	
					} else {
						printStep(k, model.evaluatedDataPoints/N, fOut.functionValue, Array.dist(w_k, w_k_prev), 
								Array.norm(fOut.gradientVector), effEtaMean);	
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
			x_t  = Array.clone(w_k);							// Set x_t to current value of w_k
			w_k_prev = Array.clone(w_k);					
			
			//TODO:
			System.setOut(outputFile);
			System.out.print(k+"\t"+m+"\t>\t\t");
			Array.print(Array.cat(Array.cat(w_k, mu_k), Array.cat(new double[4*d], fOut.functionValue)));
			System.setOut(original);
			
			// Perform m stochastic iterations before a full gradient computation takes place
			for (int t=1; t<=m; t++) {
				// Need to run a different update step if it is the first eval
				if (firstEval) { 
					model.sampleBatch(b);
					f_xt = stochasticEvaluate(x_t);
					f_wk = stochasticEvaluate(w_k);
					f_t = (f_xt.functionValue + fFull) - f_wk.functionValue;
					v_t = Array.subtract(Array.add(f_xt.gradientVector, mu_k), f_wk.gradientVector);
					u_r = Array.add(u_r, x_t);
					x_t_hist.add(Array.clone(x_t));
					dir = v_t;
					dirNorm = Array.norm(dir);
					divisor = 1;
					while (dirNorm/divisor > dirNormBound) {
						divisor *= 10;
					}
					plsOut = probLineSearch(x_t, f_t, v_t, f_xt.varF, f_xt.varDF, effEta, Array.scalarMultiply(dir, -1.0/divisor));
					firstEval = false;
					continue;
				}
				
                //TODO:
				String hessUpdate="";
                
				// Check to see if L iterations have passed to trigger Hessian
				// update; need to do this NOW before NEXT pls call
				if (t % L == 0) {
					//TODO:
//					System.out.println("Hessian update");
					
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
					
					//TODO
					hessUpdate = Double.toString(Array.dist(u_r, u_r_prev));
					
					// Resetting u_r for next evaluation
			        u_r_prev = Array.clone(u_r);
			        u_r = new double[d];
				}
				
				// Copy values from pls (correspond to NEXT position)
				x_t 	= plsOut.position;
				f_t 	= plsOut.functionValue;
				v_t 	= plsOut.gradientVector;
				effEta	= plsOut.effEta;
				
				//TODO
				gNormMean = gDecayRate*Array.norm(v_t)+(1-gDecayRate)*gNormMean;
				v_t = Array.scalarMultiply(v_t, gNormMean/Array.norm(v_t));
				
				// Update u_r with current position
				u_r = Array.add(u_r, x_t);			
				// Compute next iteration step; condition the gradient so as not to produce rapid oscilations (extreme function values)
				x_t_hist.add(Array.clone(x_t));		// Need to store the history of iteration steps				
				if (r < 1) {						// Until a single hessian correction has taken place, H_0 = I
					dir = v_t;
				} else {							// Compute the two-loop recursion product
					dir = twoLoopRecursion(v_t);
				}			
				// Bound the effective gradient update step
				dirNorm = Array.norm(dir);
				
				//TODO:
				// See if egNorm spikes
				if (dirNorm > 2.718282*dirNormMean && r>1) { 
					desiredNorm = dirNormMean;
				} else {
					desiredNorm = dirNorm;
				}
				// Bound gradient update
				desiredNorm = Math.min(desiredNorm, Math.sqrt(d)*dirNormBound);
				// Compute exponential moving average
				dirNormMean = egDecayRate*desiredNorm+(1-egDecayRate)*dirNormMean;
				divisor = dirNorm/dirNormMean;
				// Run line search
				plsOut = probLineSearch(x_t, f_t, v_t, plsOut.varF, plsOut.varDF, effEta, Array.scalarMultiply(dir, -1.0/divisor));
				
				//TODO: Write to file
				System.setOut(outputFile);
				System.out.print(k+"\t"+t+"\t\t"+hessUpdate+"\t");
				Array.print(Array.cat(Array.cat(Array.cat(x_t, v_t), new double[2*d]), Array.cat(mu_k, Array.cat(dir, dirNorm/divisor))));
				System.setOut(original);
			}
			// Choose either the last position vector x_t or a random one from the previous epoch as the starting point for the next
			if (averageSteps) {
				w_k = Array.clone(x_t_hist.get(0));
				for (int currPos=1; currPos<x_t_hist.size(); currPos++) {
					w_k = Array.add(w_k, x_t_hist.get(currPos));
				}
				w_k = Array.scalarMultiply(w_k, 1.0/((double) x_t_hist.size()));
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
	
	private GP.CompactOutput probLineSearch(double[] x0, double f0, double[] df0, 
			double var_f0, double[] var_df0, double alpha0, double[] dir) throws Exception{
		double tt, beta, effEta, zScore;
		
		// Create new GP and update at start
		beta		= Array.norm(df0);		// Compute scaling factor beta
		gp 			= new GP(0, 0, f0, beta, alpha0, dir, gpOutputFile, original, model);
		tt			= 0;
	    gp.updateGP(tt, x0, f0, df0, var_f0, var_df0);
		
	    // Take unit alpha0 step TODO: include NaN catching here
	    evaluateFunction(1.0);
	    
	    // Compute effEta using alpha0 and z-score
//	    zScore = (gp.m(0)-gp.m(1.0))/Math.sqrt(Math.max(gp.V(0), gp.V(1.0)));
	    zScore = gp.d1m(0)/Math.sqrt(gp.dVd(0));
	    
	    if (zScore>0) {
	    	effEta = effEtaMean*Math.pow(2, -zScore);
	    } else {
	    	effEta = effEtaMean*(1+(etaGrowthRate-1)*Math.exp(zScore));
	    }
	    
	    
//	    effEta = alpha0*Math.pow(etaGrowthRate, -zScore);
//	    effEta = effEtaMean*Math.pow(etaGrowthRate, -zScore);
//	    if (zScore >= 0) {
//	    	effEta = effEtaMean*1.01;
//	    } else {
//	    	effEta = effEtaMean*Math.pow(1.5, zScore);
//	    }
	    
	    // Bound alpha0 TODO: make this bound less stringent
	    effEta = Math.min(effEta, eta*10);
	    effEtaMean = (1-etaDecayRate)*effEtaMean + etaDecayRate*effEta;
//	    effEtaMean = (1-etaDecayRate)*effEtaMean + etaDecayRate*Math.min(alpha0, eta*10);
	    
	    // Return output
	    //System.out.println("newEffEta: "+effEtaMean);
	    return gp.makeOuts(1.0, effEtaMean);
	}
	
	private void evaluateFunction(double tt) throws Exception {
		double f;
		double[] x_t, df = new double[d];
		Model.CompactGradientOutput f_xt, f_wk;
		
		// Start by setting the minibatch and computing the new position
		model.sampleBatch(b);
		x_t		= gp.truePosition(tt);
		// Next, compute the reduced variance function and gradient
		f_xt	= stochasticEvaluate(x_t);
		f_wk	= stochasticEvaluate(w_k);
		f		= (f_xt.functionValue + fFull) - f_wk.functionValue;
		for (int idx=0; idx<d; idx++) {
			df[idx] = (f_xt.gradientVector[idx] + mu_k[idx]) - f_wk.gradientVector[idx];
		}
	    
	    // Test for NaN or Inf and reset
		if (Double.isNaN(f) || Double.isInfinite(f)) {
			
			System.out.println("TT: "+tt);
			Array.print(gp.T);
			System.out.println(f_xt);
			System.out.println(f_wk);
			gp.dump2file();
			
			throw new Exception("Optimizer Failure: NaN encountered during "
					+ "line search! Try using a smaller step size.");
		}
	    
	    // No numerical instability; update GP
	    gp.updateGP(tt, x_t, f, df, f_xt.varF, f_xt.varDF);
	}
}