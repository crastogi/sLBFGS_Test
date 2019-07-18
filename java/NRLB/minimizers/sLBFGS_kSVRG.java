package minimizers;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import base.*;

public class sLBFGS_kSVRG extends Minimizer{
	private boolean isVerbose;
	private int k, d, N, bs, bH, currDepth, maxEpoch, memDepth, L, eIters;
	private double eta, delta, epsilon, gradientNormBound = 2;
	private double fdHVPStepSize = 5E-3;
	private double[] rho;
	private double[][] s, y;
	private Fit fitOutput;
	private MersenneTwisterFast mtfast = new MersenneTwisterFast();
	
	// kSVRG object constructor; load basic minimization parameters. To be used 
	// for all subsequent minimizations using this object.
	public sLBFGS_kSVRG(Model model, int k, int gradientBatch, int hessianBatch,
			int maxEpoch, int hessianPeriod, int memorySize, double stepsize, 
			double epsilon, double delta, boolean isVerbose) {
		this.model		= model;
		d				= model.getNDimensions();	// Dimensionality
		N				= model.nCount;				// Number of data points
		this.k			= k;
		if (gradientBatch<2) {
			throw new IllegalArgumentException("Gradient batchsize must be >1!");
		}
		bs 				= gradientBatch;
		if (hessianBatch<1) {
			throw new IllegalArgumentException("Hessian batchsize must be >1!");
		}
		bH				= hessianBatch;			// Batch size for stochastic Hessian updates
		// The maximum number of epochs before termination
		if (maxEpoch<1) {
			throw new IllegalArgumentException("Max epochs must be >1!");
		}
		this.maxEpoch	= maxEpoch;
		// Number of iterations before the inverse hessian is updated
		if (hessianPeriod<2) {
			throw new IllegalArgumentException("Hessian period must be >1!");
		}
		L				= hessianPeriod;
		if (memorySize<=0)	{
			throw new IllegalArgumentException("Memory size must be positive!");
		}
		memDepth		= memorySize;			// Maximum memory depth

		// Set the number of iterations per epoch
		eIters = (int) Math.floor(((double) N)/((double) (k*bs)));
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
		// An optional inverse hessian regularization parameter
		if (delta<0) {
			throw new IllegalArgumentException("Delta (inverse hessian "
					+ "regularization parameter) cannot be negative!");
		}
		this.delta		= delta;				
		this.isVerbose	= isVerbose;		// Verbosity
	}
	
	public Fit doMinimize(double[] seed, String trajectoryFile) throws Exception {
		int maxBindings, min_TIdx, max_TIdx, currPhiSize, r = 0, huIters = 0;
		currDepth = 0;
		double egNorm, divisor, avgF;
		int[] theta_m_binding, theta_m1_binding, ind, batchIdx, Phi, bsInd=null;
		double[] x0, x_t, avg_x, avg_g, x_t_prev, g_t, f_it, a_it, aBar_m, xBar;
		double[] effGrad;
		// Average of path traveled in the current and previous inverse Hessian updates
		double[] u_r = new double[d], u_r_prev = new double[d];
		// Components of two-loop update
		double[] s_r = new double[d], y_r= new double[d];
		rho = new double[memDepth];
		s = new double[memDepth][d];
		y = new double[memDepth][d];
		int[][] phiBlocks = idxAlloc();
		double[][] theta_m, theta_m1;
		ArrayList<Integer> PhiArray;
		ArrayList<Integer>[] subBindings;
		Model.CompactGradientOutput fOut;

		//TODO:
		double gNormMean = 0, egNormMean = 0, desiredNorm = 1;
		double gDecayRate = .1;			// Decay rate for function gradient
		double egDecayRate= .1;			// egNorm decay rate
		PrintStream outputFile	= new PrintStream(new FileOutputStream("/Users/chaitanya/Documents/GitWorkspaces/slbfgs/java/output/RawPath_2.txt"));
		PrintStream original	= System.out;

		/*Check to see if seed is of correct length and symmetrize seed, if need
		 * be. If no seed is provided, use a random start point to break 
		 * symmetry. */
		if (seed!=null) {
			x0 = Array.clone(seed);
		} else {
			x0 = Array.scalarMultiply(Array.randomDouble(d), .05);
		}
		fitOutput = model.generateFit(seed);
		
		// Initialize parameters
		x_t		= Array.clone(x0);
		x_t_prev= Array.clone(x0);
		g_t 	= new double[d];
		theta_m = new double[1][d];
		maxBindings = 1;
		theta_m[0] = Array.clone(x0);
		theta_m_binding = new int[N];
		theta_m1_binding= new int[N];
		// aBar_m must be set to the gradient on the whole dataset
		aBar_m 	= gradientEval(x0, true).gradientVector;
		
		//TODO:
		gNormMean = Array.norm(aBar_m);
		egNormMean= gNormMean;
		
		// Print reporting line
		if (isVerbose) {
			if (xStar!=null) {
				System.out.println("    Epochs   Data Loops        Avg. Func "
						+ "Val       Distance Moved        Gradient Norm"
						+ "       Distance to x*");
			} else {
				System.out.println("    Epochs   Data Loops        Avg. Func"
						+ " Val       Distance Moved        Gradient Norm");
			}
		}
		
		// Outer loop over epochs
		for (int m=0; m<maxEpoch; m++) {
			// Create a random permutation of the points
			ind = fisherYates();
			bsInd = fisherYates();
			
			//TODO:
//			System.out.println("Start of new epoch m = "+m);
//			System.out.print("ind:\t");
//			Array.print(ind);
//			System.out.print("bsInd:\t");
//			Array.print(bsInd);
			
			// Loop over sub-ks
			for (int j=0; j<k; j++) {
				
				//TODO:
//				System.out.println("Start of new sub-k "+j+"; # of current bindings = "+maxBindings);
				
				// Zero-out averages
				avgF  = 0;
				avg_x = new double[d];
				avg_g = new double[d];
				subBindings = new ArrayList[maxBindings];
				for (int currM=0; currM<maxBindings; currM++) {
					subBindings[currM] = new ArrayList<Integer>();					
				}
				// Loop over epoch iterates
				for (int t=0; t<eIters; t++) {
//					// Update average stores
//					for (int currIdx=0; currIdx<d; currIdx++) {
//						avg_x[currIdx] += x_t[currIdx];
//						avg_g[currIdx] += g_t[currIdx];
//					}
					// Update batch idx directly
					//model.currBatchIdx = Arrays.copyOfRange(bsInd, 
					//		j*eIters*bs+t*bs, j*eIters*bs+(t+1)*bs);
					//model.currBatchSize = bs;
					model.sampleBatch(bs);
					batchIdx = model.currBatchIdx;
					
					//TODO:
//					System.out.println("For step t = "+t);
//					System.out.print("cbIdx:\t");
//					Array.print(batchIdx);
					
					fOut = stochasticEvaluate(x_t);
					avgF +=fOut.functionValue;
					f_it = fOut.gradientVector;
					// Create a new sub binding array
					for (int currM=0; currM<maxBindings; currM++) {
						subBindings[currM].clear();
					}
					// Re-assign values in subBinding array
					for (int idx=0; idx<bs; idx++) {
						subBindings[theta_m_binding[batchIdx[idx]]].add(batchIdx[idx]);
					}
					// Now run stochastic evaluation on every theta_m batch
					a_it = new double[d];
					for (int currM=0; currM<maxBindings; currM++) {
						if (subBindings[currM].size()==0) {
							continue;
						}
						a_it = Array.add(a_it, Array.scalarMultiply(stochasticEvaluate(theta_m[currM], subBindings[currM]).gradientVector, subBindings[currM].size()));
					}
					// Compute SVRG gradient
					for (int idx=0; idx<d; idx++) {
						g_t[idx] = f_it[idx] - a_it[idx]/bs + aBar_m[idx];
					}
					// Update u_r with current position
					u_r = Array.add(u_r, x_t);
					
					//TODO
					gNormMean = gDecayRate*Array.norm(g_t)+(1-gDecayRate)*gNormMean;
					g_t = Array.scalarMultiply(g_t, gNormMean/Array.norm(g_t));
					
					// Compute next iteration step; condition the gradient so as not to produce rapid oscilations (extreme function values)
					if (r < 1) {						// Until a single hessian correction has taken place, H_0 = I
						effGrad = g_t;
					} else {							// Compute the two-loop recursion product
						effGrad = twoLoopRecursion(g_t);
					}			
					egNorm = Array.norm(effGrad);
					
					if (Double.isNaN(egNorm) || Double.isInfinite(egNorm)) {
						throw new Exception("sLBFGS+kSVRG Failure: NaN encountered! Try reducing the step size...");
					}
					
					//TODO:
					// See if egNorm spikes
					if (egNorm > 2.718282*egNormMean && r>1) {
						desiredNorm = egNormMean;
					} else {
						desiredNorm = egNorm;
					}
					// Bound gradient update
					desiredNorm = Math.min(desiredNorm, d*gradientNormBound);
					// Compute exponential moving average
					egNormMean = egDecayRate*desiredNorm+(1-egDecayRate)*egNormMean;
					divisor = egNorm/egNormMean;
					
					// Use typical update step
					x_t = Array.addScalarMultiply(x_t, -eta/divisor, effGrad);
					huIters++;
					
					//TODO: test new placement here
					// Update average stores
					for (int currIdx=0; currIdx<d; currIdx++) {
						avg_x[currIdx] += x_t[currIdx];
						avg_g[currIdx] += g_t[currIdx];
					}
					
					//TODO:
					String hessUpdate="";
					
					// Check to see if L iterations have passed (triggers hessian update)
					if (huIters % L == 0) {
						// Increment the number of hessian correction pairs
						r++;
						huIters = 0;
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
						
						//TODO
						y_r = Array.addScalarMultiply(y_r, 0.0, s_r);
						
						// Store latest values of s_r and y_r
						add(s_r, y_r);
						
						//TODO
						hessUpdate = Double.toString(Array.dist(u_r, u_r_prev));
						
						// Resetting u_r for next evaluation
				        u_r_prev = Array.clone(u_r);
				        u_r = new double[d];
					}
					
					//TODO: Write to file
					System.setOut(outputFile);
					System.out.print(m+"\t"+j+"\t"+t+"\t"+hessUpdate+"\t");
					Array.print(Array.cat(Array.cat(Array.cat(x_t, g_t),Array.cat(f_it, a_it)), Array.cat(aBar_m, Array.cat(effGrad, egNorm/divisor))));
					System.setOut(original);
				}

				// Phi Array
				Phi = Arrays.copyOfRange(ind, phiBlocks[j][0], phiBlocks[j][1]);
				//Phi = Arrays.copyOfRange(ind, j*eIters*bs, (j+1)*eIters*bs);
				// Compute xBar, letting S_l = l (or eIters in this case)
				xBar = Array.scalarMultiply(avg_x, 1.0/((double) eIters));
				// Need to set theta_m1_binding, compute min/max_Tidx, etc.
				// Start by deep copying theta_m_binding and computing the min
				PhiArray = new ArrayList<Integer>(Phi.length);
				currPhiSize = Phi.length;
				theta_m1_binding = Array.clone(theta_m_binding);
				// Now update the theta_m value of the recently evaluated points
				max_TIdx = maxBindings;
				for (int idx : Phi) {
					theta_m1_binding[idx] = max_TIdx;
					PhiArray.add(idx);
				}
				min_TIdx = Array.min(theta_m1_binding);
				// Update full gradient estimate: start by computing on new location
				f_it = stochasticEvaluate(xBar, PhiArray).gradientVector;
				a_it = new double[d];
				// Compute update components from previous theta_m's on Phi. Start
				// by creating a new sub binding array
				for (int currM=0; currM<maxBindings; currM++) {
					subBindings[currM].clear();
				}
				// Re-assign values in subBinding array
				for (int idx=0; idx<Phi.length; idx++) {
					subBindings[theta_m_binding[Phi[idx]]].add(Phi[idx]);
				}
				
				// Now run stochastic evaluation on every theta_m batch
				for (int currM=0; currM<maxBindings; currM++) {
					if (subBindings[currM].size()==0) {
						continue;
					}
					a_it = Array.add(a_it, Array.scalarMultiply(stochasticEvaluate(theta_m[currM], subBindings[currM]).gradientVector, subBindings[currM].size()));
				}
				// Finally, recompute aBar_m
				for (int idx=0; idx<d; idx++) {
					aBar_m[idx] += (f_it[idx]*currPhiSize-a_it[idx])/((double) N);
				}
				
				if (isVerbose) {
					fOut = gradientEval(xBar, true);
					model.evaluatedDataPoints -= N;
					if (xStar!=null) {
						printStep(m+1, model.evaluatedDataPoints/N, avgF/eIters,
								Array.norm(Array.subtract(xBar, x_t_prev)), 
								Array.norm(aBar_m), Array.norm(Array.subtract(xBar, xStar)), fOut.functionValue, Array.norm(fOut.gradientVector));
					} else {
						printStep(m+1, model.evaluatedDataPoints/N, avgF/eIters,
								Array.norm(Array.subtract(xBar, x_t_prev)), 
								Array.norm(aBar_m), fOut.functionValue, Array.norm(fOut.gradientVector));
					}
					
					//TODO:
					System.setOut(outputFile);
					System.out.print(m+"\t"+j+"\t>\t\t");
					Array.print(Array.cat(Array.cat(xBar, aBar_m), Array.cat(new double[4*d], fOut.functionValue)));
					System.setOut(original);
				}
				x_t_prev = Array.clone(xBar);
				// Check for convergence
				if (Array.norm(aBar_m)/Math.max(1, Array.norm(xBar)) < epsilon) {
					fitOutput.recordFit(m+1, -1, -1, -1, model);
					System.out.println("Convergence criteria met.");
					model.hessianEval();
					fitOutput.storeHessian(model.getHessian());
					return fitOutput;
				}
				// Not converged; reallocate the theta_m matrix
				theta_m1 = new double[max_TIdx-min_TIdx+1][d];
				for (int currTIdx=min_TIdx; currTIdx<max_TIdx; currTIdx++) {
					theta_m1[currTIdx-min_TIdx] = theta_m[currTIdx];
				}
				theta_m1[theta_m1.length-1] = xBar;
				theta_m = theta_m1;
				for (int idx=0; idx<N; idx++) {
					theta_m_binding[idx] = theta_m1_binding[idx]-min_TIdx;
				}
				maxBindings = theta_m1.length;
			}
		}
		throw new Exception("kSVRG Failure: maximum epochs exceeded without convergence!");
	}
	
	// Randomly permute array
	private int[] fisherYates() {
		int idx, swap;
		int[] output = new int[N];

		// Init output array
		for (int i=0; i<N; i++) {
			output[i] = i;
		}
		
		// Perform Fisher-Yates shuffle
		for (int i=N-1; i>0; i--) {
			idx = mtfast.nextInt(i);
			swap = output[i];
			output[i] = output[idx];
			output[idx] = swap;
		}
		
		return output;
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
		currDepth = Math.min(currDepth+1, memDepth);
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
	
	private int[][] idxAlloc() {
		int blockSize, addPoints, currPoint=-1;
		int[][] output = new int[k][2];
		
		// First, determine the size of blocks and additional points
		blockSize = (int) Math.floor(((double) N)/((double) k));
		addPoints = N - blockSize*k;
		
		// Now enter the indices in the output array
		for (int currK=0; currK<k; currK++) {
			// Start of block
			currPoint++;
			output[currK][0] = currPoint;
			// Compute end point of block
			currPoint += blockSize-1;
			// Adjust if we have extra points that need to be accounted for
			if (addPoints>0) {
				currPoint++;
				addPoints--;
			}
			output[currK][1] = currPoint+1;		// Need +1 offset for array copy
		}
		
		return output;
	}
}