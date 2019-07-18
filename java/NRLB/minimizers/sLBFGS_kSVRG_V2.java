package minimizers;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;

import base.*;

public class sLBFGS_kSVRG_V2 extends Minimizer{
	private boolean isVerbose;
	private int k, d, N, bs, bH, currDepth, maxEpoch, memDepth, L, eIters, maxBindings;
	private double eta, delta, epsilon, gradientNormBound = 100;
	private double fdHVPStepSize = 5E-2;
	private int[] PhiCount;
	private double[] rho, PhiProb;
	private double[][] s, y;
	private ArrayList<Integer>[] PhiIdxes;
	private Fit fitOutput;
	private MersenneTwisterFast mtfast = new MersenneTwisterFast();
	
	// kSVRG object constructor; load basic minimization parameters. To be used 
	// for all subsequent minimizations using this object.
	public sLBFGS_kSVRG_V2(Model model, int k, int gradientBatch, int hessianBatch,
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
		int min_TIdx, max_TIdx, currPhiSize, r = 0, huIters = 0, group, currGroup;
		currDepth = 0;
		double egNorm, divisor, avgF;
		int[] theta_m_binding, theta_m1_binding, ind, batchIdx, Phi;
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
		ArrayList<Integer> PhiArray, groupIdxes = new ArrayList<Integer>(bs);
		ArrayList<Integer>[] subBindings;
		Model.CompactGradientOutput fOut;
				
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
		PhiCount= new int[1];
		PhiCount[0]= N;
		PhiProb = new double[1];
		PhiProb[0] = 1.0;
		PhiIdxes = new ArrayList[1];
		PhiIdxes[0] = new ArrayList<Integer>();
		for (int idx=0; idx<N; idx++) {
			PhiIdxes[theta_m_binding[idx]].add(idx);
		}
		// aBar_m must be set to the gradient on the whole dataset
		aBar_m 	= gradientEval(x0, true).gradientVector;
		
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
			
			//TODO:
//			System.out.println("Start of new epoch m = "+m);
//			System.out.print("ind:\t");
//			Array.print(ind);
			
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
				
				//TODO:
//				System.out.print("New Bindings:\t");
//				Array.print(theta_m_binding);
//				System.out.println("New Splitting:");
//				for (int currM=0; currM<maxBindings; currM++) {
//					System.out.print("M: "+currM+"\t");
//					for (int blah=0; blah<PhiIdxes[currM].size(); blah++) {
//						System.out.print(PhiIdxes[currM].get(blah)+"\t");
//					}
//					System.out.print("\n");
//				}
				
				// Loop over epoch iterates
				for (int t=0; t<eIters; t++) {
					// Update average stores
					for (int currIdx=0; currIdx<d; currIdx++) {
						avg_x[currIdx] += x_t[currIdx];
						avg_g[currIdx] += g_t[currIdx];
					}
					// Create random indices 
					group = randGroup();	// Select group
					randIdxes(group, groupIdxes);
					//TODO:
//					System.out.println("For step t = "+t+"; Group:\t"+group);
//					for (int cidx=0; cidx<bs; cidx++) {
//						System.out.print(groupIdxes.get(cidx)+"\t");
//					}
//					System.out.print("\n");
					
					// Compute a_it and f_it
					fOut = stochasticEvaluate(x_t, groupIdxes);
					avgF +=fOut.functionValue;
					f_it = fOut.gradientVector;
					a_it = stochasticEvaluate(theta_m[group], groupIdxes).gradientVector;
					
					// Compute SVRG gradient
					for (int idx=0; idx<d; idx++) {
						g_t[idx] = f_it[idx] - a_it[idx] + aBar_m[idx];
					}
					// Update u_r with current position
					u_r = Array.add(u_r, x_t);
					
					// Compute next iteration step; condition the gradient so as not to produce rapid oscilations (extreme function values)
					if (r < 1) {						// Until a single hessian correction has taken place, H_0 = I
						effGrad = g_t;
					} else {							// Compute the two-loop recursion product
						effGrad = twoLoopRecursion(g_t);
					}			
					// Bound the effective gradient update step
					egNorm = Array.norm(effGrad);
					divisor = 1;
					while (egNorm/divisor > gradientNormBound) {
						divisor *= 10;
					}
					// Use typical update step
					x_t = Array.addScalarMultiply(x_t, -eta/divisor, effGrad);
					huIters++;
					
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
						// Store latest values of s_r and y_r
						add(s_r, y_r);
						
						// Resetting u_r for next evaluation
				        u_r_prev = Array.clone(u_r);
				        u_r = new double[d];
					}
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
				// Reset phiCounts; this is inefficient
				PhiCount= new int[maxBindings];
				PhiProb = new double[maxBindings];
				PhiIdxes = new ArrayList[maxBindings];
				for (int currM=0; currM<maxBindings; currM++) {
					PhiIdxes[currM] = new ArrayList<Integer>();					
				}
				for (int idx=0; idx<N; idx++) {
					currGroup = theta_m_binding[idx];
					PhiCount[currGroup]++;
					PhiIdxes[currGroup].add(idx);
				}
				PhiProb[0] = PhiCount[0]/((double) N);
				for (int currM=1; currM<maxBindings; currM++) {
					PhiProb[currM] = PhiProb[currM-1] + PhiCount[currM]/((double) N);					
				}
				PhiProb[maxBindings-1] = 1.0;
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
	
	// Generate random group
	private int randGroup() {
		// Start by generating a random double
		double randD = mtfast.nextDouble();

		// Return group
		for (int currGroup=maxBindings-2; currGroup>=0; currGroup--) {
			if (PhiProb[currGroup]<=randD) {
				return (currGroup+1);
			}
		}
		return 0;
	}
	
	// Generate random ordering for a group given batchsize, pass by reference
	private void randIdxes(int group, ArrayList<Integer> input) {
		HashSet<Integer> h = new HashSet<Integer>();	
		
		// Create list of candidates
		for (int i=0; i<k; i++) {
			h.add(mtfast.nextInt(PhiCount[group]));
		}
		while(h.size()<k) {
			h.add(mtfast.nextInt(PhiCount[group]));
		}
		
		// Copy into input
		input.clear();
		Iterator<Integer> i = h.iterator(); 
        while (i.hasNext()) {
        	input.add(PhiIdxes[group].get(i.next()));
        }
	}
}