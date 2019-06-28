package minimizers;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import base.*;

public class kSVRG extends Minimizer{
	private boolean isVerbose, useADAM, usek2;
	private int k, d, N, bs, M, eIters;
	private double eta, epsilon, gradientNormBound = 100;
	private Fit fitOutput;
	private MersenneTwisterFast mtfast = new MersenneTwisterFast();
	
	// NOTE: ADAM parameters below
	private double b1 = .9, b2 = .999, e = 1E-8;
	
	// kSVRG object constructor; load basic minimization parameters. To be used 
	// for all subsequent minimizations using this object.
	public kSVRG(Model model, int k, boolean usek2, boolean isStochastic, 
			int gradientBatch, int maxEpoch, double stepsize, double epsilon, 
			boolean useADAM, boolean isVerbose) {
		this.model		= model;
		d				= model.getNDimensions();	// Dimensionality
		N				= model.nCount;				// Number of data points
		this.k			= k;
		this.usek2		= usek2;
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
		this.useADAM	= useADAM;			// Should ADAM be used
		this.isVerbose	= isVerbose;		// Verbosity
	}
	
	public Fit doMinimize(double[] seed, String trajectoryFile) throws Exception {
		if (usek2) {
			return run_k2SVRG(seed, trajectoryFile);
		} else {
			return run_kSVRG(seed, trajectoryFile);
		}
	}
	
	private Fit run_kSVRG(double[] seed, String trajectoryFile) throws Exception {
		int totIters = 0, maxBindings, min_TIdx, max_TIdx, currPhiSize;
		double gNorm, divisor, v_tSum, a_t;
		int[] theta_m_binding, batchIdx, theta_m1_binding;
		double[] x0, x_t, avg_x, avg_g, x_t_prev, g_t, m_t, v_t, f_it, a_it;
		double[] aBar_m, xBar;
		double[][] theta_m, theta_m1;
		HashSet<Integer> Phi;
		ArrayList<Integer> PhiArray;
		ArrayList<Integer>[] subBindings;
		
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
		m_t 	= new double[d];
		v_t 	= new double[d];
		theta_m = new double[1][d];
		maxBindings = 1;
		theta_m[0] = Array.clone(x0);
		theta_m_binding = new int[N];
		theta_m1_binding= new int[N];
		// aBar_m must be set to the gradient on the whole dataset
		aBar_m 	= gradientEval(x0, true).gradientVector;
		
		// Print reporting line
		if (isVerbose) {
			if (xStar!=null) {
				System.out.println("    Epochs   Data Loops       Distance "+
						"Moved        Gradient Norm       Distance to x*");
			} else {
				System.out.println("    Epochs   Data Loops       "+
						"Distance Moved        Gradient Norm");
			}
		}
		
		// Outer loop over epochs
		for (int m=0; m<M; m++) {
			// Clear phi
			Phi = new HashSet<Integer>(eIters*bs);
			avg_x = new double[d];
			avg_g = new double[d];
			subBindings = new ArrayList[maxBindings];
			for (int currM=0; currM<maxBindings; currM++) {
				subBindings[currM] = new ArrayList<Integer>();
			}
			// Loop over epoch iterates
			for (int t=0; t<eIters; t++) {
				// Update average stores
				for (int currIdx=0; currIdx<d; currIdx++) {
					avg_x[currIdx] += x_t[currIdx];
					avg_g[currIdx] += g_t[currIdx];
				}
				totIters++;						// For ADAM update
				// Set batchsize and store indices in Phi
				model.sampleBatch(bs);
				batchIdx = model.currBatchIdx;	// Pass by reference
				f_it = stochasticEvaluate(x_t).gradientVector;
				// Optimize for bs = 1
				if (bs==1) {
					// Simple operations
					Phi.add(batchIdx[0]);
					a_it = stochasticEvaluate(theta_m[theta_m_binding[batchIdx[0]]]).gradientVector;
					// Compute SVRG gradient
					for (int idx=0; idx<d; idx++) {
						g_t[idx] = f_it[idx] - a_it[idx] + aBar_m[idx];
					}
				} else {
					// Begin by adding all new points
					for (int idx=0; idx<bs; idx++) {
						Phi.add(batchIdx[idx]);
					}
					a_it = new double[d];
					// Create a new sub binding array
					for (int currM=0; currM<maxBindings; currM++) {
						subBindings[currM].clear();
					}
					// Re-assign values in subBinding array
					for (int idx=0; idx<bs; idx++) {
						subBindings[theta_m_binding[batchIdx[idx]]].add(batchIdx[idx]);
					}
					// Now run stochastic evaluation on every theta_m batch
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
			// Compute xBar, letting S_l = l (or eIters in this case)
			xBar = Array.scalarMultiply(avg_x, 1.0/((double) eIters));
			// Need to set theta_m1_binding, compute min/max_Tidx, etc.
			// Start by deep copying theta_m_binding and computing the min
			PhiArray = new ArrayList<Integer>(Phi);
			currPhiSize = Phi.size();
			theta_m1_binding = Array.clone(theta_m_binding);
			// Now update the theta_m value of the recently evaluated points
			max_TIdx = maxBindings;
			for (int idx : PhiArray) {
				theta_m1_binding[idx] = max_TIdx;
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
			for (int idx=0; idx<PhiArray.size(); idx++) {
				subBindings[theta_m_binding[PhiArray.get(idx)]].add(PhiArray.get(idx));
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
				if (xStar!=null) {
					printStep(m+1, model.evaluatedDataPoints/N, 
							Array.norm(Array.subtract(xBar, x_t_prev)), 
							Array.norm(aBar_m), Array.norm(Array.subtract(xBar, xStar)));
				} else {
					printStep(m+1, model.evaluatedDataPoints/N,
							Array.norm(Array.subtract(xBar, x_t_prev)), 
							Array.norm(aBar_m));
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
		}
		throw new Exception("kSVRG Failure: maximum epochs exceeded without convergence!");
	}
	
	// k2SVRG here
	private Fit run_k2SVRG(double[] seed, String trajectoryFile) throws Exception {
		int totIters = 0, maxBindings, min_TIdx, max_TIdx, currPhiSize;
		double gNorm, divisor, v_tSum, a_t;
		int[] theta_m_binding, theta_m1_binding, ind, batchIdx, Phi, bsInd=null;
		double[] x0, x_t, avg_x, avg_g, x_t_prev, g_t, m_t, v_t, f_it, a_it;
		double[] aBar_m, xBar;
		double[][] theta_m, theta_m1;
		ArrayList<Integer> PhiArray;
		ArrayList<Integer>[] subBindings;
		
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
		m_t 	= new double[d];
		v_t 	= new double[d];
		theta_m = new double[1][d];
		maxBindings = 1;
		theta_m[0] = Array.clone(x0);
		theta_m_binding = new int[N];
		theta_m1_binding= new int[N];
		// aBar_m must be set to the gradient on the whole dataset
		aBar_m 	= gradientEval(x0, true).gradientVector;
		
		// Print reporting line
		if (isVerbose) {
			if (xStar!=null) {
				System.out.println("    Epochs   Data Loops       Distance "+
						"Moved        Gradient Norm       Distance to x*");
			} else {
				System.out.println("    Epochs   Data Loops       "+
						"Distance Moved        Gradient Norm");
			}
		}
		
		// Outer loop over epochs
		for (int m=0; m<M; m++) {
			// Create a random permutation of the points
			ind = fisherYates();
			if (bs>1) {
				bsInd = fisherYates();
			}
			// Loop over sub-ks
			for (int j=0; j<k; j++) {
				// Zero-out averages
				avg_x = new double[d];
				avg_g = new double[d];
				subBindings = new ArrayList[maxBindings];
				for (int currM=0; currM<maxBindings; currM++) {
					subBindings[currM] = new ArrayList<Integer>();					
				}
				
//				System.out.println("In sub-k "+j+" for epoch m="+m);
				
				// Loop over epoch iterates
				for (int t=0; t<eIters; t++) {
//					System.out.println("Current t iterate: "+t);
//					Array.print(x_t);
					// Update average stores
					for (int currIdx=0; currIdx<d; currIdx++) {
						avg_x[currIdx] += x_t[currIdx];
						avg_g[currIdx] += g_t[currIdx];
					}
					totIters++;						// For ADAM update
					// Optimize for bs = 1
					if (bs==1) {
						model.sampleBatch(1);
						f_it = stochasticEvaluate(x_t).gradientVector;
						a_it = stochasticEvaluate(theta_m[theta_m_binding[model.currBatchIdx[0]]]).gradientVector;
						// Compute SVRG gradient
						for (int idx=0; idx<d; idx++) {
							g_t[idx] = f_it[idx] - a_it[idx] + aBar_m[idx];
						}
					} else {
//						System.out.println("Updating current batch idx");
						// Update batch idx directly
						model.currBatchIdx = Arrays.copyOfRange(bsInd, 
								j*eIters*bs+t*bs, j*eIters*bs+(t+1)*bs);
						batchIdx = model.currBatchIdx;
						model.currBatchSize = bs;
//						System.out.print("bIdx:\t");
//						Array.print(batchIdx);
//						System.out.println(model.currBatchSize);
//						System.out.print("x_t:\t");
//						Array.print(x_t);
						f_it = stochasticEvaluate(x_t).gradientVector;
//						System.out.print("f_it:\t");
//						Array.print(f_it);
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

				// Phi Array
				Phi = Arrays.copyOfRange(ind, j*eIters*bs, (j+1)*eIters*bs);
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
					if (xStar!=null) {
						printStep(m+1, model.evaluatedDataPoints/N, 
								Array.norm(Array.subtract(xBar, x_t_prev)), 
								Array.norm(aBar_m), Array.norm(Array.subtract(xBar, xStar)));
					} else {
						printStep(m+1, model.evaluatedDataPoints/N,
								Array.norm(Array.subtract(xBar, x_t_prev)), 
								Array.norm(aBar_m));
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
}