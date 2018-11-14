package base;

import java.util.ArrayList;

public abstract class Model {
	public int nDim = 0, N = 0, evaluatedDataPoints = 0;
	public String fName = "N/A";
	public int[] currBatchIdx;
	private ArrayList<int[]> storedBatches = new ArrayList<int[]>();
	private MersenneTwisterFast mtfast = new MersenneTwisterFast();
	
	public int getNFeatures() {
		return nDim;
	}

	public abstract void setParams(double[] position);
	
	public abstract double[] getParams();
 
	public abstract double functionEval() throws Exception;
	
	public abstract CompactGradientOutput gradientEval() throws Exception;
	
	public abstract CompactGradientOutput stochasticGradientEval() throws Exception;
	
	public abstract CompactGradientOutput getGradient();
	
	public Fit generateFit(double[] seed) {
		return new Fit(fName, nDim, seed);
	}
	
	public void sampleBatch(int k) {
		int nBatches, randIdx, temp;
		int[] shuffledIdx, tempBatch;
		
		if (k>N) {
			throw new IllegalArgumentException("Cannot sample a batch that is larger than the dataset!");
		}
		// Ensure that batches exist and that they are of the expected length
		if (storedBatches.size()!=0 && storedBatches.get(0).length==k) {
			currBatchIdx = storedBatches.get(0);
			storedBatches.remove(0);
			return;
		} 
		// Else, generate a new set of batches
		nBatches = (int) Math.floor( ((double) N)/((double) k) );
		shuffledIdx = new int[N];
		
		// Begin by populating the array
		for (int i=0; i<N; i++) {
			shuffledIdx[i] = i;
		}
		
		// Now shuffle
		for (int i=N-1; i>0; i--) {
			randIdx = mtfast.nextInt(i+1);
			// Swap at index
			temp = shuffledIdx[randIdx];
			shuffledIdx[randIdx] = shuffledIdx[i];
			shuffledIdx[i] = temp;
		}
		
		// Create batches and store
		for (int i=0; i<nBatches; i++) {
			tempBatch = new int[k];
			for (int j=0; j<k; j++) {
				tempBatch[j] = shuffledIdx[i*k+j];
			}
			storedBatches.add(tempBatch);
		}
		
		// set idx and exit
		currBatchIdx = storedBatches.get(0);
		storedBatches.remove(0);
		return;
	}
	
	public class CompactGradientOutput {
		public double functionValue = 0;
		public double[] gradientVector;

		public CompactGradientOutput(double functionValue, double[] gradientVector) {
			this.functionValue	= functionValue;
			this.gradientVector	= gradientVector;
		}
		
		public CompactGradientOutput(double[] gradientVector) {
			this.gradientVector	= gradientVector;
		}
	}
}
