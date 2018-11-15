package base;

import java.util.HashSet;
import java.util.Iterator;

public abstract class Model {
	
	public double SampleTime = 0;
	public double StochasticTime = 0;
	
	public int nDim = 0, N = 0, evaluatedDataPoints = 0;
	public String fName = "N/A";
	public int[] currBatchIdx;
	private MersenneTwisterFast mtfast = new MersenneTwisterFast();
	
	public int getNFeatures() {
		return nDim;
	}

	public abstract void setParams(double[] position);
	
	public abstract double[] getParams();
 
	public abstract CompactGradientOutput evaluate() throws Exception;
		
	public abstract CompactGradientOutput stochasticEvaluate() throws Exception;
	
	public Fit generateFit(double[] seed) {
		return new Fit(fName, nDim, seed);
	}
	
	public void sampleBatch(int k) {
		double tStart = System.nanoTime();
		
		int idx = 0;
		currBatchIdx = new int[k];
		HashSet<Integer> h = new HashSet<Integer>();
		
		for (int i=0; i<k; i++) {
			h.add(mtfast.nextInt(N));
		}
		while(h.size()<k-1) {
			h.add(mtfast.nextInt(N));
		}
		Iterator<Integer> i = h.iterator(); 
        while (i.hasNext()) {
        	currBatchIdx[idx] = i.next();
        	idx++;
        }
        
		SampleTime += (System.nanoTime()-tStart)/1E9;
		
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
