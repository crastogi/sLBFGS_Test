package base;

import java.util.HashSet;
import java.util.Iterator;

public abstract class Model {
	public int nDim = 0, N = 0, evaluatedDataPoints = 0;
	public String fName = "N/A";
	public int[] currBatchIdx;
	public int[] currSuperBatchIdx = null;
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
	
	public void setSuperBatch(int k) {
		int idx = 0;
		currSuperBatchIdx = new int[k];
		HashSet<Integer> h = new HashSet<Integer>();
		
		for (int i=0; i<k; i++) {
			h.add(mtfast.nextInt(N));
		}
		while(h.size()<k) {
			h.add(mtfast.nextInt(N));
		}
		Iterator<Integer> i = h.iterator(); 
        while (i.hasNext()) {
        	currSuperBatchIdx[idx] = i.next();
        	idx++;
        }
		return;
	}
	
	public void sampleBatch(int k) {
		int idx = 0, upperBound = 0;
		currBatchIdx = new int[k];
		HashSet<Integer> h = new HashSet<Integer>();	
		
		if (currSuperBatchIdx!=null) {
			upperBound = currSuperBatchIdx.length;
			if (k>=upperBound) {
				currBatchIdx = Array.clone(currSuperBatchIdx);
				return;
			}
			for (int i=0; i<k; i++) {
				h.add(currSuperBatchIdx[mtfast.nextInt(upperBound)]);
			}
			while(h.size()<k) {
				h.add(currSuperBatchIdx[mtfast.nextInt(upperBound)]);
			}
		} else {
			for (int i=0; i<k; i++) {
				h.add(mtfast.nextInt(N));
			}
			while(h.size()<k) {
				h.add(mtfast.nextInt(N));
			}
		}
		Iterator<Integer> i = h.iterator(); 
        while (i.hasNext()) {
        	currBatchIdx[idx] = i.next();
        	idx++;
        }
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
