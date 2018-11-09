package test;

import base.Fit;
import base.Model;

public class MinimizerFit extends Fit {
	private static final long serialVersionUID = -6491336134655349557L;
	public int functionType;							//The type of function this Functions Object will represent
	public int nDim;									//The number of dimensions for this space
	public double[] seed;


	public MinimizerFit(int functionType, int nDim, double[] seed) {		//Function constructor
		this.functionType	= functionType;				//Stores values for the type of function
		this.nDim			= nDim;						//And the number of dimensions
		this.seed			= seed;
	}
	
	@Override
	public boolean equals(Object o) {
		return false;
	}

	@Override
	public int hash() {
		return 0;
	}

	@Override
	public void merge(Fit newFit) {
		//do nothing
	}

	@Override
	public void recordFit(int fitSteps, int functionCalls, double fitTime, double trainLikelihood, Model input) {
		this.fitSteps = fitSteps;
		this.functionCalls = functionCalls;
		this.fitTime = fitTime;
		this.trainLikelihood = trainLikelihood;
		this.finalPosition = input.getPositionVector();
	}

	@Override
	public void recordErrorBars(int nullVectors, double[] input) {
		//do nothing
	}
}
