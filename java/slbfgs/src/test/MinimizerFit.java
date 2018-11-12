package test;

import base.Fit;
import base.Model;

public class MinimizerFit extends Fit {
	public int functionType;							//The type of function this Functions Object will represent
	public int nDim;									//The number of dimensions for this space
	public double[] seed;

	public MinimizerFit(int functionType, int nDim, double[] seed) {		//Function constructor
		this.functionType	= functionType;				//Stores values for the type of function
		this.nDim			= nDim;						//And the number of dimensions
		this.seed			= seed;
	}

	@Override
	public void recordFit(int fitSteps, int functionCalls, double fitTime, double likelihood, Model input) {
		this.fitSteps = fitSteps;
		this.functionCalls = functionCalls;
		this.fitTime = fitTime;
		this.likelihood = likelihood;
		this.finalPosition = input.getPositionVector();
	}
}
