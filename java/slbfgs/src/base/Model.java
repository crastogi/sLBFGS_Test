package base;

public abstract class Model {
	public int maxFitRetries	= 5;	
	
	//The user must check to see if o is the right data type
	public abstract void replaceData(Object o);
	
	public abstract double likelihoodNormalizer();
	
	public abstract double maxLikelihood();
	
	public abstract int getNFeatures();

	//setParams does NOT SYMMETRIZE input
	public abstract void setParams(double[] position);
	
	public abstract double[] getPositionVector();

	//None of the following 5 functions utilize or obey symmetries. 
	public abstract double functionEval() throws Exception;
	
	public abstract CompactGradientOutput gradientEval() throws Exception;
	
	public abstract CompactGradientOutput getGradient();
	
	public abstract Fit generateFit(double[] seed);
	
	//Central difference method
	//startLoc is FULL dimensionality (NOT symmetry reduced)
	public double[] gradientFiniteDifferences(double[] startLoc, double stepSize) {
		int totFeatures		= this.getNFeatures();					//Operate in the FULL dimensionality space
		double forwardDifference, reverseDifference;
		double[] baseBetas, modBetas;
		double[] fdGradient	= new double[totFeatures];
		
		baseBetas = Array.clone(startLoc);
		//compute function values
		for (int i=0; i<totFeatures; i++) {
			try {
				modBetas			= Array.clone(baseBetas);
				modBetas[i]			+= stepSize;
				forwardDifference	= functionEval();
				modBetas[i]			-= 2*stepSize;
				reverseDifference	= functionEval();
				fdGradient[i]		= (forwardDifference-reverseDifference)/(2*stepSize);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		setParams(startLoc);
		return fdGradient;
	}
	
	//Central difference method
	public double[][] hessianFiniteDifferences(double[] startLoc, double stepSize) {
		int totFeatures 			= this.getNFeatures();
		double[] baseBetas, modBetas;
		double[][] forwardDifference= new double[totFeatures][totFeatures];
		double[][] reverseDifference= new double[totFeatures][totFeatures];
		double[][] fdHessian		= new double[totFeatures][totFeatures];
        CompactGradientOutput output;
		
		baseBetas = Array.clone(startLoc);
		//compute gradients
		for (int i=0; i<totFeatures; i++) {
			try {
				modBetas			= Array.clone(baseBetas);
				modBetas[i]			+= stepSize;
				setParams(modBetas);
				output				= gradientEval();
				//Handle the case where the gradient has not been defined
				if (output==null)	throw new UnsupportedOperationException(
						"The function gradientEval() has not been defined!");
				forwardDifference[i]= output.gradientVector;
				modBetas[i]			-= 2*stepSize;
				setParams(modBetas);
				output				= gradientEval();
				reverseDifference[i]= output.gradientVector;
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		//Find finite differences (forward/backward FD)
		for (int i=0; i<totFeatures; i++) {
			for (int j=0; j<totFeatures; j++) {
				fdHessian[i][j] = (forwardDifference[j][i]-reverseDifference[j][i])/(4*stepSize) + 
									(forwardDifference[i][j]-reverseDifference[i][j])/(4*stepSize);
			}
		}
		setParams(startLoc);
		return fdHessian;
	}
	
	public class CompactGradientOutput {
		public double functionValue		= 0;
		public double[] gradientVector;
		public double[][] hessian;

		public CompactGradientOutput(double functionValue, double[] gradientVector) {
			this.functionValue	= functionValue;
			this.gradientVector	= gradientVector;
		}
		
		public CompactGradientOutput(double[] gradientVector) {
			this.gradientVector	= gradientVector;
		}
		
		public CompactGradientOutput(double functionValue, double[] gradientVector, double[][] hessian) {
			this.functionValue	= functionValue;
			this.gradientVector	= gradientVector;
			this.hessian		= hessian;
		}
		
		public CompactGradientOutput(double[] gradientVector, double[][] hessian) {
			this.gradientVector	= gradientVector;
			this.hessian		= hessian;
		}
		
		public CompactGradientOutput(double[][] hessian) {
			this.hessian		= hessian;
		}
	}
}
