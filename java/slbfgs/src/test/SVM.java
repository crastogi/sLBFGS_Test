package test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

import base.Array;
import base.Model;

public class SVM extends Model {
	public double lambda;								// Regularizer for svm
	public double[] theta;								// The vector storing the parameters for current evaluation
	public double[] classes;							// The class of every data point
	public double[][] data;								// The data for every point

	public SVM(double lambda) {						// Load data
		int maxCols;
        String csvFile = "./src/featurized_df.csv";
        String line;
        double[] temp;
        String[] split;
        ArrayList<double[]> parsed = new ArrayList<double[]>();

        // Read in data file
        try (BufferedReader br = new BufferedReader(new FileReader(csvFile))) {
        	// Ignore first line
        	line = br.readLine();
        	maxCols = line.split(",").length;
            // Read, split, and ensure that the remaining lines have the same number of columns
        	// Store the remaining values in a matrix
        	while ((line = br.readLine()) != null) {
            	split = line.split(",");
            	if (split.length != maxCols) {
            		throw new Exception("Variable number of columns!");
            	}
            	temp = new double[maxCols];
            	for (int i=0; i<maxCols; i++) {
            		temp[i] = Double.parseDouble(split[i]);
            	}
            	parsed.add(temp);
            }
        	// Store values in data
        	N = 10000; //32560; //10000 
        	data = new double[N][maxCols-1];
        	classes = new double[N];
        	for (int i=0; i<N; i++) {
        		for (int j=0; j<maxCols-1; j++) {
        			data[i][j] = parsed.get(i)[j];
        		}
        		classes[i] = parsed.get(i)[maxCols-1];
        	}
        	nDim = maxCols-1;
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        this.lambda = lambda;
        theta = new double[nDim];
        fName = "Simple SVM";
	}

	// Computes both the likelihood and the gradient
	public CompactGradientOutput evaluate() {
		double alpha, y;
		double likelihood = 0, temp;
		double x[];
		double gradient[] = new double[nDim];
		
		for (int i=0; i<N; i++) {
			x = data[i];
			y = classes[i];
			alpha = Array.dotProduct(x, theta);
			temp = Math.max(0.0, 1 - y*alpha);
			likelihood += temp*temp;
			if (y*alpha < 1.0) {
				gradient = Array.addScalarMultiply(gradient, -1, Array.scalarMultiply(data[i], y*(1 - y*alpha)));
			}
		}
		likelihood /= 2*N;											// Normalize by the number of data points
		likelihood += lambda*Array.dotProduct(theta, theta)/2;		// Add regularization
		gradient = Array.scalarMultiply(gradient, 1/((double) N));
		gradient = Array.addScalarMultiply(gradient, lambda, theta);
		evaluatedDataPoints += N;
		
		return (new CompactGradientOutput(likelihood, gradient));
	}
	
	public CompactGradientOutput stochasticEvaluate() {
		double alpha, y, varLik;
		double likelihood = 0, temp, lik2 = 0;
		double x[];
		double gradient[] = new double[nDim];
		double varGrad[] = new double[nDim];
		
		for (int i : currBatchIdx) {
			x = data[i];
			y = classes[i];
			alpha = Array.dotProduct(x, theta);
			temp = Math.max(0.0, 1 - y*alpha);
			temp = temp*temp;
			likelihood += temp;
			lik2 += temp*temp;
			if (y*alpha < 1.0) {
				// Unwrapped gradient update
				temp = y*(1 - y*alpha);
				for (int j=0; j<nDim; j++) {
					gradient[j] += -data[i][j]*temp;
					varGrad[j] += data[i][j]*data[i][j]*temp*temp;
				}
//				gradient = Array.addScalarMultiply(gradient, -1, Array.scalarMultiply(data[i], temp));
			}
		}
		varLik = (lik2-likelihood*likelihood/currBatchIdx.length)/(4*currBatchIdx.length*(currBatchIdx.length-1));
		for (int j=0; j<nDim; j++) {
			varGrad[j] = (varGrad[j]-gradient[j]*gradient[j]/currBatchIdx.length)/(currBatchIdx.length*(currBatchIdx.length-1));
		}
		likelihood /= 2*currBatchIdx.length;						// Normalize by the number of data points
		likelihood += lambda*Array.dotProduct(theta, theta)/2;		// Add regularization
		gradient = Array.scalarMultiply(gradient, 1/((double) currBatchIdx.length));
		gradient = Array.addScalarMultiply(gradient, lambda, theta);
		evaluatedDataPoints += currBatchIdx.length;
		return (new CompactGradientOutput(likelihood, gradient, varLik, varGrad));
	}
	
	public CompactGradientOutput stochasticEvaluate(ArrayList<Integer> idx) {
		double alpha, y;
		double likelihood = 0, temp;
		double x[];
		double gradient[] = new double[nDim];
		
		for (int i=0; i<idx.size(); i++) {
			x = data[idx.get(i)];
			y = classes[idx.get(i)];
			alpha = Array.dotProduct(x, theta);
			temp = Math.max(0.0, 1 - y*alpha);
			likelihood += temp*temp;
			if (y*alpha < 1.0) {
				gradient = Array.addScalarMultiply(gradient, -1, Array.scalarMultiply(data[idx.get(i)], y*(1 - y*alpha)));
			}
		}
		likelihood /= 2*idx.size();						// Normalize by the number of data points
		likelihood += lambda*Array.dotProduct(theta, theta)/2;		// Add regularization
		gradient = Array.scalarMultiply(gradient, 1/((double) idx.size()));
		gradient = Array.addScalarMultiply(gradient, lambda, theta);
		evaluatedDataPoints += idx.size();
		return (new CompactGradientOutput(likelihood, gradient));
	}
	
	public CompactGradientOutput evaluate(double[] theta) {			//Overloaded operator for easy function calling
		setParams(theta);
		return evaluate();
	}
	
	@Override
	public void setParams(double[] position) {
		theta = Array.clone(position);
	}

	@Override
	public double[] getParams() {
		return Array.clone(theta);
	}
}
