package test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

import base.Array;
import base.Model;

public class SVMTest extends Model {
	public double lambda;								// Regularizer for svm
	public double[] theta;								// The vector storing the parameters for current evaluation
	public double[] classes;							// The class of every data point
	public double[][] data;								// The data for every point

	public SVMTest(double lambda) {						// Load data
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
        	N = 10000; //parsed.size();
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
		
		return (new CompactGradientOutput(likelihood, gradient));
	}
	
	public CompactGradientOutput stochasticEvaluate() {
		double tStart = System.nanoTime();
		
		double alpha, y;
		double likelihood = 0, temp;
		double x[];
		double gradient[] = new double[nDim];
		
		for (int i : currBatchIdx) {
			x = data[i];
			y = classes[i];
			alpha = Array.dotProduct(x, theta);
			temp = Math.max(0.0, 1 - y*alpha);
			likelihood += temp*temp;
			if (y*alpha < 1.0) {
				gradient = Array.addScalarMultiply(gradient, -1, Array.scalarMultiply(data[i], y*(1 - y*alpha)));
			}
		}
		likelihood /= 2*currBatchIdx.length;						// Normalize by the number of data points
		likelihood += lambda*Array.dotProduct(theta, theta)/2;		// Add regularization
		gradient = Array.scalarMultiply(gradient, 1/((double) currBatchIdx.length));
		gradient = Array.addScalarMultiply(gradient, lambda, theta);
		
		StochasticTime += (System.nanoTime()-tStart)/1E9; 
		
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
