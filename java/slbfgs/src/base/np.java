package base;

import Jama.Matrix;
import java.util.Arrays;

public class np {
	
	// Matrix transpose
	public static double[][] t(double[][] input) {
		double[][] output = new double[input[0].length][input.length];
		
		for (int x=0; x<input.length; x++) {
			for (int y=0; y<input[0].length; y++) {
				output[y][x] = input[x][y];
			}
		}
		return output;
	}
	
	// Vector dot product
	public static double dot(double[]a, double[]b) {
		double output = 0;
		for (int i=0; i<a.length; i++) {
			output += a[i]*b[i];
		}
		return output;
	}
	
	// Element-wise vector square
	public static double[] square(double[] a) {
		double[] output = new double[a.length];
		for (int i=0; i<a.length; i++) {
			output[i] = a[i]*a[i];
		}
		return output;
	}
	
	// Negate all elements in vector
	public static double[] neg(double[] a) {
		double[] output = new double[a.length];
		for (int i=0; i<a.length; i++) {
			output[i] = -1*a[i];
		}
		return output;
	}
	
	// Replace with Array.cat
/*	public static double[] append(double[] a, double b) {
		double[] output = new double[a.length + 1];
		for (int i=0; i<a.length; i++) {
			output[i] = a[i];
		}
		output[a.length] = b;
		return output;
	}
	
	public static double[] append(double[] a, double[] b) {
		double[] output = new double[a.length + b.length];
		for (int i=0; i<a.length; i++) {
			output[i] = a[i];
		}
		for (int i=0; i<b.length; i++) {
			output[i+a.length] = b[i];
		}
		return output;
	} */
	
	// Create vector or array of zeros (no need for this)
/*	public static double[] zeros(int n){
		double[] output = new double[n];
		return output;
	} 
	
	
	public static double[][] zeros(int n, int m){
		double[][] output = new double[n][m];
		return output;
	} */
	
	// Create vector of ones
	public static double[] ones(int n){
		double[] output = new double[n];
		for (int i=0; i<n; i++) {
			output[i] = 1;
		}
		return output;
	}

	public static double[] array(double a) {
		double[] output = new double[1];
		output[0] = a; 
		return output;
	}
	
	public static double[] minimum(double[] a, double[] b) {
		double[] output = new double[a.length];
		for (int i=0; i<a.length; i++) {
			output[i] = Math.min(a[i], b[i]);
		}
		return output;
	}
	
	public static double[] minimum(double a, double[] b) {
		double[] output = new double[b.length];
		for (int i=0; i<b.length; i++) {
			output[i] = Math.min(a, b[i]);
		}
		return output;
	}
	
	public static double[] pow(double[] a, double pow) {
		double[] output = new double[a.length];
		for (int i=0; i<a.length; i++) {
			output[i] = Math.pow(a[i], pow);
		}
		return output;
	}
	
	public static double[] sum(double[] a, double b) {
		double[] output = new double[a.length];
		for (int i=0; i<a.length; i++) {
			output[i] = a[i] + b;
		}
		return output;
	}
	
	public static double[] sum(double[] a, double[] b) {
		double[] output = new double[a.length];
		for (int i=0; i<a.length; i++) {
			output[i] = a[i] + b[i];
		}
		return output;
	}
	
	public static double[] absolute(double[] a) {
		double[] output = new double[a.length];
		for (int i=0; i<a.length; i++) {
			output[i] = Math.abs(a[i]);
		}
		return output;
	}
	
	public static double[] mult(double[] a, double m) {
		double[] output = new double[a.length];
		for (int i=0; i<a.length; i++) {
			output[i] = a[i]*m;
		}
		return output;
	}
	
	public static double[] mult(double[] a, double[] b) {
		double[] output = new double[a.length];
		for (int i=0; i<a.length; i++) {
			output[i] = a[i]*b[i];
		}
		return output;
	}

	public static double[] less(double[] a, double[] b) {
		double[] output = new double[a.length];
		for (int i=0; i<a.length; i++) {
			if(a[i] < b[i]) {
				output[i] = 1;
			} else {
				output[i] = 0;
			}
		}
		return output;
	}
	
	public static double[] notless(double[] a, double[] b) {
		double[] output = new double[a.length];
		for (int i=0; i<a.length; i++) {
			if(a[i] < b[i]) {
				output[i] = 0;
			} else {
				output[i] = 1;
			}
		}
		return output;
	}
	
	public static double[][] diag(double[] d) {
		double[][] output = new double[d.length][d.length];
		for(int i = 0; i < d.length; i++) {
			output[i][i] = d[i];
		}
		return output;
	}
	
	public static double[][] block(double[][] a, double[][] b, double[][] c, double[][] d) {
		int n = a.length + c.length;
		int m = a[0].length + b[0].length;
		double[][] output = new double[n][m];
		
		for(int i = 0; i < a.length; i++) {
			for(int j = 0; j < a[0].length; j++) {
				output[i][j] = a[i][j];
			}
		}
		
		for(int i = 0; i < a.length; i++) {
			for(int j = 0; j < b[0].length; j++) {
				output[i][a[0].length + j] = b[i][j];
			}
		}
		
		for(int i = 0; i < c.length; i++) {
			for(int j = 0; j < c[0].length; j++) {
				output[a.length + i][j] = c[i][j];
			}
		}
		
		for(int i = 0; i < d.length; i++) {
			for(int j = 0; j < d[0].length; j++) {
				output[a.length + i][c[0].length + j] = d[i][j];
			}
		}
		return output;
	}

	public static double[][] add(double[][] a, double[][] b) {
		// TODO Auto-generated method stub
		int n = a.length;
		int m = a[0].length;
		
		double[][] output = new double[n][m];
		
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < m; j++) {
				output[i][j] = a[i][j] + b[i][j];
			}
		}
		
		return output;
	}
	
	public static double[][] matmul(double[][] A, double[][] B){
		int n = A.length;
		int m = B[0].length;
		int x = A[0]. length;
		double[][] output = new double[n][m];
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < m; i++) {
				for(int k = 0; k < x; k++) {
					output[i][j] += A[i][k]*B[k][j];
				}
			}
		}
		return output;
	}
	
	public static double dot(double[] B, double[][] A){
		if(A[0].length>1 || A.length != B.length) {
			System.out.println("Error: Matrix multiplication dimensions");
		}
		double output = 0;
		for (int i=0; i<B.length; i++) {
			output += A[i][0]*B[i];
		}
		return output;
	}
	
	public static double dot(double[] B, Matrix Am){
		double[][] A = Array.clone(Am.getArray());
		if(A[0].length>1 || A.length != B.length) {
			System.out.println("Error: Matrix multiplication dimensions");
		}
		double output = 0;
		for (int i=0; i<B.length; i++) {
			output += A[i][0]*B[i];
		}
		return output;
	}


	public static Matrix getMatrix(double[] array_in) {
		// TODO Auto-generated method stub
		double[][] array_out = new double[array_in.length][1];    
	    for (int i=0; i<array_in.length; i++) {
	    	array_out[i][0] = array_in[i];
	    }
	    Matrix matrix_out = new Matrix(array_out);
	    return matrix_out;
	}
	
	public static double[] min(double[] a) {
		double[] output = new double[2];
		double min = a[0];
		double argmin = 0;
		for (int i=1; i<a.length; i++) {
			if(a[i] < min) {
				min = a[i];
				argmin = i;
			}
		}
		output[0] = min;
		output[1] = argmin;
		return output;
	}

	// Append array to the bottom of a matrix
	public static double[][] stack(double[][] a, double[] b) {
		int lenA = (a==null) ? 0 : a.length;
		int lenB = (a==null) ? b.length : a[0].length;
		double[][] output = new double[lenA + 1][lenB];
		for (int i=0; i<lenA; i++) {
			for (int j=0; j<lenB; j++) {
				output[i][j] = a[i][j];
			}
		}
		for (int j=0; j<lenB; j++) {
			output[lenA][j] = b[j];
		}
		return output;
	}

	// Sort array
	public static double[] sort(double[] a) {
		// TODO Auto-generated method stub
		double[] output = Array.clone(a);
		Arrays.sort(output);
		return output;
	}

	public static int argmin(double[] a) {
		double min = a[0];
		int argmin = 0;
		for (int i=1; i<a.length; i++) {
			if(a[i] < min) {
				min = a[i];
				argmin = i;
			}
		}
		return argmin;
	}
	
	public static double max(double[] a) {
		double max = a[0];
		for (int i=1; i<a.length; i++) {
			if(a[i] > max) {
				max = a[i];
			}
		}
		return max;
	}
	
	public static int argmax(double[] a) {
		double max = a[0];
		int argmax = 0;
		for (int i=1; i<a.length; i++) {
			if(a[i] > max) {
				max = a[i];
				argmax = i;
			}
		}
		return argmax;
	}
}
