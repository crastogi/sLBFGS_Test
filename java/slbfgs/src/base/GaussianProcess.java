package base;
import java.util.Arrays;

import Jama.*;

public class GaussianProcess {
	int N = 0;
	public double c1, c2, f0, beta, sigmaf, sigmadf, alpha0, offset = 10;
	public double m0 = 0, dm0 = 0, V0 = 0, Vd0 = 0, dVd0 = 0;
	Matrix A, G;
	public double[] T, ordT, F, Y, dY_projected, Sigmaf;
	public double[][] X, dY, Sigmadf;

	// Basic constructor; requires input of starting function value and scaling
	// This allows the GP to be 'naturally scaled'
	public GaussianProcess(double f0, double beta, double alpha0) {
		this.f0		= f0;
		this.beta	= beta;
		this.alpha0	= alpha0;
		c1 			= .0001;		// Default values
		c2 			= .99;
	}
	
	// Change default c1, c2 values
	public GaussianProcess(double c1, double c2, double f0, double beta, 
			double alpha0) {
		this.f0		= f0;
		this.beta	= beta;
		this.alpha0	= alpha0;
		this.c1 	= c1;
		this.c2 	= c2;
	}
	
	public void updateGP(double t, double[] x, double f, double[] df, 
			double dfProj, double var_f, double[] var_df, double var_dfProj) {
		double tempSigma;
		
		// Store tested positions, function and adj. function values
		T			= Array.cat(T, t);
		X			= np.stack(X, x);
		F			= Array.cat(F, f);
		Y			= Array.cat(Y, (f-f0)/beta);
		
		// Store gradient and projected gradient values
		dY			= np.stack(dY, df);
		dY_projected= Array.cat(dY_projected, dfProj/beta);
		
		// Store noise values
		Sigmaf		= Array.cat(Sigmaf, var_f);
		Sigmadf		= np.stack(Sigmadf, var_df);
		// Update Noise, if needed
		tempSigma	= Math.sqrt(var_f)/beta;
		if (tempSigma > sigmaf) {
			sigmaf 	= tempSigma;
		}
		tempSigma	= Math.sqrt(var_dfProj)/beta;
		if (tempSigma > sigmadf) {
			sigmadf	= tempSigma;
		}
		
		// Increment total storage
		N++;
		
		// Create sorted T
		ordT = Array.clone(T);
		Arrays.sort(ordT);
		
		// Build GP
		buildGP();
	}
	
	// Helper function to build GP
	private void buildGP() {
		// build Gram matrix
	    double[][] kTT   = new double[N][N];
	    double[][] kdTT  = new double[N][N]; 
	    double[][] dkdTT = new double[N][N];
	    
	    for (int i=0; i<N; i++) {
			for (int j=0; j<N; j++) {
				kTT[i][j]   = k(T[i],  T[j]);
				kdTT[i][j]  = kd(T[i], T[j]);
				dkdTT[i][j] = dkd(T[i], T[j]);
			}
	    }
	    
	    // build noise matrix
	    double[] Sig = new double[2*N];
	    for(int j=0; j<N; j++) {
	    	Sig[j] = Math.pow(sigmaf,2);
	    }
	    for(int j=N; j<2*N; j++) {
	    	Sig[j] = Math.pow(sigmadf, 2);
	    }
	    
	    // build Gram matrix
	    double[][] G_array = np.add(np.diag(Sig), np.block(kTT,kdTT,np.t(kdTT),dkdTT));
	    
	    G = new Matrix(G_array);
	    
	    double[] Y_array = Array.cat(Y, dY_projected);
	    double[][] Y_matrix = new double[Y_array.length][1];
	    
	    for (int i=0; i<Y_array.length; i++) {
	    	Y_matrix[i][0] = Y_array[i];
	    }
	    
	    Matrix Y_m = new Matrix(Y_matrix);    
	    A = G.solve(Y_m);
	    
	    m0   = m(0.);
	    dm0  = d1m(0.);
	    V0   = V(0.);
	    Vd0  = Vd(0.);
	    dVd0 = dVd(0.);
	}

	public double k(double a, double b) {
		double k = (1./3.) * Math.pow(Math.min(a+offset,b+offset),3) + 0.5 * Math.abs(a-b) * Math.pow(Math.min(a+offset,b+offset),2);
		return k;
	}
	
	public double[] k(double[] A, double[] B) {
		double[] k = np.sum(np.mult(np.pow(np.minimum(np.sum(A, offset),np.sum(B, offset)),3), 1./3.), np.mult(np.mult(np.absolute(np.sum(A, np.neg(B))),np.pow(np.minimum(np.sum(A, offset),np.sum(B, offset)),2)),0.5));
		return k;
	}
	
	public double[] kd(double[] a, double[] b) {
	    double[] kd;
	    kd = np.sum(np.mult(np.less(a,b),np.mult(np.pow(np.sum(a,offset),2),1./2.)), 
	    		np.mult(np.notless(a,b) , np.sum(np.mult(np.sum(a,offset),np.sum(b,offset)) ,np.mult(np.pow(np.sum(b,offset), 2), -0.5))));
	    return kd;
	}
	
	
	public double kd(double a, double b) {
	    double kd;
	    if(a<b) {
	    	kd = Math.pow(a+offset,2)/2.;
	    } else {
	    	kd = (a+offset)*(b+offset) - 0.5 * Math.pow(b+offset,2);
	    }
	    return kd;
	}
	
	private double dkd(double a, double b) {
		return Math.min(a+offset,b+offset);
	}
	
	private double[] dkd(double a, double[] b) {
		double[] dkd =  new double[b.length];
		for (int i=0; i<b.length; i++) {
			dkd[i] = dkd(a,b[i]);
		}
		return dkd;
	}
	
	public double[] k(double a, double[] b) {
		double[] k =  new double[b.length];
		for (int i=0; i<b.length; i++) {
			k[i] = k(a,b[i]);
		}
		return k;
	}
	
	public double[] kd(double a, double[] b) {
		double[] kd =  new double[b.length];
		for (int i=0; i<b.length; i++) {
			kd[i] = kd(a,b[i]);
		}
		return kd;
	}
	
	public double[] dk(double a, double[] b) {
		double[] dk =  new double[b.length];
		for (int i=0; i<b.length; i++) {
			dk[i] = dk(a,b[i]);
		}
		return dk;
	}
	
	public double dk(double a, double b) {
		double dk =  0;
		if(a > b) {
			dk = Math.pow(b+offset,2)/2.;
		} else {
			dk = (a+offset)*(b+offset) - 0.5 * Math.pow(a+offset,2);
		}		
		return dk;
	}
	
	public double[] ddk(double a, double[] b) {
		double[] ddk =  new double[b.length];
		for (int i=0; i<b.length; i++) {
			if(a <= b[i]){
				ddk[i] = b[i] - a;
			}
		}
		return ddk;
	}
	
	public double[] ddkd(double a, double[] b) {
		double[] ddkd =  new double[b.length];
		for (int i=0; i<b.length; i++) {
			if(a <= b[i]){
				ddkd[i] = 1;
			}
		}
		return ddkd;
	}
	
	public double[] dddk(double a, double[] b) {
		double[] dddk =  new double[b.length];
		for (int i=0; i<b.length; i++) {
			if(a <= b[i]){
				dddk[i] = -1;
			}
		}
		return dddk;
	}
	
	public double m(int idx) {
		return m(T[idx]);
	}
	
	public double m(double t) {
		double m = np.dot(Array.cat(k(t, T), kd(t, T)), A);
		return m;
	}
	
	public double d1m(int idx) {
		return d1m(T[idx]);
	}
	
	public double d1m(double t) {
		double d1m = np.dot(Array.cat(dk(t, T), dkd(t, T)),A);
		return d1m;
	}
	
	public double d2m(double t) {
		double d2m = np.dot(Array.cat(ddk(t, T), ddkd(t, T)),A);
		return d2m;
	}
	
	public double d3m(double t) {
		double d3m = np.dot(Array.cat(dddk(t, T), new double[N]),A);
		return d3m;
	}
	
	public double V(double t) {
		Matrix second = G.solve(np.getMatrix(Array.cat(k(t, T) ,kd(t, T))));
		double V = k(t,t)  - np.dot(Array.cat(k(t, T) ,  kd(t, T)), second);
		return V;
	}
	
	public double Vd(double t) {
		Matrix second = G.solve(np.getMatrix(Array.cat(dk(t, T) ,  dkd(t, T))));
		double Vd = kd(t,t)  - np.dot(Array.cat(k(t, T) ,  kd(t, T)), second);
		return Vd;
	}
	
	public double dVd(double t) {
		Matrix second = G.solve(np.getMatrix(Array.cat(dk(t, T) ,  dkd(t, T))));
		double dVd = dkd(t,t)  - np.dot(Array.cat(dk(t, T) ,  dkd(t, T)),second);
		return dVd;
	}
	
	public double V0f(double t) {
		Matrix second = G.solve(np.getMatrix(Array.cat(k(t, T) , kd(t, T))));
		double V0f = k(0,t)  - np.dot(Array.cat(k(0, T) ,  kd(0, T)),second);
		return V0f;
	}
	
	public double Vd0f(double t) {
		Matrix second = G.solve(np.getMatrix(Array.cat(k(t, T) , kd(t, T))));
		double Vd0f = dk(0,t)  - np.dot(Array.cat(dk(0, T) ,  dkd(0, T)),second);
		return Vd0f;
	}
	
	public double V0df(double t) {
		Matrix second = G.solve(np.getMatrix(Array.cat(dk(t, T) , dkd(t, T))));
		double V0df = kd(0,t)  - np.dot(Array.cat(k(0, T) ,  kd(0, T)),second);
		return V0df;
	}
	
	public double Vd0df(double t) {
		Matrix second = G.solve(np.getMatrix(Array.cat(dk(t, T) , dkd(t, T))));
		double Vd0df = dkd(0,t)  - np.dot(Array.cat(dk(0, T) ,  dkd(0, T)),second);
		return Vd0df;
	}
	
	public double probWolfe(int idx) {
		return probWolfe(T[idx]);
	}
	
	public double probWolfe(double t) {	    
		double p;
	    // marginal for Armijo condition
	    double ma = m0 - m(t) + c1*t*dm0;
	    double Vaa = V0 + Math.pow(c1*t,2) * dVd0 + V(t) + 2*(c1*t*(Vd0 - Vd0f(t)) - V0f(t));

	    // marginal for curvature condition
	    double mb  = d1m(t) - c2*dm0;
	    double Vbb = Math.pow(c2,2)*dVd0 - 2*c2*Vd0df(t) + dVd(t);
	    
	    // covariance between conditions
	    double Vab = -c2*(Vd0 + c1*t*dVd0) + V0df(t) + c2*Vd0f(t) + c1*t*Vd0df(t) - Vd(t);
	    
	    // deterministic evaluations
	    if (Vaa < 1e-9 && Vbb < 1e-9) {
	    	if(ma >= 0 && mb >=0) {
	    		p = 1.;
	    	} else {
	    		p = 0.;
	    	}
	    	return p;
	    }
	    
	    double rho = Vab / Math.sqrt(Vaa * Vbb);
	    
	    if (Math.abs(rho) > 1) {
	    	rho = Math.signum(rho) * 1.;
	    }
	    if (Vaa <= 0 || Vbb <= 0) {
	    	p  = 0;
	    	return p;
	    }
	    
	    double upper = (2*c2*(Math.abs(dm0)+2*Math.sqrt(dVd0))-mb)/Math.sqrt(Vbb);
	    
	    double xl = -ma/Math.sqrt(Vaa);
	    double xu = Double.POSITIVE_INFINITY;
	    double yl = -mb/Math.sqrt(Vbb);
	    double yu = upper;

	    p = bvn.cdf(xu, yu, rho) + bvn.cdf(xl, yl, rho) - bvn.cdf(xu, yl, rho) - bvn.cdf(xl, yu, rho);
	    return p;
	}
}
