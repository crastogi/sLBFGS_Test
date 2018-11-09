package test;

import base.Array;

public class LBFGSMinimizers {
	public int nDim, maxM, maxIterations, maxLSIterations, nLSEvals, nFunctionEvals;
	private int type;
	public int counter=0, counter2=0;
	public double epsilon, uTol, c1, c2, stepAlphaMin, stepAlphaMax;
	private boolean isBracketed;
	private double fCurr, alphaMin, alphaMax, alphaL, fL, gL, fT, gT, alphaU, fU, gU;
	private double[] xCurr, gCurr, pCurr;
	private LBFGSFunctions f;
	private LBFGSFunctions.CompactOutput fOut;
	public double reportGNorm, reportFVal, reportIters;
	
	//LBFGS object constructor; load basic minimization parameters. To be used for all subsequent minimizations using this object.
	public LBFGSMinimizers(LBFGSFunctions f, int maxM, double uTol, double c1, double c2,
			double stepAlphaMin, double stepAlphaMax, double epsilon, int maxIterations, int maxLSIterations, int type) {
		this.f				= f;				//Define the function that needs to be minimized.
		nDim				= f.nDim;			//Get the number of dimensions from the Function object
		if (maxM<=0)	throw new IllegalArgumentException("Maximum memory depth must be positive!");
		this.maxM			= maxM;				//Set the maximum memory depth.
		if (uTol<0)		throw new IllegalArgumentException("Machine tolerance limit must be positive!");
		this.uTol			= uTol;				//Set the machine precision (line search)
		if (c1<0||c1>1)	throw new IllegalArgumentException("The parameter c1 must be between 0 and 1!");
		this.c1				= c1;				//c1 in the strong Wolfe conditions (line search)
		if (c2<=c1)		throw new IllegalArgumentException("The parameter c2 must be larger than c1 in the strong Wolfe conditions!");
		if (c2>1)		throw new IllegalArgumentException("The parameter c2 must be less than 1!");
		this.c2				= c2;				//c2 in the strong Wolfe conditions (line search)
		if (stepAlphaMin<0)	
						throw new IllegalArgumentException("The line search step lower bound must be positive!");
		this.stepAlphaMin	= stepAlphaMin;		//Lower bound for a step in the line search
		if (stepAlphaMax<stepAlphaMin)
						throw new IllegalArgumentException("The line search step upper bound must be greater than the lower bound!");
		this.stepAlphaMax	= stepAlphaMax;		//Upper bound for a step in the line search
		this.epsilon		= epsilon;			//Set the accuracy with which the solution needs to be found.
		if (maxIterations<2)	throw new IllegalArgumentException("The maximum number of minimization steps must be greater than 1!");
		this.maxIterations	= maxIterations;	//Maximum number of minimization steps before convergence failure is declared
		if (maxLSIterations<2)	throw new IllegalArgumentException("The maximum number of line search steps must be greater than 1!");
		this.maxLSIterations= maxLSIterations;	//Maximum number of line search function calls before convergence failure is declared.
		this.type = type;
	}
	
	public void minimize(double[] xInput) throws Exception {
		int k 			= 0;
		nFunctionEvals	= 0;
		int idx;
		double fNext, alphaCurr, gammaCurr, beta;
		double[] xNext, gNext, q;
		double[] h0Curr	= new double[nDim];
		pCurr			= new double[nDim];
		double[] alpha	= new double[maxM];
		double[] r		= new double[nDim];
		double[][] s	= new double[maxM][nDim];
		double[][] y	= new double[maxM][nDim];
		
		xCurr	= Array.clone(xInput);
		fOut	= f.evaluate(xCurr);
		fCurr	= fOut.value;
		gCurr	= fOut.gradient;
		nFunctionEvals++;
		if (Array.norm(gCurr)/Math.max(1, Array.norm(xCurr)) < epsilon) {
			throw new Exception("Already at minimum!");
		}
		h0Curr	= Array.ones(nDim);						//Create a new hk0 array with ones on the diagonal

		for (int i=0; i<nDim; i++) {					//Store the first search direction
			pCurr[i] = -gCurr[i]*h0Curr[i];
		}
		alphaCurr	= (type==1) ? lineSearch(1/Array.norm(gCurr)) : mcsearch(1/Array.norm(gCurr)); //lineSearch(1/Array.norm(gCurr));
		xNext		= addScalarMultiply(xCurr, alphaCurr, pCurr);
		fNext		= fOut.value;
		gNext		= fOut.gradient;
		s[0]		= Array.subtract(xNext, xCurr);
		y[0]		= Array.subtract(gNext, gCurr);
		k++;
		
		while(k < maxIterations) {
			xCurr		= xNext;
			fCurr		= fNext;
			gCurr		= gNext;
			gammaCurr	= Array.dotProduct(s[0], y[0])/Array.dotProduct(y[0], y[0]);
			h0Curr		= Array.scalarMultiply(Array.ones(nDim), gammaCurr);
			
			q	= Array.clone(gCurr);
			for (int i=k-1; i>=Math.max(0, k-maxM); i--) {
				idx = k-1-i;
				alpha[idx] = Array.dotProduct(s[idx], q)/Array.dotProduct(y[idx], s[idx]);
				q = addScalarMultiply(q, -alpha[idx], y[idx]);
			}
			for (int i=0; i<nDim; i++) {
				r[i] = h0Curr[i]*q[i];
			}
			for (int i=Math.max(0, k-maxM); i<=k-1; i++) {
				idx = k-1-i;
				beta = Array.dotProduct(y[idx], r)/Array.dotProduct(y[idx], s[idx]);
				r = addScalarMultiply(r, (alpha[idx]-beta), s[idx]);
			}
			
			pCurr		= Array.scalarMultiply(r, -1.0);
			alphaCurr	= (type==1) ? lineSearch(1) : mcsearch(1); //lineSearch(1);
			xNext		= addScalarMultiply(xCurr, alphaCurr, pCurr);
			fNext		= fOut.value;
			gNext		= fOut.gradient;
			if (Array.norm(gNext)/Math.max(1, Array.norm(xNext)) <= epsilon) {
//				Array.print(xNext);
//				System.out.println("F Value:\t"+fNext+"\tGNorm:\t"+Array.norm(gNext));
//				System.out.println("Iterations:\t"+(k+1)+"\t\t\tTotal Function Calls:\t"+nFunctionEvals);
				reportFVal = fNext;
				reportGNorm= Array.norm(gNext);
				reportIters= k+1;
				return;
			}
			s			= cycleDown(s, Array.subtract(xNext, xCurr));
			y			= cycleDown(y, Array.subtract(gNext, gCurr));
			k++;
//			System.out.println("F Value: "+fNext+" GNorm: "+Array.norm(gNext)+" Iterations: "+k+" currAlpha: "+alphaCurr);
		}
		throw new Exception("LBFGS Failure: maximum iterations exceeded!");
	}
	
	private double mcsearch(double alphaT) throws Exception {
		boolean useModifiedFunc;
		double currIntWidth, prevIntWidth, f0, g0, fTTest;
		double[] x0;
		
		//Check to see if guess alphaT is within bounds
		if (alphaT<stepAlphaMin || alphaT>stepAlphaMax) {
			throw new Exception("Line search failure: initial step alpha guess out of bounds!");
		}
		//Initialize variables and see if search direction is a descent direction
		x0				= xCurr;
		f0				= fCurr;
		g0				= Array.dotProduct(gCurr, pCurr);
		if (g0 > 0) {
			throw new Exception("Line search failure: the search direction is not a descent direction!");
		}
		alphaL			= 0;
		fL				= f0;
		gL				= g0;
		alphaU			= 0;
		fU				= f0;
		gU				= g0;
		currIntWidth	= stepAlphaMax-stepAlphaMin;
		prevIntWidth	= 2*currIntWidth;
		useModifiedFunc	= true;
		isBracketed		= false;
		nLSEvals		= 0;
		
		//Create first interval, and ensure that alphaT lies within the acceptable range of alpha values
		alphaT			= Math.min(Math.max(alphaT, stepAlphaMin), stepAlphaMax);
		alphaMin		= alphaL;
		alphaMax		= alphaT + 4*(alphaT-alphaL);
		
		while (nLSEvals<=maxLSIterations) {
			//Evaluate function
			fOut= f.evaluate(addScalarMultiply(x0, alphaT, pCurr));
			fT	= fOut.value;
			gT	= Array.dotProduct(fOut.gradient, pCurr);
			nLSEvals++;
			nFunctionEvals++;
			
			//Control for NaN: rescale pCurr to be of unit length
			if (Double.isNaN(fT)) {
				System.out.println("NaN ERROR!");
				pCurr	= Array.normalize(pCurr);
				fOut	= f.evaluate(addScalarMultiply(x0, alphaT, pCurr));
				fT		= fOut.value;
				gT		= Array.dotProduct(fOut.gradient, pCurr);
				nLSEvals++;
				nFunctionEvals++;
			}
			
			//Check for convergence
			fTTest = f0+alphaT*c1*g0;
			if (alphaT==stepAlphaMax && (fT-fTTest<=0 && gT-c1*g0<0)) { 			// \psi(\alpha_t) \leq 0 and \psi^'(\alpha_t) < 0
				throw new Exception("Line search failure: search terminated, step alpha at maximum!");
			} else if (alphaT==stepAlphaMin && (fT-fTTest>0 || gT-c1*g0>=0)) {		// \psi(\alpha_t) > 0 and \psi^'(\alpha_t) \geq 0
				throw new Exception("Line search failure: search terminated, optimal step alpha is lower than minimum!");
			} if (fT-fTTest<=0 && Math.abs(gT)<=c2*Math.abs(g0)) {
				return alphaT;														//converged under strong Wolfe conditions
			}
			
			//Check if modified update should be used
			if (useModifiedFunc && fT-fTTest<=0 && gT>=0) {							// \psi(\alpha_t) \leq 0 and \phi^'(\alpha_t)>0
				useModifiedFunc = false;
			}
			//Generate a new safeguarded alphaT and interval I
			if (useModifiedFunc) {
				alphaT = mcstep(alphaL, fL-alphaL*c1*g0, gL-c1*g0, alphaT, fT-alphaT*c1*g0, gT-c1*g0, 
						alphaU, fU-alphaU*c1*g0, gU-c1*g0);
			} else {
				alphaT = mcstep(alphaL, fL, gL, alphaT, fT, gT, alphaU, fU, gU);
			}
			
			//Ensure that the interval I decreases sufficiently using a bisection prediction  
			if (isBracketed && (Math.abs(alphaU-alphaL)>=.66*prevIntWidth) ) {
				alphaT = alphaL + (alphaU-alphaL)/2;
			}
			alphaT = Math.min(Math.max(alphaT, stepAlphaMin), stepAlphaMax);
			//Define bounds of new interval
			if (isBracketed) {
				alphaMin = Math.min(alphaL, alphaU);
				alphaMax = Math.max(alphaL, alphaU);
			} else {
				alphaMin = alphaL;
				alphaMax = alphaT+4*(alphaT-alphaL);
			}
			prevIntWidth = currIntWidth;
			currIntWidth = alphaMax-alphaMin;
			//Check to see that the interval is larger than uTol (machine precision)
			if (isBracketed && currIntWidth<=uTol*alphaMax) {
				throw new Exception("Line search failure: The relative width of the interval of uncertainty is less than uTol!");
			}
			//If unusual termination is about to occur, let alphaT be the lowest point found so far
			if ( (isBracketed && (alphaT<=alphaMin || alphaT>=alphaMax)) || nLSEvals>=maxLSIterations-1) {
				System.out.println("UNUSUAL TERMINATION!");
				alphaT = (alphaU+alphaL)/2;
			}
		}
		throw new Exception("Line search failure: number of line search function calls exceeded maximum!");
	}
	
	private double mcstep(double alphaLM, double fLM, double gLM, double alphaTM, double fTM, double gTM, 
			double alphaUM, double fUM, double gUM) throws Exception {
		boolean isBound = false;
		double d1, d2, s, p, q, cubic, quadratic, stpf;
		double derivativeSign = Math.signum(gTM)*Math.signum(gLM);
		
		if (gLM*(alphaTM-alphaLM)>=0.0) {
			throw new Exception("Interval cannot contain a minimizer!"); //Violates Theorem 2.1
		}
		
		//Case 1: Higher function value. Corresponds with case U1
		if(fTM>fLM) {
			isBracketed = true;
			d1 = gLM + gTM + 3*(fLM-fTM)/(alphaTM-alphaLM);
			s = Math.max(Math.abs(d1), Math.max(Math.abs(gLM), Math.abs(gTM)));
			d2 = s*Math.sqrt((d1/s)*(d1/s) - (gLM/s)*(gTM/s));
			if (alphaTM<alphaLM) d2 = -d2;
			p = d1 + d2 - gLM;
			q = gTM + 2*d2 - gLM;
			cubic = alphaLM + p/q*(alphaTM-alphaLM);
			quadratic = alphaLM + (alphaTM-alphaLM)*(gLM/( ( (fLM-fTM)/(alphaTM-alphaLM) ) +gLM) )/2;
			if ( Math.abs(cubic-alphaLM) < Math.abs(quadratic-alphaLM) ) {
				stpf = cubic;
			} else {
				stpf = cubic + (quadratic-cubic)/2;
			}
		} //Case 2: lower function value with derivatives of opposite sign. Corresponds with case U3
		else if(derivativeSign<0) {	
			isBracketed = true;
			d1 = gLM + gTM + 3*(fLM-fTM)/(alphaTM-alphaLM);
			s = Math.max(Math.abs(d1), Math.max(Math.abs(gLM), Math.abs(gTM)));
			d2 = s*Math.sqrt((d1/s)*(d1/s) - (gLM/s)*(gTM/s));
			if (alphaTM>alphaLM) d2 = -d2;
			p = d1 + d2 - gTM;
			q = gLM + 2*d2 - gTM;
			cubic = alphaTM + p/q*(alphaLM - alphaTM);
			quadratic = alphaTM + (gTM/(gTM-gLM))*(alphaLM-alphaTM);
			if (Math.abs(cubic-alphaTM) > Math.abs(quadratic-alphaTM)) {
				stpf = cubic;
			} else {
				stpf = quadratic;
			}
		} //Case 3: lower function value, derivatives of the same sign, with decreasing magnitude. Case U2
		else if(Math.abs(gTM)<Math.abs(gLM)) {	
			isBound = true;
			d1 = gLM + gTM + 3*(fLM-fTM)/(alphaTM-alphaLM);
			s = Math.max(Math.abs(d1), Math.max(Math.abs(gLM), Math.abs(gTM)));
			d2 = s*Math.sqrt(Math.max(0, ((d1/s)*(d1/s) - (gLM/s)*(gTM/s)) ));
			if (alphaTM>alphaLM) d2 = -d2;
			p = d1 + d2 - gTM;
			q = gLM + 2*d2 - gTM;
			if (p/q<0 && d2!=0) {		//If cubic tends to infinity in the direction of the step
				cubic = alphaTM + p/q*(alphaLM-alphaTM);
			} 
			else if (alphaTM>alphaLM) {	//If the cubic estimation will be out of bounds, restrict it
				cubic = alphaMax;
			} else {
				cubic = alphaMin;
			}
			quadratic = alphaTM + (gTM/(gTM-gLM))*(alphaLM-alphaTM);
			if(isBracketed) {
				if(Math.abs(alphaTM-cubic) < Math.abs(alphaTM-quadratic)) {
					stpf = cubic;
				} else {
					stpf = quadratic;
				}
			} //Since the interval has not been bracketed, along with case U2
			//the following conditions satisfy safeguarding conditions 2.2+2.3 and 2.4+2.5
			else {
				if (alphaTM>alphaLM) {
					stpf = alphaMax;
				} else {
					stpf = alphaMin;
				}
			}
		} //Case 4: Lower function value, derivatives have same sign without decreasing magnitude. Case U2
		else {
			if (isBracketed) {
				d1 = gUM + gTM + 3*(fTM-fUM)/(alphaUM-alphaTM);
				s = Math.max(Math.abs(d1), Math.max(Math.abs(gUM), Math.abs(gTM)));
				d2 = 2*Math.sqrt((d1/s)*(d1/s) - (gUM/s)*(gTM/s));
				if (alphaTM>alphaUM) d2 = -d2;
				p = d1 + d2 - gTM;
				q = gUM + 2*d2 - gTM;
				cubic = alphaTM + p/q*(alphaUM-alphaTM);
				stpf = cubic;
			} //Since the interval has not been bracketed, along with case U2
			//the following conditions satisfy safeguarding conditions 2.2+2.3 and 2.4+2.5
			else if (alphaTM>alphaLM) {
				stpf = alphaMax;
			} else {
				stpf = alphaMin;
			}
		}
		
		//Update interval of uncertainty (Updating algorithm)
		if (fTM>fLM) {
			alphaU = alphaTM;
			fU = fTM;
			gU = gTM;
		} else {
			if (derivativeSign<0) {
				alphaU = alphaLM;
				fU = fLM;
				gU = gLM;
			}
			alphaL = alphaTM;
			fL = fTM;
			gL = gTM;
		}
		
		//Compute new safeguarded step
		stpf = Math.min(Math.max(stpf, alphaMin), alphaMax);
		if (isBracketed && isBound) {	//Modified version of bounds in case3?
			if (alphaUM>alphaLM) {
				stpf = Math.min(alphaLM + 0.66*(alphaUM-alphaLM), stpf);
			} else {
				stpf = Math.max(alphaLM + 0.66*(alphaUM-alphaLM), stpf);
			}
		}
		return stpf;
	}
	
	private double lineSearch(double alphaNew) throws Exception {
		double fNew, gNew, alphaOld, fOld, gOld, f0, g0;
		double[] x0;
		
		if (alphaNew<stepAlphaMin || alphaNew>stepAlphaMax) {			//Check to see if guess is within bounds
			throw new Exception("Line search failure: initial step alpha guess out of bounds!");
		}
		alphaOld= 0;
		x0		= xCurr;
		fOld	= fCurr;
		f0		= fCurr;
		gOld	= Array.dotProduct(gCurr, pCurr);
		g0		= gOld;
		nLSEvals= 1;
		if (gOld>0) {								//Ensure search direction is a descent direction
			throw new Exception("Line search failure: the search direction is not a descent direction!");
		}
		
		while (nLSEvals<=maxLSIterations) {
			fOut = f.evaluate(addScalarMultiply(x0, alphaNew, pCurr));
			fNew = fOut.value;
			gNew = Array.dotProduct(fOut.gradient, pCurr);
			nFunctionEvals++;
			
			if ( (fNew>fOld+c1*alphaNew*gOld) || ((fNew>=fOld) && nLSEvals>1) ) {
				return zoom(alphaOld, fOld, gOld, alphaNew, fNew, gNew, x0, f0, g0);
			} else if (Math.abs(gNew)<=-c2*g0) {
				return alphaNew;
			} else if (gNew>=0) {
				return zoom(alphaNew, fNew, gNew, alphaOld, fOld, gOld, x0, f0, g0);
			}
			alphaOld= alphaNew;
			fOld	= fNew;
			gOld	= gNew;
			alphaNew= 10*alphaNew;
			if (alphaNew>stepAlphaMax || alphaNew<stepAlphaMin) {
				throw new Exception("Line search failure: step alpha guess is out of bounds!");
			}
			nLSEvals++;
		}
		throw new Exception("Line search failure: number of line search function calls exceeded maximum!");		
	}
	
	public double zoom(double alphaLow, double fLow, double gLow, double alphaHigh, double fHigh, double gHigh, 
			double[] x0, double f0, double g0) throws Exception {
		double d1, d2, s, alphaJ, fJ, gJ;
		
		while (nLSEvals<=maxLSIterations) {
			//Cubic interpolation
			d1	= gHigh + gLow + 3*(fLow - fHigh)/(alphaHigh-alphaLow);
			s	= Math.max(Math.abs(d1), Math.max(Math.abs(gLow), Math.abs(gHigh)));
			d2	= s*Math.sqrt( ((d1/s)*(d1/s) - (gLow/s)*(gHigh/s)) );
			if (alphaHigh<alphaLow)	d2 = -d2;
			alphaJ = alphaLow + (alphaHigh-alphaLow)*(d1+d2-gLow)/(-gLow+gHigh+2*d2);
			if (alphaJ>stepAlphaMax || alphaJ<stepAlphaMin) {
				throw new Exception("Line search failure: step alpha guess is out of bounds! alphaJ: "+alphaJ);
			}
			
			fOut= f.evaluate(addScalarMultiply(x0, alphaJ, pCurr));
			fJ	= fOut.value;
			gJ	= Array.dotProduct(fOut.gradient, pCurr);
			nLSEvals++;
			nFunctionEvals++;
			if( (fJ>f0+c1*alphaJ*g0) || (fJ>=fLow) ) {
				alphaHigh	= alphaJ;
				fHigh		= fJ;
				gHigh		= gJ;
			} else {
				if (Math.abs(gJ)<=-c2*g0) {
					return alphaJ;
				} else if (gJ*(alphaHigh-alphaLow)>=0) {
					alphaHigh	= alphaLow;
					fHigh		= fLow;
					gHigh		= gLow;
				}
				alphaLow= alphaJ;
				fLow	= fJ;
				gLow	= gJ;
			}
		}
		throw new Exception("Line search failure: number of line search function calls exceeded maximum!");	
	}
	
	public double[][] cycleDown(double[][] input, double[] newRow) {	//cycle rows in matrix downwards and add a new row on top
		double[][] output = new double[maxM][nDim];
		
		for (int i=maxM-1; i>0; i--) {
			for (int j=0; j<nDim; j++) {
				output[i][j] = input[i-1][j];
			}
		}
		for (int i=0; i<nDim; i++) {
			output[0][i] = newRow[i];
		}
		return output;
	}
	
	public double[] addScalarMultiply(double[] x, double c, double[] y) {
		double[] output = new double[nDim];
		
		for (int i=0; i<nDim; i++) {
			output[i] = x[i]+c*y[i];
		}
		return output;
	}
}