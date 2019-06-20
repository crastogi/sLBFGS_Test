package utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Random;

import base.Array;
import base.Sequence;
import base.Shape;
import model.MultinomialFit;
import model.MultinomialResults;
import model.Round0Model;

public class WangLandauDOS{
	static String R0Results		= "/Users/chaitanya/Documents/Research/SELEX/MultinomialPaper/ModelComparisons/DeepBind/AUROC/FittedModels/CEBPd/CEBPD_TAATGA20NCG_Y_R0.dat";
	static int R0k				= 2;
	static int k, windows, varLen;
	static double nsBinding		= 0;
	static long[] splitSeq;
	static ScoringObject so;
	static Random generator;
	static char[] validBases	= {'A', 'C', 'G', 'T'};
	
	public static void main(String[] args) {
		int nInit, minCount, nDivs, I0, In, iters;
		int nChars, modelIdx, nBases;
		double Emin, Emax, divSize, tol, fCrit, E0, En, a, f=1;
		String currLine, homeDir, resultsPath, probesPath;
		String lFlank, rFlank, x0, xn;
		boolean[] RH;
		int[] H;
		double[] S;
		Shape shapeModel;
		generator			= new Random();
		Round0Model R0 		= new Round0Model(R0Results, 20, R0k, true, null, null);
		BufferedReader br;
		MultinomialFit fit;
		MultinomialResults r;
		
		//Parse command line arguments
		if (args.length<1) {	
			throw new IllegalArgumentException("Please provide arguments: homeDir, modelRelativePath, modelIndex, nBases, "
					+ "varLen, lFlank, rFlank, Emin, Emax, divSize, nInit, tol, minCount, fCrit");
		}
		
		homeDir = args[0];
		resultsPath = homeDir + args[1];
		modelIdx = Integer.parseInt(args[2]);
		nBases = Integer.parseInt(args[3]);
		varLen = Integer.parseInt(args[4]);
		lFlank = args[5];
		rFlank = args[6];
		Emin = Double.parseDouble(args[7]);
		Emax = Double.parseDouble(args[8]);
		divSize = Double.parseDouble(args[9]);
		nInit = Integer.parseInt(args[10]);
		tol = Double.parseDouble(args[11]);
		minCount = Integer.parseInt(args[12]);
		fCrit = Double.parseDouble(args[13]);
		
		//Load model
		r = new MultinomialResults(resultsPath);
		fit	= r.getFit(modelIdx);
		if (!fit.isMulti) {
			k = fit.k;
		} else {
			k = fit.ks[0];
		}
		r.l = k;
		nChars = varLen+lFlank.length()+rFlank.length();
		windows= nChars-k+1;
		splitSeq = new long[windows];
		if (fit.isNSBinding) {
			nsBinding = Math.exp(fit.finalPosition[fit.finalPosition.length-1])/(92+nBases-k+1)*windows;
		}
		System.out.println(Math.log(nsBinding));
		fit.isFlank = false;
		fit.isNSBinding = false;
		shapeModel = null;
		
		//Create scoring object
		so = new ScoringObject(r, modelIdx, -1, false, false, shapeModel, R0);
		
		//Initialize variables
		nDivs = (int) ((Emax-Emin)/divSize + 1);
		H = new int[nDivs];
		RH = new boolean[nDivs];
		S = new double[nDivs];
		
		x0 = "";
		for (int i=0; i<varLen; i++) {
			x0 = x0+"A";
		}
		E0 = scoreSeq(lFlank+x0+rFlank);
		I0 = (int) Math.ceil((E0+Emax)/divSize);
		
		for (int i=0; i<(nDivs*nInit); i++) {
			// propose new state
			xn = randomSequence(x0);
			En = scoreSeq(lFlank+xn+rFlank);
			In = (int) Math.ceil((En+Emax)/divSize);
			// calculate acceptance
			a = Math.exp(S[I0] - S[In]);
			if (generator.nextDouble() <= a) {
				// accept current proposal
				x0 = xn;
				E0 = En;
				I0 = In;
			}
			// perform updates
			S[I0] += f;
			H[I0] += 1;
		}
		// Set reference histogram and reset H, S
		for (int i=0; i<nDivs; ++i) {
			RH[i] = (H[i]>0) ? true : false;
			H[i] = 0;
			S[i] = 0;
		}
		System.err.println("Reference Histogram Computed for n = "+nDivs*nInit+" samples.");
		
		// Begin DOS computation
		x0 = randomSequence(x0);
		E0 = scoreSeq(lFlank+x0+rFlank);
		I0 = (int) Math.ceil((E0+Emax)/divSize);
		iters = 0;
		
		while (f > tol) {
			// propose new state
			xn = randomSequence(x0);
			En = scoreSeq(lFlank+xn+rFlank);
			In = (int) Math.ceil((En+Emax)/divSize);
			// calculate acceptance
			a = Math.exp(S[I0] - S[In]);
			if (generator.nextDouble() <= a) {
				// accept current proposal
				x0 = xn;
				E0 = En;
				I0 = In;
			}
			// perform updates
			S[I0] += f;
			H[I0] += 1;
			iters += 1;
			if (min(H, RH) > minCount) {
				if (min(H, RH) >= fCrit*mean(H, RH)) {
					System.err.println("Flatness criteria reached for f = "+f+" in "+iters+" iterations.");
					f = f/2;
					iters = 0;
					for (int i=0; i<nDivs; ++i) {
						H[i] = 0;
					}
				}
			}
		}
		System.err.println("Density of States computed. Entropy S is given below:");
		for (int i=0; i<nDivs; ++i) {
			System.out.printf("%.3f\t%f\n",Emin+i*divSize,S[i]);
		}
	}
	
	private static double scoreSeq(String seq) {
		for (int j=0; j<windows; j++) {
			splitSeq[j] = (new Sequence(seq.substring(j, j+k), 0, k)).getValue();
		}
		return (Math.log(Array.sum(so.seqEval(splitSeq))+nsBinding));
	}
	
	private static String randomSequence(String inSeq) {
		int changePos = generator.nextInt(varLen);
		char currNuc  = inSeq.charAt(changePos);
		char newNuc   = validBases[generator.nextInt(4)];
		
		while (currNuc==newNuc) {
			newNuc = validBases[generator.nextInt(4)];
		}
		
		char[] temp   = inSeq.toCharArray();
		temp[changePos] = newNuc;
		
		return String.valueOf(temp);
	}
	
	private static double mean(int[] H, boolean[] RH) {
		double out = 0, count = 0;
		
		for (int i=0; i<H.length; i++) {
			if (RH[i]) {
				out += H[i];
				count += 1;
			}
		}
		return (out/count);
	}

	private static int min(int[] H, boolean[] RH) {
		int minval = Integer.MAX_VALUE;
		
		for (int i=0; i<H.length; i++) {
			if (RH[i]) {
				if ( H[i] < minval) {
					minval = H[i];
				}
			}
		}
		return minval;
	}
}
