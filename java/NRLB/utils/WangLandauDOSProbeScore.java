package utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Random;

import base.Sequence;
import base.Shape;
import model.MultinomialResults;
import model.Round0Model;

public class WangLandauDOSProbeScore{
	static int R0k				= 2;
	static int k, windows, varLen;
	static double nsBinding		= 0;
	static long[] splitSeq;
	static ScoringObject so;
	static Random generator;
	static char[] validBases	= {'A', 'C', 'G', 'T'};
	
	public static void main(String[] args) {
		int nInit, minCount, nDivs, I0, In, iters, modelIdx, currIdx, currCount;
		long x0Long, xnLong, currSeq;
		double Emin, Emax, divSize, tol, fCrit, E0, En, a, f=1, Z, currE;
		double Zeff = 0, maxValue = 0;
		String currLine, homeDir, resultsPath, R0Path, probesPath;
		String x0, xn;
		boolean[] RH;
		int[] H;
		long[] probeCounts, obsCounts;
		double[] S;
		Shape shapeModel;
		generator			= new Random();
		Round0Model R0;
		BufferedReader br;
		MultinomialResults r;		
		
		//Parse command line arguments
		if (args.length<1) {	
			throw new IllegalArgumentException("Please provide arguments: homeDir, modelRelativePath, modelIndex, R0modelRelativePath, "
					+ "Emin, Emax, divSize, nInit, tol, minCount, fCrit, probesPath");
		}
		
		homeDir = args[0];
		resultsPath = homeDir + args[1];
		modelIdx = Integer.parseInt(args[2]);
		R0Path = homeDir + args[3];
		Emin = Double.parseDouble(args[4]);
		Emax = Double.parseDouble(args[5]);
		divSize = Double.parseDouble(args[6]);
		nInit = Integer.parseInt(args[7]);
		tol = Double.parseDouble(args[8]);
		minCount = Integer.parseInt(args[9]);
		fCrit = Double.parseDouble(args[10]);
		probesPath = homeDir + args[11];
		
		//Load model and create new scoring object
		r = new MultinomialResults(resultsPath);
		varLen = r.l;
		shapeModel = null;
		R0 = new Round0Model(R0Path, varLen, R0k, true, null, null);
		so = new ScoringObject(r, modelIdx, -1, false, false, shapeModel, R0);
		Z = so.getZ();
				
		//Initialize variables
		nDivs = (int) ((Emax-Emin)/divSize + 1);
		H = new int[nDivs];
		RH = new boolean[nDivs];
		S = new double[nDivs];
		probeCounts = new long[nDivs];
		obsCounts = new long[nDivs];
		
		x0 = "";
		for (int i=0; i<varLen; i++) {
			x0 = x0+"A";
		}
		x0Long = new Sequence(x0, 0, varLen).getValue();
		E0 = Math.log10(R0.getSeqProb(x0Long)*so.seqEval(x0Long)/Z);
		I0 = (int) Math.ceil((E0+Emax)/divSize);
		
		for (int i=0; i<(nDivs*nInit); i++) {
			// propose new state
			xn = randomSequence(x0);
			xnLong = new Sequence(xn, 0, varLen).getValue();
			En = Math.log10(R0.getSeqProb(xnLong)*so.seqEval(xnLong)/Z);
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
		x0Long = new Sequence(x0, 0, varLen).getValue();
		E0 = Math.log10(R0.getSeqProb(x0Long)*so.seqEval(x0Long)/Z);
		I0 = (int) Math.ceil((E0+Emax)/divSize);
		iters = 0;
		
		while (f > tol) {
			// propose new state
			xn = randomSequence(x0);
			xnLong = new Sequence(xn, 0, varLen).getValue();
			En = Math.log10(R0.getSeqProb(xnLong)*so.seqEval(xnLong)/Z);
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
		System.err.println("Density of States computed; evaluating probability.");
		
		//First find smallest value in S that is nonzero:
		for (int i=1; i<nDivs; i++) {
			if ((S[i]>maxValue) && RH[i]) {
				maxValue = S[i];
			}
		}
		//Rescale Entropy S to DOS:
		for (int i=0; i<nDivs; i++) {
			if (RH[i]) {
				S[i] = Math.exp(S[i] - maxValue+.01);
			}
		}
		//Compute an effectiveZ via trapezoidal rule:
		for (int i=1; i<nDivs; i++) {
			Zeff += (S[i]+S[i-1])/2*divSize;
		}
		//Normalize S into a probability density:
		for (int i=0; i<nDivs; i++) {
			S[i] /= Zeff;
		}
		
		System.err.println("Computing observed sequencing rate.");
		
		//Compute total number of probes in each bin
		String[] splitString;
		for (int i=1; i<nDivs; i++) {
			probeCounts[i] = (long) ((S[i]+S[i-1])/2*divSize*Math.pow(4, varLen));
		}
		
		//Loop over all probes in file and sort output
		try {
			br = new BufferedReader(new FileReader(probesPath));
			while ((currLine = br.readLine())!=null) {
				splitString = currLine.split("\t");
				currCount = Integer.parseInt(splitString[2]);
				if (currCount>0) {
					currSeq = new Sequence(splitString[0], 0, varLen).getValue();
					currE = Math.log10(R0.getSeqProb(currSeq)*so.seqEval(currSeq)/Z);
					currIdx = (int) Math.ceil((currE+Emax)/divSize);
					obsCounts[currIdx] += currCount;
				}
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		System.err.println("Binning complete.");
		for (int i=0; i<nDivs; ++i) {
			System.out.printf("%.3f\t%d\t%d\n",Emin+i*divSize,probeCounts[i],obsCounts[i]);
		}
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
