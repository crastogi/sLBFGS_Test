package utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import base.Array;
import base.Sequence;
import base.Shape;
import model.MultinomialFit;
import model.MultinomialResults;
import model.Round0Model;

public class ModelProfiler2{
	static String R0Results		= "/Users/chaitanya/Documents/Research/SELEX/MultinomialPaper/ModelComparisons/DeepBind/AUROC/FittedModels/CEBPd/CEBPD_TAATGA20NCG_Y_R0.dat";
	static int R0k				= 2;
	
	public static void main(String[] args) {
		int k, nChars, windows, modelIdx, nBases;
		double nsBinding = 0;
		long[] splitSeq;
		String currLine, homeDir, resultsPath, probesPath;
		Shape shapeModel;
		Round0Model R0 		= new Round0Model(R0Results, 20, R0k, true, null, null);
		BufferedReader br;
		MultinomialFit fit;
		MultinomialResults r;
		ScoringObject so;

		//Parse command line arguments
		if (args.length<1) {	
			throw new IllegalArgumentException("Please provide arguments: homeDir, modelRelativePath, modelIndex, probesPath, nBases");
		}
		homeDir = args[0];
		resultsPath = homeDir + args[1];
		modelIdx = Integer.parseInt(args[2]);
		probesPath = homeDir + args[3];
		nBases = Integer.parseInt(args[4]);
		
		//Load model
		r = new MultinomialResults(resultsPath);
		fit	= r.getFit(modelIdx);
		if (!fit.isMulti) {
			k = fit.k;
		} else {
			k = fit.ks[0];
		}
		if (fit.isNSBinding) {
			nsBinding = Math.exp(fit.finalPosition[fit.finalPosition.length-1])/(92+nBases-k+1);
		}
		r.l = k;
		fit.isFlank = false;
		fit.isNSBinding = false;
		shapeModel = null;
		
		//Create scoring object
		so = new ScoringObject(r, modelIdx, -1, false, false, shapeModel, R0);
		
		//Read in probes and score
		try {
			br = new BufferedReader(new FileReader(probesPath));
			currLine = br.readLine();
			nChars = currLine.length();
			windows= nChars-k+1;
			nsBinding = nsBinding*windows;
			splitSeq = new long[windows];
			
			System.out.println(args[3]);
			//Print the score of the first sequence
			for (int j=0; j<windows; j++) {
				splitSeq[j] = (new Sequence(currLine.substring(j, j+k), 0, k)).getValue();
			}
			System.out.println(Math.log(Array.sum(so.seqEval(splitSeq))+nsBinding));
			
			//Print the remainder of the sequences 
			while ((currLine = br.readLine())!=null) {
				for (int j=0; j<windows; j++) {
					splitSeq[j] = (new Sequence(currLine.substring(j, j+k), 0, k)).getValue();
				}
				System.out.println(Math.log(Array.sum(so.seqEval(splitSeq))+nsBinding));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
