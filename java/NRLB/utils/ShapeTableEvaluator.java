package utils;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import base.Sequence;
import base.Shape;
import model.MultinomialResults;
import model.Round0Model;

public class ShapeTableEvaluator {
	static int l			= 16;     //length of SELEX-seq probes
	static int k			= 10;     //length of kmers used in tables
	static String basePath, R0Results, resultsPath, kmerPath, shapePath, outputPath;
	static int R0k;
	
	public static void main(String[] args) {
		double sum, Z;
		String currLine;
		ScoringObject so;
		Shape shapeModel;
		MultinomialResults r;
		BufferedReader br;
		Round0Model R0;
		long[] tempSeqs;
		ArrayList<Long> seqs	= new ArrayList<Long>();
		PrintStream outputFile, original = System.out;
		
		//Read in arguments: base path, R0Results, R0k, resultsPath, kmerPath, shapePath, outputPath 
		if (args.length==7) {
			basePath	= args[0];			//Base file path
			R0Results	= basePath+args[1];		//Round 0 Model relative path
			R0k		= Integer.parseInt(args[2]);	//Length of Round 0 model
			resultsPath	= basePath+args[3];		//NRLB model relative path
			kmerPath	= basePath+args[4];		//10mer table relative path
			shapePath	= basePath+args[5];		//Shape table relative path
			outputPath	= basePath+args[6];		//Relative path for output files
		} else {
			throw new IllegalArgumentException("Improper number of arguments supplied!");
		}
		//Load R0 model, shape tables, and NRLB model
		R0		= new Round0Model(R0Results, l, R0k, true, null, null);
		shapeModel	= new Shape(shapePath, true);
		r		= new MultinomialResults(resultsPath);		

		//Read in all kmers to be evaluated
		try {
			br = new BufferedReader(new FileReader(kmerPath));
			while ((currLine = br.readLine())!=null) {
				seqs.add((new Sequence(currLine, 0, k)).getValue());
			}
			br.close();
			System.out.println("Finished loading kmers.");
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//Loop over different fits
		for (int currFitIdx = 0; currFitIdx<r.nFits(); currFitIdx++) {
			//Load fit and create new scoring object
			so = new ScoringObject(r, currFitIdx, -1, false, false, shapeModel, R0);
			Z = so.getZ();
			//Loop over all kmers and print counts to a file
			try {
				outputFile = new PrintStream(new FileOutputStream(outputPath+currFitIdx+".txt"));					
				//Loop over all probes
				System.setOut(outputFile);
				for (int i=0; i<seqs.size(); i++) {
					sum = 0;
					tempSeqs = seqBuilder(seqs.get(i));
					for (long currSeq : tempSeqs) {
						sum += R0.getSeqProb(currSeq)*so.seqEval(currSeq)/Z;
					}
					sum *= r.totCounts;
					System.out.println((new Sequence(seqs.get(i), k)).getString()+"\t"+sum);
				}
				System.setOut(original);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	public static long[] seqBuilder(long input) {
		int maxIter		= (int) Math.pow(4, l-k)*(l-k+1);
		int groupIter	= (int) Math.pow(4, l-k);
		int addLeft, addRight, currIter = 0;
		long base, lMask, rMask;
		long[] output	= new long[maxIter];
		
		for (int currOffset=0; currOffset<=l-k; currOffset++) {
			addLeft = currOffset;
			addRight= l-k-addLeft;
			lMask	= (long) Math.pow(4, addLeft)-1;
			rMask	= (long) Math.pow(4, addRight)-1;
			base	= input << 2*addRight;
			lMask	= lMask << 2*addRight;
			for (long i=0; i<groupIter; i++) {
				output[currIter] = (((i & lMask) << 2*(k)) | base) | (i & rMask);
				currIter++;
			}
		}
		return(output);
	}
}
