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

public class ShapeModelTester {
	static int l				= 16;
	static int k				= 10;
	static String R0Results		= "/vega/hblab/users/cr2166/SELEX/chaitanya/Round0/scr_CCACGTC.dat";
	static int R0k				= 6;	
	static String resultsPath	= "/vega/hblab/users/cr2166/SELEX/chaitanya/scripts/full_shape_sweep/output/";
	static String kmerPath		= "/vega/hblab/users/cr2166/SELEX/chaitanya/scripts/full_shape_sweep/Scr_R1_kmers.txt";
	static String shapePath		= "/vega/hblab/users/cr2166/SELEX/chaitanya/scripts/ShapeFolder/Orthogonal/";
	static String outputPath	= "/vega/hblab/users/cr2166/SELEX/chaitanya/scripts/full_shape_sweep/output/";
	
	public static void main(String[] args) {
		double sum, Z;
		String currLine;
		ScoringObject so;
		Shape shapeModel;
		MultinomialResults r;
		BufferedReader br;
		Round0Model R0 		= new Round0Model(R0Results, l, R0k, true, null, null);
		long[] tempSeqs;
		ArrayList<Long> seqs = new ArrayList<Long>();
		PrintStream outputFile, original = System.out;
		
		//First read in all kmers to be evaluated
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
		
		//Loop over files
		for (int currFileIdx = 1; currFileIdx<=1; currFileIdx++) {
			//Load new results and shape tables
			r = new MultinomialResults(resultsPath+"Scr-orthogonal.dat");
			shapeModel = new Shape(shapePath, true);
			//Loop over different fits
			for (int currFitIdx = 0; currFitIdx<r.nFits(); currFitIdx++) {
				//Load fit and create new scoring object
				so = new ScoringObject(r, currFitIdx, -1, false, false, shapeModel, R0);
				Z = so.getZ();
				//Loop over all kmers and print counts to a file
				try {
					outputFile = new PrintStream(new FileOutputStream(outputPath+"Orthogonal_"+currFitIdx+".txt"));					
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
