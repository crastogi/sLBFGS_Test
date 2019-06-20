package utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import base.Sequence;
import base.Shape;
import model.MultinomialFit;
import model.MultinomialResults;
import model.Round0Model;

public class ShapeModelTester2 {
	static int l				= 16;
	static int k				= 10;
	static int R0k				= 6;
	static int fitIdx;
	static String resultsPath	= "/vega/hblab/users/cr2166/SELEX/chaitanya/scripts/full_shape_sweep/output/";
	static String shapePath		= "/vega/hblab/users/cr2166/SELEX/chaitanya/scripts/ShapeFolder/";
	static String R0ResultsPath	= "/vega/hblab/users/cr2166/SELEX/chaitanya/Round0/";
	static String kmerPath		= "/vega/hblab/users/cr2166/SELEX/chaitanya/scripts/full_shape_sweep/Scr_R1_kmers.txt";
	static boolean[] isNS;
	static ArrayList<double[]> betas;
	
	/**
	 * 
	 * Input param order: 
	 * 		args[0]: RELATIVE path to *.csv file containing new betas
	 * 		args[1]: Index of fit
	 * 		args[2]: RELATIVE path to shape folder
	 * 		args[3]: (optional) RELATIVE path to R0 model
	 * 
	 */
	
	public static void main(String[] args) {
		double sum, Z;
		String currLine;
		ScoringObject so;
		Shape shapeModel;
		MultinomialFit f;
		MultinomialResults baseResults;
		BufferedReader br;
		Round0Model R0;
		long[] tempSeqs;
		ArrayList<Long> seqs = new ArrayList<Long>();
		
		//Read in command line arguments
		if (args.length==3 || args.length==4) {
			resultsPath += args[0];
			fitIdx = Integer.parseInt(args[1]);
			shapePath += args[2];
			// R0 model fixed to Scr model
			R0ResultsPath += (args.length==4) ? "scr_CCACGTC.dat" : args[3];
		} else {
			throw new IllegalArgumentException("Improper number of arguments supplied!");
		}
		//Load R0 model, shape model
		R0 = new Round0Model(R0ResultsPath, l, R0k, true, null, null);
		shapeModel = new Shape(shapePath, true);
		//Bootstrap results for quick manipulation
		baseResults = new MultinomialResults("/vega/hblab/users/cr2166/SELEX/chaitanya/scripts/full_shape_sweep/output/Scr-orthogonal.dat");
		
		//First parse CSV containing file reads	
		try {
			readConfigFile(resultsPath);
		} catch (Exception e1) {
			e1.printStackTrace();
		}
		
		//Next read in all kmers to be evaluated
		try {
			br = new BufferedReader(new FileReader(kmerPath));
			while ((currLine = br.readLine())!=null) {
				seqs.add((new Sequence(currLine, 0, k)).getValue());
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
				
		//Perform inference for a given fit index. First modify fit result
		f = baseResults.getFit(0);
		f.isNSBinding = isNS[fitIdx];
		f.finalPosition = betas.get(fitIdx);
		
		//Load fit and create new scoring object
		so = new ScoringObject(baseResults, 0, -1, false, false, shapeModel, R0);
		Z = so.getZ();
		System.out.println("Z:\t"+Z);
		System.out.println("totCounts:+\t"+baseResults.totCounts);
		//Loop over all kmers and print counts to a file
		for (int i=0; i<seqs.size(); i++) {
			sum = 0;
			tempSeqs = seqBuilder(seqs.get(i));
			for (long currSeq : tempSeqs) {
				sum += R0.getSeqProb(currSeq)*so.seqEval(currSeq);
			}
			System.out.println((new Sequence(seqs.get(i), k)).getString()+"\t"+sum);
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
	
	private static void readConfigFile(String configFile) throws Exception {
		int nFits;
		String currLine;
		double[] tempBetas;
		String[] temp, temp2;
		ArrayList<String> config= new ArrayList<String>(100);
		BufferedReader br		= new BufferedReader(new FileReader(configFile));
		
		//First load all non-commented lines and skip first line
		while ((currLine = br.readLine())!=null) {
			if (currLine.startsWith("//")) {
				continue;
			}
			config.add(currLine);
		}
		br.close();
		
		//Number of fits to compute
		nFits = config.size();
		
		//Create isNS and beta parameters output file
		isNS = new boolean[nFits];
		betas = new ArrayList<double[]>();
		
		//Loop over fits
		for (int currFit=0; currFit<nFits; currFit++) {
			//split line based on comma
			temp = (config.get(currFit)).split(",");
			//Store isNS
			isNS[currFit] = Boolean.parseBoolean(temp[7]);
			//extract betas, including NSB (if included)
			temp2 = Arrays.copyOfRange(temp, 25, temp.length);
			tempBetas = new double[temp2.length];
			for (int i=0; i<temp2.length; i++) {
				tempBetas[i] = Double.parseDouble(temp2[i]);
			}
			betas.add(tempBetas);
		}
	}
}
