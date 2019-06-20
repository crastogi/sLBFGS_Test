package utils;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

import base.Array;
import base.Sequence;
import base.Shape;
import model.Data;
import model.MultinomialFit;
import model.MultinomialResults;
import model.Round0Model;

public class ModelProfiler {
	static boolean kmerAnalysis = true;
	static boolean modeAnalysis = false;
	static boolean positionBias = false;
	static boolean divergence	= false;
	static boolean singleRound	= true;
	static int l				= 16;
	static long nSamples		= 500*1000*1000;
	static String R0Results		= "./Round0/null.dat";
	static int R0k				= 1;
	static String resultsPath	= "C:/Users/Chaitanya/Desktop/HT.dat";
	static int modelIndex		= 0;
	static String dataPath		= "F:/Data/Datasets/tmp/pbflExd.r1.v16.Pbfl-Exd.16mer1.1.16.dat";
	static double lowerBound	= -7;
	static double upperBound	= -0;
	static double divSize		= .001;
	static String outputPath	= "C:/Users/Chaitanya/Desktop/Scr_HTSELEX.tsv";
	static String kmerPath		= "C:/Users/Chaitanya/Desktop/Mode Profiles/KmerTables/Scr_R1_kmers.txt";
	static String shapePath		= "C:/Users/Chaitanya/Desktop/ShapeFolder/Orthogonal/";
	static int k				= 10;
	
	public static void main(String[] args) {
		int nModes, divs, maxCount = 0, totCount = 0, currBin, currRatio, overflow;
		long stepSize, minBin = 0, nAvailable, currProbe;
		double nsBinding = 0, Z, affSum, score, total;
		int[] counts;
		long[] probes;
		double[] aff;
		double[][] output = null;
		double[][] empDist;
		Data data;
		Shape shapeModel;
		MultinomialFit fit;
		MultinomialResults r= new MultinomialResults(resultsPath);
		Random generator	= new Random();
		FileUtilities reader= new FileUtilities(l);
		Round0Model R0 		= new Round0Model(R0Results, l, R0k, true, null, null);
		ScoringObject[] so;
		
		//Load model
		fit	= r.getFit(modelIndex);
		if (fit.isMulti) {
			nModes = fit.ks.length;
		} else {
			nModes = 1;
		}
		if (fit.isNSBinding) {
			nsBinding = Math.exp(fit.finalPosition[fit.finalPosition.length-1]);
		}
		if (fit.isShape) {
			shapeModel = new Shape(shapePath, true);
		} else {
			shapeModel = null;
		}
		//Compute the number of divisions
		divs = (int) Math.ceil((upperBound-lowerBound)/divSize)+1;
		
		//Run mode analysis for only one round?
		if (!kmerAnalysis) {
		if (modeAnalysis & singleRound) {
			//Create scoring objects
			so = new ScoringObject[nModes];
			Z  = 0;
			for (int m=0; m<nModes; m++) {
				so[m] = new ScoringObject(r, modelIndex, m, false, false, shapeModel, R0);
				Z += so[m].getZ();
			}
			Z += nsBinding;
			//Read in probes
			data	= reader.readSeqFile(dataPath, false, 0, 0, null);
			System.out.println("Total number of probes in file: "+data.nCount);
			counts	= data.counts;
			probes	= data.probes;
			
			//Create output matrix
			output		= new double[divs][nModes+1];
			aff			= new double[nModes+1];
			aff[nModes]	= nsBinding;
			
			//Score all probes within the dataset
			for (int i=0; i<probes.length; i++) {
				//Score each mode
				for (int m=0; m<nModes; m++) {
					aff[m] = so[m].seqEval(probes[i])-nsBinding;
				}
				affSum = Array.sum(aff);
				currBin = (int) (Math.floor((Math.log10(R0.getSeqProb(probes[i])*affSum/Z)-lowerBound)/divSize)+1);

				for (int j=0; j<=nModes; j++) {
					output[currBin][j] += counts[i]*aff[j]/affSum;
				}
			}
			System.out.println("Finished scoring dataset.");
		} else {
			//Create scoring object
			so = new ScoringObject[1];
			so[0] = new ScoringObject(r, modelIndex, -1, false, false, shapeModel, R0);
			Z = so[0].getZ();
			//Read in probes
			data	= reader.readSeqFile(dataPath, false, 0, 0, null);
			System.out.println("Total number of probes in file: "+data.nCount);
			counts	= data.counts;
			probes	= data.probes;
			
			if (divergence) {
				//Create output matrix
				output		= new double[divs][2];
				
				if (l>16) {
					//Score all probes within the dataset
					for (int i=0; i<probes.length; i++) {
						currBin = (int) (Math.floor((Math.log10(R0.getSeqProb(probes[i])*so[0].seqEval(probes[i])/Z)-lowerBound)/divSize)+1);
						output[currBin][1] += counts[i];
					}
					System.out.println("Finished scoring dataset.");
					//Now estimate universe
					stepSize= (long) Math.floor(Math.pow(4, l)/nSamples);
					overflow= (int) (Math.pow(4, l) - ((double) nSamples)*((double) stepSize));
					total	= 0;
					for (long i=0; i<nSamples; i++) {
						if (overflow>0) {
							currProbe = (long) (minBin + Math.floor(generator.nextDouble()*(stepSize+1)));
							minBin += stepSize+1;
							overflow--;
						} else {
							currProbe = (long) (minBin + Math.floor(generator.nextDouble()*stepSize));
							minBin += stepSize;
						}
						//Score each mode
						score = Math.log10(R0.getSeqProb(currProbe)*so[0].seqEval(currProbe)/Z);
						currBin = (int) (Math.floor((score-lowerBound)/divSize)+1);
						output[currBin][0] += score;
						total += score;
					}
					//Normalize
					for (int i=0; i<divs; i++) {
						output[i][0] /= total;
					}
					System.out.println("Finished sampling universe.");
				} else {
					//First score all probes of length l
					for (long i=0; i<Math.pow(4, R0.getL()); i++) {
						score = Math.log10(R0.getSeqProb(i)*so[0].seqEval(i)/Z);
						currBin = (int) (Math.floor((score-lowerBound)/divSize)+1);
						output[currBin][0] += score;
					}
					System.out.println("Finished scoring entire universe.");
					//Now score all probes within the dataset
					for (int i=0; i<probes.length; i++) {
						currBin = (int) (Math.floor((Math.log10(R0.getSeqProb(probes[i])*so[0].seqEval(probes[i])/Z)-lowerBound)/divSize)+1);
						output[currBin][1] += counts[i];
					}
					System.out.println("Finished scoring dataset.");
				}
			} else if (positionBias) {
				//Create output matrix
				output		= new double[divs][101];
				
				//Score all probes within the dataset
				for (int i=0; i<probes.length; i++) {
					//Score each mode
					aff		 = so[0].perWindowSeqEval(probes[i]);
					affSum	 = Array.sum(aff);
					currBin  = (int) (Math.floor((Math.log10((affSum+nsBinding))-lowerBound)/divSize)+1); //R0.getSeqProb(probes[i])*(affSum+nsBinding)/Z
					currRatio= (int) (Math.floor((Array.max(aff)/affSum)/0.01)+1);
					output[currBin][currRatio]++;
				}
				System.out.println("Finished scoring dataset.");
			} else {
				//Find the highest count value
				for (int count : counts) {
					totCount += count;
					if (count > maxCount) {
						maxCount = count;
					}
				}
				//Create output matrix
				output		= new double[divs][maxCount+1];
				
				if (l>16) {
					//Score all probes within the dataset
					for (int i=0; i<probes.length; i++) {
						//Score each mode
						currBin = (int) (Math.floor((Math.log10(R0.getSeqProb(probes[i])*so[0].seqEval(probes[i])/Z)-lowerBound)/divSize)+1);
						output[currBin][counts[i]]++;
					}
					System.out.println("Finished scoring dataset.");
					//Now estimate universe
					empDist	= new double[divs][1];
					stepSize= (long) Math.floor(Math.pow(4, l)/nSamples);
					overflow= (int) (Math.pow(4, l) - ((double) nSamples)*((double) stepSize));
					for (long i=0; i<nSamples; i++) {
						if (overflow>0) {
							currProbe = (long) (minBin + Math.floor(generator.nextDouble()*(stepSize+1)));
							minBin += stepSize+1;
							overflow--;
						} else {
							currProbe = (long) (minBin + Math.floor(generator.nextDouble()*stepSize));
							minBin += stepSize;
						}
						//Score each mode
						currBin = (int) (Math.floor((Math.log10(R0.getSeqProb(currProbe)*so[0].seqEval(currProbe)/Z)-lowerBound)/divSize)+1);
						empDist[currBin][0]++;
					}
					//Create frequency distribution based on samples
					nAvailable = (long) Math.pow(4, l) - totCount;
					for (int i=0; i<divs; i++) {
						output[i][0] = (long) (empDist[i][0]/nSamples*nAvailable);
					}
					System.out.println("Finished sampling universe.");
				} else {
					//First score all probes of length l
					for (long i=0; i<Math.pow(4, R0.getL()); i++) {
						//Score each mode
						currBin = (int) (Math.floor((Math.log10(R0.getSeqProb(i)*so[0].seqEval(i)/Z)-lowerBound)/divSize)+1);
						output[currBin][0]++;
					}
					System.out.println("Finished scoring entire universe.");
					//Now score all probes within the dataset
					for (int i=0; i<probes.length; i++) {
						//Score each mode
						currBin = (int) (Math.floor((Math.log10(R0.getSeqProb(probes[i])*so[0].seqEval(probes[i])/Z)-lowerBound)/divSize)+1);
						output[currBin][counts[i]]++;
						output[currBin][0]--;
					}
					System.out.println("Finished scoring dataset.");
				}
			}
		}
		
		//Lastly, print the table to an output file
		try {
			PrintStream outputFile	= new PrintStream(new FileOutputStream(outputPath));
			PrintStream original	= System.out;
			
			System.setOut(outputFile);
			for (int i=0; i<divs; i++) {
				Array.print(output[i]);
			}
			System.setOut(original);
		} catch (FileNotFoundException e) {
			System.out.println("Cannot create trajectory file at this "
					+ "location: "+outputPath);
			e.printStackTrace();
		}
		System.out.println("Finished printing tables.");
		} else {
			//Kmer analysis
			double sum;
			String currLine;
			BufferedReader br;
			long[] tempSeqs;
			ArrayList<Long> seqs = new ArrayList<Long>();
			PrintStream outputFile, original = System.out;
			
			//Create scoring object
			so = new ScoringObject[1];
			so[0] = new ScoringObject(r, modelIndex, -1, false, false, shapeModel, R0);
			Z = so[0].getZ();
			System.out.println(R0.getZ());
			System.out.println(Z);
						
			//First read in all kmers
			try {
				outputFile = new PrintStream(new FileOutputStream(outputPath));
				br = new BufferedReader(new FileReader(kmerPath));
				while ((currLine = br.readLine())!=null) {
					seqs.add((new Sequence(currLine, 0, k)).getValue());
				}
				br.close();
				System.out.println("Finished loading kmers.");
				
				//Loop over all probes
				System.setOut(outputFile);
				for (int i=0; i<seqs.size(); i++) {
					sum = 0;
					tempSeqs = seqBuilder(seqs.get(i));
					for (long currSeq : tempSeqs) {
						sum += R0.getSeqProb(currSeq)*so[0].seqEval(currSeq)/Z;
					}
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
