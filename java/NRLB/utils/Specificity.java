package utils;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Random;

import base.Array;
import base.Sequence;
import base.Shape;
import main.SimpleKmerCount;
import model.*;

public class Specificity {
	public static String dir		= "/Users/chaitanya/Documents/Research/SELEX/Biological Analysis/";
	public static String R0Dir		= "./Round0/";
	public static String resultsc1	= dir+"/HM Dataset/Fits/Lab.dat";
	public static String resultsc2  = dir+"/HM Dataset/Fits/Scr.dat";
	public static String resultsc3	= dir+"/HM Dataset/Fits/UbxIVa.dat";
	public static int modelIndexc1	= 19;
	public static int modelIndexc2	= 20;	
	public static int modelIndexc3	= 19;
	public static int modeIndexc1	= 0;
	public static int modeIndexc2	= 0;
	public static int modeIndexc3	= 0;
	public static String chrDir		= dir+"/SVB Enhancer/Genomes/dm3/";
	public static String outPath	= "/Users/chaitanya/Desktop/sample2";
	public static double divSize	= 1.0/101.0;
	public static long refSeq		= (new Sequence("AAGATGATTTATGACCTT", 0, 18)).getValue(); //AAGATGATTTATGACCTT
	public static boolean sample	= true;
	public static long nSamples		= 100*1000*1000;
	private static final double s60 = Math.sqrt(3)/2.0;
	
	public static void main(String[] args) {
		int totBases, k, divs, overflow, xbin, ybin;
		long currSeq, stepSize, minBin = 0;
		double ref1, ref2, ref3, a1, a2, a3, total, x, y, max;
		FA chr;
		Round0Model R0;
		MultinomialFit currFit;
		DataOutputStream os, fs;
		ScoringObject c1, c2, c3;
		Random generator		= new Random();
		MultinomialResults mr1	= new MultinomialResults(resultsc1);
		MultinomialResults mr2	= new MultinomialResults(resultsc2);
		MultinomialResults mr3	= new MultinomialResults(resultsc3);
		int[] status			= new int[1];
		String[] parsedPath;
		File[] chrPaths			= (new File(chrDir)).listFiles();
		long[][] binnedOutput;
		double[][] avgAffinity;
		
		//'fix' selected fits
		currFit = mr1.getFit(modelIndexc1);				
		for (int i=0; i<4; i++) {
			currFit.finalPosition[i] = 0;
		}
		currFit.finalPosition = Array.cat(Array.reverse(Arrays.copyOfRange(currFit.finalPosition, 0, 72)), 0);
		currFit = mr3.getFit(modelIndexc3);
		for (int i=0; i<4; i++) {
			currFit.finalPosition[17*4+i] = 0;
		}
		currFit = mr2.getFit(modelIndexc2);
		currFit.finalPosition = Array.cat(Array.reverse(Array.cycleRight(Arrays.copyOfRange(currFit.finalPosition, 0, 72), 4)), 0.0);
		
		//Load Round0 model from results file
		parsedPath	= mr1.r0ModelPath.split("/");
		R0			= new Round0Model(R0Dir+parsedPath[parsedPath.length-1], mr1.l, 
						mr1.r0k, mr1.r0Flank, mr1.lFlank, mr1.rFlank);
		//Create sliding window systems
		c1 = new ScoringObject(mr1, modelIndexc1, modeIndexc1, true, false, null, R0);
		c2 = new ScoringObject(mr2, modelIndexc2, modeIndexc2, true, false, null, R0);
		c3 = new ScoringObject(mr3, modelIndexc3, modeIndexc3, true, false, null, R0);
		if (mr1.getFit(modelIndexc1).isMulti) {
			k = mr1.getFit(modelIndexc1).ks[modeIndexc1];
		} else {
			k = mr1.getFit(modelIndexc1).k;
		}
		//Compute the number of divisions and create output matrix (dependent on l)
		divs 		= (int) Math.ceil(1.0/divSize)+1;
		binnedOutput= new long[divs][divs];
		avgAffinity = new double[divs][divs];
		
		//Compute relative affinity offsets
		ref1 = 1; //c1.seqEval(refSeq); //0.0467658000249295;	//c1.seqEval(refSeq)
		ref2 = 1; //c2.seqEval(refSeq); //0.0458173168909696;	//c2.seqEval(refSeq)
		ref3 = 1; //c3.seqEval(refSeq);
		System.out.println(ref1);
		System.out.println(ref2);
		System.out.println(ref3);
		
		//See if probes should be sampled or if the genome should be used
		if (sample) {
			k = 17;
			stepSize= (long) Math.floor(Math.pow(4, k)/nSamples);
			overflow= (int) (Math.pow(4, k) - ((double) nSamples)*((double) stepSize));
			//Loop over samples
			for (long i=0; i<nSamples; i++) {
//				//Generate next sample
				if (overflow>0) {
					currSeq = (long) (minBin + Math.floor(generator.nextDouble()*(stepSize+1)));
					minBin += stepSize+1;
					overflow--;
				} else {
					currSeq = (long) (minBin + Math.floor(generator.nextDouble()*stepSize));
					minBin += stepSize;
				}
				currSeq <<= 2;
				//Compute relative affinity scores for each of the three classes
				a1 = c1.seqEval(currSeq)/ref1;
				a2 = c2.seqEval(currSeq)/ref2;
				a3 = c3.seqEval(currSeq)/ref3;
				
				//Normalize so they fall on the simplex
				total = a1 + a2 + a3;
				max = Math.max(a1, Math.max(a2, a3));
				a1 /= total;
				a2 /= total;
				//Transform simplex to equilateral triangle
				x  = .5*a1+a2;
				y  = s60*a1;
				//Bin value
				xbin = (int) Math.floor((x/divSize));
				ybin = (int) Math.floor((y/divSize));
//				System.out.println(a1+"\t"+a2+"\t"+a3+"\t"+total+"\t"+xbin+"\t"+ybin);
				binnedOutput[ybin][xbin]++;
				avgAffinity[ybin][xbin] += max;
			}
		} else {
			double tStart = System.currentTimeMillis();
			//Evaluate all sequences in genome
			for (File currChr : chrPaths) {
				if (currChr.isFile()) {
					chr = new FA(chrDir+currChr.getName());
//					if (!chr.name.equals("chrX")) {
//						continue;
//					}
					System.out.println("Currently Processing "+chr.name);
					//Evaluate chromosome
					totBases = chr.nBases;
					for (int start=0; start<totBases-k+1; start++) {
						currSeq = chr.getLong(start, start+k-1, status);
						//evaluate sequence if it does not contain an N, repeats accepted
						if (status[0]!=-1) {
							//Compute relative affinity scores for each of the three classes
							a1 = c1.seqEval(currSeq)/ref1;
							a2 = c2.seqEval(currSeq)/ref2;
							a3 = c3.seqEval(currSeq)/ref3;
							
							//Normalize so they fall on the simplex
							total = a1 + a2 + a3;
							max = Math.max(a1, Math.max(a2, a3));
							a1 /= total;
							a2 /= total;
							//Transform simplex to equilateral triangle
							x  = a2 + .5*a1;
							y  = s60*a1;
							//Bin value
							xbin = (int) Math.floor((x/divSize));
							ybin = (int) Math.floor((y/divSize));
							binnedOutput[ybin][xbin]++;
							avgAffinity[ybin][xbin] += max;
						}
					}
				}
			}
			System.out.println(System.currentTimeMillis()-tStart);	
		}
		System.out.println("Finished scoring sequences.");
		
		//evaluate avgAffinity
		for (int i=0; i<divs; i++) {
			for (int j=0; j<divs; j++) {
				avgAffinity[i][j] /= binnedOutput[i][j];
			}
		}

		//Print binned output to file
		try {
			PrintStream outputFile	= new PrintStream(new FileOutputStream(outPath+"_count.tsv"));
			PrintStream original	= System.out;
			
			System.setOut(outputFile);
			for (int i=divs-1; i>=0; i--) {
				Array.print(binnedOutput[i]);
			}
			System.setOut(original);
			outputFile = new PrintStream(new FileOutputStream(outPath+"_aff.tsv"));
			System.setOut(outputFile);
			for (int i=divs-1; i>=0; i--) {
				Array.print(avgAffinity[i]);
			}
			System.setOut(original);
		} catch (FileNotFoundException e) {
			System.out.println("Cannot create trajectory file at this "
					+ "location: "+outPath);
			e.printStackTrace();
		}
		System.out.println("Finished printing tables.");
	}
}
