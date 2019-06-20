package utils;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;

import base.Array;
import base.Sequence;
import base.Shape;
import model.Data;
import model.MultinomialResults;
import model.Round0Model;
import model.SingleModeSW;

public class ShapeProfiler {
	public static String dir		= "/Users/chaitanya/Documents/Research/SELEX/Biological Analysis/";
	public static String R0Dir		= "./Round0/";
	public static String shapeDir	= "./../Shape/";
	public static String results	= dir+"/HM Dataset/Fits/Scr_Monomer.dat";
	public static int modelIndex	= 23;
	public static int modeIndex		= -1;
	public static String[] shape	= new String[]{"MGW"};
	public static String dataFile	= "/Users/chaitanya/Documents/Magic Briefcase/exdLab.exdAbdb.L.1.barcodeCCAGCTG.v2.low.1.16.dat";
	public static int l				= 16;
	public static boolean isFilter	= true;
	public static int maxBaseCount	= 12;
	public static int maxBaseRepeat	= 12;
	public static String outDir		= dir+"/Shape/Scr_Monomer";
	public static double divSize	= 0.1;
	public static double lowDivStep	= 1.5;
	public static int lowDivs		= 30;
	public static boolean isLog		= false;
	public static boolean useDataset= false;
	public static int minBinCount	= 100;
	public static int maxBinCount	= 1000;
	public static String bestMotif	= "TTTAATTACTT";
	public static double normalizer	= 0.4164444137896402;
	
	@SuppressWarnings("unchecked")
	public static void main(String[] args) {
		boolean unfinished				= true;
		int k, nShapeClasses, nDivs, currBin, nextOffset, counter=0, currLim=4;
		long currSeq, bestSeq, mask, addMask;
		double lowRange, highRange, lowBinOffset;
		Data dataset;
		Shape shapeModel, evalShapeModel;
		Random generator				= new Random();
		Round0Model R0;
		SingleModeSW sw;
		ScoringObject so;
		FileUtilities reader			= new FileUtilities(l);
		MultinomialResults mr			= new MultinomialResults(results);
		ArrayList<Long> tempSeq;
		int[] position, binCounts;
		long[] testSeqs;
		double[] weight, binWeight;
		String[] parsedPath;
		Object[] output;
		HashSet<Long>[] probeSet, lowProbeSet;
		ArrayList<Long>[] probes;
		double[][] tempShapes;
		double[][][] shapes, binShapes;
		PrintStream outputFile;
		
		
		//Load Round0 model from results file
		parsedPath	= mr.r0ModelPath.split("/");
		R0			= new Round0Model(R0Dir+parsedPath[parsedPath.length-1], 
						mr.l, mr.r0k, mr.r0Flank, mr.lFlank, mr.rFlank);
		//Create new sliding window system
		if (mr.getFit(modelIndex).isShape) {
			shapeModel = new Shape(shapeDir, true);
			shapeModel.setFeatures(mr.getFit(modelIndex).shapes);
		} else {
			shapeModel = null;
		}
		evalShapeModel = new Shape(shapeDir, true);
		evalShapeModel.setFeatures(shape);
		nShapeClasses  = evalShapeModel.nShapeFeatures();
		
		if (useDataset) {
			so = new ScoringObject(mr, modelIndex, modeIndex, shapeModel, R0, 
					evalShapeModel);
			k  = so.getK();

			//Load data files
			dataset = reader.readSeqFile(dataFile, isFilter, maxBaseCount, 
					maxBaseRepeat, null);
			
			//Compute profiles
			System.out.println("Computing profiles.");
			output	= so.shapeEval(dataset.probes);
			weight	= (double[]) output[0];
			position= (int[]) output[1];
			shapes	= (double[][][]) output[2];
			System.out.println("Number of elements: "+weight.length+".");
			
			//Transform weights to log scale if necessary
			if (isLog) {
				for (int i=0; i<weight.length; i++) {
					weight[i] = Math.log(weight[i]);
				}
			}
			
			//Compute Affinity Range and Divisions
			if (isLog) {
				lowRange	= Double.POSITIVE_INFINITY;
				highRange	= Double.NEGATIVE_INFINITY;
				for (int i=0; i<weight.length; i++) {
					if (weight[i] < lowRange) {
						lowRange = weight[i];
					} else if (weight[i] > highRange) {
						highRange = weight[i];
					}
				}
			} else {
				highRange	= 1;
				lowRange	= 0;
			}
			nDivs	= (int) Math.ceil((highRange-lowRange)/divSize)+1;
			
			//Bin, clean profiles, average profiles, and print data to file
			System.out.println("Binning.");
			binCounts = new int[nDivs];
			binShapes = new double[nDivs][nShapeClasses][k];
			//Bin shapes
			for (int i=0; i<weight.length; i++) {
				currBin = (int) (Math.floor((weight[i]-lowRange)/divSize));
				binCounts[currBin]++;
				for (int j=0; j<nShapeClasses; j++) {
					for (int m=0; m<k; m++) {
						binShapes[currBin][j][m] += shapes[i][j][position[i]+m];
					}
				}
			}
			//Merge last bin
			nDivs--;
			binCounts[nDivs-1] += binCounts[nDivs];
			for (int j=0; j<nShapeClasses; j++) {
				for (int m=0; m<k; m++) {
					binShapes[nDivs-1][j][m] += binShapes[nDivs][j][m];
				}
			}
			//Average shapes
			for (int i=0; i<nDivs; i++) {
				if (binCounts[i]==0)	continue;
				for (int j=0; j<nShapeClasses; j++) {
					for (int m=0; m<k; m++) {
						binShapes[i][j][m] /= binCounts[i];
					}
				}
			}
			//Assign bin weights
			binWeight = new double[nDivs];
			for (int i=0; i<nDivs; i++) {
				binWeight[i] = lowRange+divSize*(i+0.5);
			}
		} else {		
			so			= new ScoringObject(mr, modelIndex, modeIndex, true, 
							false, shapeModel, R0);
			k			= so.getK();
			sw			= new SingleModeSW(evalShapeModel, R0, k, k, false, 0, 
							"AA", "AA", false, false);
			lowRange	= 0;
			nDivs		= (int) Math.ceil(1.0/divSize)+1;
			lowBinOffset= -Math.log(divSize)/Math.log(lowDivStep)+lowDivs;
			
			//Create unique list of binned sequences
			System.out.println("Generating sequences.");
			mask		= (long) Math.pow(4, k) - 1;
			addMask 	= currLim - 1;
			binCounts	= new int[nDivs+lowDivs];
			bestSeq		= (new Sequence(bestMotif, 0, k)).getValue();
			currSeq		= bestSeq;
			probeSet	= (HashSet<Long>[]) new HashSet[nDivs];
			lowProbeSet	= (HashSet<Long>[]) new HashSet[lowDivs];
			for (int i=0; i<nDivs; i++) {
				probeSet[i] = new HashSet<Long>();
			}
			for (int i=0; i<lowDivs; i++) {
				lowProbeSet[i]	= new HashSet<Long>();
			}
			
			//Generate random sequences until bin conditions are satisfied
			testSeqs	= mutator(bestSeq, k);
			//Dispatch & Bin
			weight		= so.seqEval(testSeqs);
			weight		= Array.scalarMultiply(weight, 1/normalizer);
			for (int i=0; i<testSeqs.length; i++) {
				if (weight[i]<=divSize) {
					currBin = (int) (Math.floor(Math.log(weight[i])/Math.log(lowDivStep)+lowBinOffset));
					if (currBin<0)	currBin = 0;
					lowProbeSet[currBin].add(testSeqs[i]);
				} else {
					currBin	= (int) (Math.floor((weight[i])/divSize));
					probeSet[currBin].add(testSeqs[i]);	
				}
			}
			for (int i=0; i<lowDivs; i++) {
				binCounts[i] = lowProbeSet[i].size();
			}
			for (int i=0; i<nDivs; i++) {
				binCounts[lowDivs+i] = probeSet[i].size();
			}
			while (unfinished) {
				//Find a random sequence within bin and mutate it
				for (int i=0; i<lowDivs; i++) {
					if (binCounts[i]>maxBinCount || binCounts[i]==0) continue;
					tempSeq = new ArrayList<Long>(lowProbeSet[i]);
					testSeqs= mutator(tempSeq.get(generator.nextInt(tempSeq.size())), k);
					//Dispatch & Bin
					weight	= so.seqEval(testSeqs);
					weight	= Array.scalarMultiply(weight, 1/normalizer);
					if (weight[i]<=divSize) {
						currBin = (int) (Math.floor(Math.log(weight[i])/Math.log(lowDivStep)+lowBinOffset));
						if (currBin<0)	currBin = 0;
						lowProbeSet[currBin].add(testSeqs[i]);
					} else {
						currBin	= (int) (Math.floor((weight[i])/divSize));
						probeSet[currBin].add(testSeqs[i]);	
					}
				}
				for (int i=1; i<nDivs; i++) {
					if (binCounts[lowDivs+i]>maxBinCount || binCounts[lowDivs+i]==0) continue;
					tempSeq	= new ArrayList<Long>(probeSet[i]);
					testSeqs= mutator(tempSeq.get(generator.nextInt(tempSeq.size())), k);
					//Dispatch & Bin
					weight	= so.seqEval(testSeqs);
					weight	= Array.scalarMultiply(weight, 1/normalizer);
					if (weight[i]<=divSize) {
						currBin = (int) (Math.floor(Math.log(weight[i])/Math.log(lowDivStep)+lowBinOffset));
						if (currBin<0)	currBin = 0;
						lowProbeSet[currBin].add(testSeqs[i]);
					} else {
						currBin	= (int) (Math.floor((weight[i])/divSize));
						probeSet[currBin].add(testSeqs[i]);	
					}
				}
				
				//Generate sequences with random mutations away from optimal
				testSeqs = new long[100000];
				for (int i=0; i<100000; i++) {
					nextOffset = 2*generator.nextInt(k);
					if (generator.nextDouble()<.1) {
						currSeq = bestSeq;
					}
					currSeq &= ~(addMask << nextOffset);
					currSeq	= (currSeq | (((long) generator.nextInt(currLim)) << nextOffset)) & mask;
					testSeqs[i] = currSeq;
				}
				//Dispatch & Bin
				weight = so.seqEval(testSeqs);
				weight = Array.scalarMultiply(weight, 1/normalizer);
				for (int i=0; i<100000; i++) {
					if (weight[i]<=divSize) {
						currBin = (int) (Math.floor(Math.log(weight[i])/Math.log(lowDivStep)+lowBinOffset));
						if (currBin<0)	currBin = 0;
						lowProbeSet[currBin].add(testSeqs[i]);
					} else {
						currBin	= (int) (Math.floor((weight[i])/divSize));
						probeSet[currBin].add(testSeqs[i]);	
					}				
				}
				
				//Test counts for termination
				unfinished = false;
				for (int i=0; i<lowDivs; i++) {
					binCounts[i] = lowProbeSet[i].size();
					if (binCounts[i]>1E6) {
						binCounts[i]=0;
						lowProbeSet[i] = new HashSet<Long>();
					}
				}
				for (int i=0; i<nDivs; i++) {
					binCounts[lowDivs+i] = probeSet[i].size();
					if (binCounts[lowDivs+i]>1E6) {
						binCounts[lowDivs+i]=0;
						probeSet[lowDivs+i] = new HashSet<Long>();
					}
				}
				for (int i=0; i<lowDivs+nDivs; i++) {
					if (binCounts[i]<maxBinCount) {
						unfinished = true;
					}
				}				
				//Control the size of the random sequence being added
				counter++;
				if (counter==10) {
					currLim = (currLim)*4;
					addMask= currLim - 1;
					counter = 0;
					if (currLim > 2E9 || currLim <= 0) {
						break;
					}
				}
			}
						
			//Convert HashSets to ArrayLists and merge last bin
			nDivs		= nDivs + lowDivs - 2;
			probes		= (ArrayList<Long>[]) new ArrayList[nDivs];
			binWeight	= new double[nDivs];
			binCounts	= new int[nDivs];
			probeSet[probeSet.length-2].addAll(probeSet[probeSet.length-1]);
			for (int i=0; i<lowDivs; i++) {
				probes[i]	= new ArrayList<Long>(lowProbeSet[i]);
				binWeight[i]= Math.pow(lowDivStep, -lowBinOffset+i+.5);
			}
			for (int i=lowDivs+1; i<nDivs+1; i++) {
				probes[i-1]		= new ArrayList<Long>(probeSet[i-lowDivs]);
				binWeight[i-1]	= divSize*(i-lowDivs+.5); 
			}
			for (int i=0; i<nDivs; i++) {
				binCounts[i]	= probes[i].size();
			}
			
			//Reduce bins to appropriate size by randomly selecting probes
			for (int i=0; i<nDivs; i++) {
				if (binCounts[i]>maxBinCount) {
					tempSeq = new ArrayList<Long>();
					for (int j=0; j<maxBinCount; j++){
						currBin = generator.nextInt(binCounts[i]);
						tempSeq.add(probes[i].get(currBin));
						probes[i].remove(currBin);
						binCounts[i]--;
					}
					probes[i] = tempSeq;
				}
			}
			System.out.println("Final probe counts per bin: [lowest to highest]");
			binCounts = new int[nDivs];
			for (int i=0; i<nDivs; i++) {
				binCounts[i] = probes[i].size();
				System.out.println(binCounts[i]);
			}
			
			//Create shape profiles; start with a new sliding window object
			binShapes = new double[nDivs][nShapeClasses][k-4];
			for (int i=0; i<nDivs; i++) {
				if (binCounts[i]==0)	continue;
				for (int j=0; j<binCounts[i]; j++) {
					tempShapes = sw.swShapeProfile(probes[i].get(j));
					for (int o=0; o<nShapeClasses; o++) {
						for (int p=2; p<k-2; p++) {
							binShapes[i][o][p-2] += tempShapes[o][p];
						}
					}
				}
				//Average shapes
				for (int j=0; j<nShapeClasses; j++) {
					for (int m=2; m<k-2; m++) {
						binShapes[i][j][m-2] /= binCounts[i];
					}
				}
			}
		}
		
		//Save profiles
		System.out.println("Printing profiles.");
		for (int j=0; j<nShapeClasses; j++) {
			try {
				outputFile = new PrintStream(new FileOutputStream(outDir+"_shape"+j+".txt"));
				System.setOut(outputFile);	
				for (int i=0; i<nDivs; i++) {
					System.out.print(binWeight[i]+"\t"+binCounts[i]+"\t");
					Array.print(binShapes[i][j]);
				}
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
	}
	
	private static long[] mutator(long seedMotif, int k) {
		long currSeq;
		long[] output = new long[k*4];
		
		for (int i=0; i<k; i++) {
			currSeq = seedMotif;
			currSeq &= ~((long) 3 << 2*i);
			for (int j=0; j<4; j++) {
				output[4*i+j] = (currSeq | ((long) j << 2*i));
			}
		}
		
		return output;
	}
}
