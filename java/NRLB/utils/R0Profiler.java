package utils;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.PrintStream;
import java.util.Random;

import base.Array;
import base.ArraysMergerHeap;
import base.CountObject;
import base.MarkovModelInfo;
import base.MarkovModelOption;
import base.Sequence;
import config.ExperimentReference;
import main.SELEX;
import main.SimpleKmerCount;
import model.Data;
import model.Round0Model;

public class R0Profiler {
	static int l				= 26;
	static boolean sampleTest	= true;
	static long nSamples		= 500*1000*1000;
	static String R0Results		= "./Round0/p53_unmeth.dat";
	static String R0DataPath	= "F:/Data/Datasets/tmp/p53-unmeth-R0.R0.0.26.dat";
	static int R0k				= 4;
	static boolean isLog		= true;
	static double lowerBound	= -20;
	static double upperBound	= -7.3;
	static double divSize		= .001;
	static String outputPath	= "F:/Data/ATF_Test.tsv";
	static String workingDir	= "F:/Data/Datasets/tmp/";
	static boolean testMarkov	= false;
	static String r0ConfigFile	= "./config/Scr-R0.xml";
	static String sequencingRun	= "exdUbx.exdScr.0";
	static String testSample	= "barcodeCCACGTC.v1";
	static String trainSample	= "barcodeCCAGCTG.v1";
	static int modelLength		= 6;
	static String modelMethod	= MarkovModelOption.DIVISION;
	static float[] modelProbs 	= null;
	static int[] modelCounts 	= null;
	static long modelTotalCount	= 0;
	static String modelObjPath;
	
	public static void main(String[] args) {
		int maxCount = 0, totCount = 0, divs, currBin, countOffset;
		int overflow;
		long stepSize, currProbe, minBin = 0, nAvailable;
		double probeVal;
		Data R0Data;
		ExperimentReference trainingSet, testingSet;
		FileUtilities reader= new FileUtilities(16);
		Round0Model R0 		= new Round0Model(R0Results, l, R0k, true, null, null);
		Random generator	= new Random();
		int[] counts;
		long[] probes;
		double[] empDist;
		long[][] output;
		
		if (testMarkov) {
			SELEX.setWorkingDirectory(workingDir);
			SELEX.loadConfigFile(r0ConfigFile);
			trainingSet	= SELEX.getExperimentReference(sequencingRun, trainSample, 0);
			testingSet	= SELEX.getExperimentReference(sequencingRun, testSample, 0);
			loadMarkovModel(trainingSet, testingSet, modelLength, modelMethod);
		}

		//Read in R0 probes
		R0Data = reader.readSeqFile(R0DataPath, false, 0, 0, null);
		System.out.println("Total number of probes in file: "+R0Data.nCount);
		counts = R0Data.counts;
		probes = R0Data.probes;
		
		//Find the highest count value
		for (int count : counts) {
			totCount += count;
			if (count > maxCount) {
				maxCount = count;
			}
		}
		//Compute the number of divisions and create output matrix (dependent on l)
		divs	= (int) Math.ceil((upperBound-lowerBound)/divSize)+1;
		//Complete enumeration is an issue for l>16
		if (l>16) {
			if (sampleTest) {
				output		= new long[divs][maxCount+1];
				countOffset	= 0;
			} else {
				output		= new long[divs][maxCount];
				countOffset	= -1;
			}
//			//Score all probes within the dataset
			if (isLog) {
				if (testMarkov) {
					for (int i=0; i<probes.length; i++) {
						probeVal= SimpleKmerCount.getPredictedCount(modelLength, modelTotalCount, 
								new Sequence(i, l), modelCounts, modelProbs);
						currBin = (int) (Math.floor((Math.log10(probeVal)-lowerBound)/divSize)+1);
						output[currBin][counts[i]+countOffset]++;
					}
				} else {
					for (int i=0; i<probes.length; i++) {
						probeVal= R0.getSeqProb(probes[i]);
						currBin = (int) (Math.floor((Math.log10(probeVal)-lowerBound)/divSize)+1);
						output[currBin][counts[i]+countOffset]++;
					}
				}
			} else {
				if (testMarkov) {
					for (int i=0; i<probes.length; i++) {
						probeVal= SimpleKmerCount.getPredictedCount(modelLength, modelTotalCount, 
								new Sequence(i, l), modelCounts, modelProbs);
						currBin = (int) (Math.floor((probeVal-lowerBound)/divSize)+1);
						output[currBin][counts[i]+countOffset]++;
					}
				} else {
					for (int i=0; i<probes.length; i++) {
						probeVal= R0.getSeqProb(probes[i]);
						currBin = (int) (Math.floor((probeVal-lowerBound)/divSize)+1);
						output[currBin][counts[i]+countOffset]++;
					}
				}
			}
			System.out.println("Finished scoring dataset.");
			if (sampleTest) {
				empDist	= new double[divs];
				stepSize= (long) Math.floor(Math.pow(4, l)/nSamples);
				overflow= (int) (Math.pow(4, l) - ((double) nSamples)*((double) stepSize));
				if (isLog) {
					if (testMarkov) {
						for (long i=0; i<nSamples; i++) {
							if (overflow>0) {
								currProbe = (long) (minBin + Math.floor(generator.nextDouble()*(stepSize+1)));
								minBin += stepSize+1;
								overflow--;
							} else {
								currProbe = (long) (minBin + Math.floor(generator.nextDouble()*stepSize));
								minBin += stepSize;
							}
							probeVal= SimpleKmerCount.getPredictedCount(modelLength, modelTotalCount, 
									new Sequence(currProbe, l), modelCounts, modelProbs);
							currBin = (int) (Math.floor((Math.log10(probeVal)-lowerBound)/divSize)+1);
							empDist[currBin]++;
						}
					} else {
						for (long i=0; i<nSamples; i++) {
							if (overflow>0) {
								currProbe = (long) (minBin + Math.floor(generator.nextDouble()*(stepSize+1)));
								minBin += stepSize+1;
								overflow--;
							} else {
								currProbe = (long) (minBin + Math.floor(generator.nextDouble()*stepSize));
								minBin += stepSize;
							}
							probeVal= R0.getSeqProb(currProbe);
							currBin = (int) (Math.floor((Math.log10(probeVal)-lowerBound)/divSize)+1);
							empDist[currBin]++;
						}
					}
				} else {
					if (testMarkov) {
						for (long i=0; i<nSamples; i++) {
							if (overflow>0) {
								currProbe = (long) (minBin + Math.floor(generator.nextDouble()*(stepSize+1)));
								minBin += stepSize+1;
								overflow--;
							} else {
								currProbe = (long) (minBin + Math.floor(generator.nextDouble()*stepSize));
								minBin += stepSize;
							}
							probeVal= SimpleKmerCount.getPredictedCount(modelLength, modelTotalCount, 
									new Sequence(currProbe, l), modelCounts, modelProbs);
							currBin = (int) (Math.floor((probeVal-lowerBound)/divSize)+1);
							empDist[currBin]++;
						}
					} else {
						for (long i=0; i<nSamples; i++) {
							if (overflow>0) {
								currProbe = (long) (minBin + Math.floor(generator.nextDouble()*(stepSize+1)));
								minBin += stepSize+1;
								overflow--;
							} else {
								currProbe = (long) (minBin + Math.floor(generator.nextDouble()*stepSize));
								minBin += stepSize;
							}
							probeVal= R0.getSeqProb(currProbe);
							currBin = (int) (Math.floor((probeVal-lowerBound)/divSize)+1);
							empDist[currBin]++;
						}
					}
				}
				//Create frequency distribution based on samples
				nAvailable = (long) Math.pow(4, l) - totCount;
				for (int i=0; i<divs; i++) {
					output[i][0] = (long) (empDist[i]/nSamples*nAvailable);
				}
				System.out.println("Finished sampling universe.");
			}
		} else {
			output	= new long[divs][maxCount+1];
			//First score all probes of length l
			if (isLog) {
				if (testMarkov) {
					for (long i=0; i<Math.pow(4, R0.getL()); i++) {
						probeVal= SimpleKmerCount.getPredictedCount(modelLength, modelTotalCount, 
								new Sequence(i, l), modelCounts, modelProbs);
						currBin = (int) (Math.floor((Math.log10(probeVal)-lowerBound)/divSize)+1);
						output[currBin][0]++;
					}
				} else {
					for (long i=0; i<Math.pow(4, R0.getL()); i++) {
						probeVal= R0.getSeqProb(i);
						currBin	= (int) (Math.floor((Math.log10(probeVal)-lowerBound)/divSize)+1);
						output[currBin][0]++;
					}	
				}					
			} else {
				if (testMarkov) {
					for (long i=0; i<Math.pow(4, R0.getL()); i++) {
						probeVal= SimpleKmerCount.getPredictedCount(modelLength, modelTotalCount, 
								new Sequence(i, l), modelCounts, modelProbs);
						currBin = (int) (Math.floor((probeVal-lowerBound)/divSize)+1);
						output[currBin][0]++;
					}
				} else {
					for (long i=0; i<Math.pow(4, R0.getL()); i++) {
						probeVal= R0.getSeqProb(i);
						currBin	= (int) (Math.floor((probeVal-lowerBound)/divSize)+1);
						output[currBin][0]++;
					}
				}
			}
			System.out.println("Finished scoring entire universe.");
			//Next score all probes within the dataset
			if (isLog) {
				if (testMarkov) {
					for (int i=0; i<probes.length; i++) {
						probeVal= SimpleKmerCount.getPredictedCount(modelLength, modelTotalCount, 
								new Sequence(i, l), modelCounts, modelProbs);
						currBin = (int) (Math.floor((Math.log10(probeVal)-lowerBound)/divSize)+1);
						output[currBin][0]--;
						output[currBin][counts[i]]++;
					}
				} else {
					for (int i=0; i<probes.length; i++) {
						probeVal= R0.getSeqProb(probes[i]);
						currBin = (int) (Math.floor((Math.log10(probeVal)-lowerBound)/divSize)+1);
						output[currBin][0]--;
						output[currBin][counts[i]]++;
					}
				}
			} else {
				if (testMarkov) {
					for (int i=0; i<probes.length; i++) {
						probeVal= SimpleKmerCount.getPredictedCount(modelLength, modelTotalCount, 
								new Sequence(i, l), modelCounts, modelProbs);
						currBin = (int) (Math.floor((probeVal-lowerBound)/divSize)+1);
						output[currBin][0]--;
						output[currBin][counts[i]]++;
					}
				} else {
					for (int i=0; i<probes.length; i++) {
						probeVal= R0.getSeqProb(probes[i]);
						currBin = (int) (Math.floor((probeVal-lowerBound)/divSize)+1);
						output[currBin][0]--;
						output[currBin][counts[i]]++;
					}
				}
			}	
			System.out.println("Finished scoring dataset.");
		}
		//Lastly, print the table to an output file
		try {
			PrintStream outputFile	= new PrintStream(new FileOutputStream(outputPath));
			PrintStream original	= System.out;
			
			System.setOut(outputFile);
			Array.print(output);
			System.setOut(original);
		} catch (FileNotFoundException e) {
			System.out.println("Cannot create trajectory file at this "
					+ "location: "+outputPath);
			e.printStackTrace();
		}
		System.out.println("Finished printing tables.");
	}
	
	public static void loadMarkovModel(ExperimentReference trainSet, ExperimentReference testSet, 
			int modelLength, String modelMethod) {
		//Load Markov Model
		MarkovModelInfo mm = SELEX.getMarkModelInfo(trainSet, testSet,  modelLength, null, modelMethod);
		//Markov Model properties
		modelObjPath 			= mm.getMarkovObjPath();
		String modelCountPath 	= mm.getMarkovCountsPath();
		modelTotalCount 		= mm.getMarkovLengthTotalCount();

		try {
			System.out.println("Reading Markov Prob file:" + modelObjPath);
			FileInputStream fos = new FileInputStream(modelObjPath);
			ObjectInputStream oos = new ObjectInputStream(fos);
			modelProbs = (float[]) oos.readObject();

			System.out.println("Reading Markov Count file:" + modelCountPath);
			ArraysMergerHeap.MinHeapNode mmNode = new ArraysMergerHeap.MinHeapNode(modelCountPath);
			CountObject obj = null;
			modelCounts = new int[1 << (2 * (modelLength))];

			while ((obj = mmNode.peek()) != null) {
				Sequence seq = (Sequence) obj.getKey();
				modelCounts[(int) seq.getValue()] = obj.getCount();
				mmNode.pop();
			}
			oos.close();
		} catch(Exception ex) {
			throw new RuntimeException(ex);
		}
		System.out.println("Finished loading Markov Model.");
	}
}
