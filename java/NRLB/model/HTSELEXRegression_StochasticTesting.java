package model;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.ExecutionException;

import base.*;
import base.Model.CompactGradientOutput;
import minimizers.*;

public class HTSELEXRegression_StochasticTesting {
	static boolean isVerbose, testShifts, storeHessian, errorBars, growMotif;
	static boolean hasSymmetry, printPSAM, printSeeds, printRaw, simpleSymmetry; 
	static boolean fitR0, multiMode, uniformK, crossValidate, useSeed;
	static int l, r0k, maxK, nThreads, lbfgsMem, lbfgsMaxIters, nShifts, nModes, round, nFolds;
	static double lambda = 0, lbfgsConvergence, testRatio, nsbSeed;
	static String nucSymmetry, dinucSymmetry, lFlank, rFlank, outputDir;
	static String dataDir, experimentName, trajectoryFileName, configFile;
	static boolean[] useDinuc, useSymmetry;
	static int[] flankLengths, r1TestKs, r0TestKs, startKs, maxKs;
	static double[] nucSeed;
	
	public static void main(String[] args) {
		boolean isFlank, isDinuc, hitMax;
		boolean notAllComplete = true, isShape = false, isNSBinding = true;
		int currFlank, startK;
		double likelihood, adjL;
		Shape shapeModel			= null;
		boolean[] activeModes;
		int[] currKs				= null;
		double[] seed, currPos;
		String[] nucSym = null, dinucSym = null;
		Fit[] currFits				= null;
		Data[] R0data;
		Round0Fit R0Fit;
		Minimizer minimizer;
		Round0Model R0Train, R0Test, R0Model;
		Round0Results R0Results;
		MultiModeModel mmModel;
		MultinomialFit bestFit		= null;
		MultinomialModel model;
		MultinomialResults results;
		ArrayList<Object[]> R1Data, R1DataTrainOnly;
		
		//Read and process config file, if it has been provided
		if (args.length<1) {
			throw new IllegalArgumentException("No configuration file defined!");
		}
		configFile = args[0];
		try {
			readConfigFile(configFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		//First fit R0 Model if requested
		if (fitR0) {
			R0data		= fileParserR0(dataDir+experimentName+".tsv");
			R0Results	= new Round0Results(l, "LBFGS", 0.0, 0.0, 0.0, false, 
						lbfgsMem, lbfgsMaxIters, lbfgsConvergence, true);
			R0Results.defineTrainSet(experimentName, "NA", 
					dataDir+experimentName+".tsv", R0data[0].nCount, lFlank, rFlank);
			R0Results.defineTestSet(true, experimentName, "NA", "NA", R0data[2].nCount);
			//Build Round0 Model with cross validation. Continue to fit until all have converged
			System.out.println("Currently fitting Round0 Model.");
			while (notAllComplete) {
				for (int k : r0TestKs) {
					System.out.println("Currently testing model for k = "+k);
					R0Train	= new Round0Model(R0data[1], k, true);
					R0Test  = new Round0Model(R0data[2], k, true);
					minimizer = new LBFGS(R0Train, lbfgsMem, lbfgsConvergence, 
							lbfgsMaxIters, true, false, false, isVerbose);					
					try {
						R0Fit	= (Round0Fit) minimizer.minimize(null, trajectoryFileName)[0];
						R0Test.setParams(R0Fit.betas);
						likelihood = R0Test.functionEval();
						adjL	= (likelihood-R0Test.maxLikelihood())*R0Test.likelihoodNormalizer();
						R0Train.replaceData(R0data[0]);
						R0Fit = (Round0Fit) minimizer.minimize(null, trajectoryFileName)[0];
						R0Fit.recordCrossValidate(likelihood, adjL);						
						R0Results.addFit(R0Fit);
					} catch (Exception e) {
						e.printStackTrace();
						System.out.println("Error while fitting Round0 model. Restarting process.");
						break;
					}
				}
				notAllComplete = false;
			}
			//Write results to file
			try {
				R0Results.store(outputDir+experimentName+"_R0", true);
				R0Results.print(outputDir+experimentName+"_R0", true, true);
				System.out.println("Finished fitting Round0 Model.");
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		//Set Round0 Model and read in Round 1 Probes
		System.out.println("Fitting Round1 Model.");
		R0Model	= new Round0Model(outputDir+experimentName+"_R0.dat", l, r0k, true, lFlank, rFlank);
		R1Data	= fileParserR1(dataDir+experimentName+".tsv", R0Model);
		
		//Define objects and fit on Round 1
		results	= new MultinomialResults(l, "LBFGS", 0.0, 0.0, 0.0, false, 
				lbfgsMem, lbfgsMaxIters, lbfgsConvergence, true, configFile);
		results.defineDataset(experimentName, "NA", dataDir+experimentName+".tsv", 
				lFlank, rFlank, ((Data) R1Data.get(0)[0]).nCount, false, 0, 0);
		results.defineR0Model(R0Model);
		results.autoSave(outputDir+"/"+experimentName+"_R1");
		
		if (multiMode) {
			//We do not want to perform cross-validation on intermediate steps
			R1DataTrainOnly = new ArrayList<Object[]>();
			R1DataTrainOnly.add(R1Data.get(0));
			if (uniformK) {				
				for (int flankLength : flankLengths) {
					for (int k : r1TestKs) {
						//Test to see if there are any windows that can accomodate said 
						//configuration; reject if k is greater than the total number of
						//windows or if l+f is greater than the capacity of a long
						if (l+2*flankLength-k+1<1 || l+2*flankLength > 32) {
							continue;	
						}
						System.out.println("Testing k "+k);
						//Define start and end K's
						startKs = new int[nModes];
						for (int i=0; i<nModes; i++) {
							startKs[i] = k;
						}
						
						isFlank = true;
						isDinuc = Array.containsTrue(useDinuc);
						System.out.println("Performing sequential fitting with starting flank length "+flankLength+" and k "+k);
						//Build base model, ensure isFlank is always true
						startK   = startKs[0];
						currFlank= flankLength;
						currKs   = new int[]{startK};
						nucSym   = null;
						//
						if (hasSymmetry && useSymmetry[0]) {
							nucSym = Array.cat(nucSym, evenOdd(startK));
						} else {
							nucSym = Array.cat(nucSym, "null");
						}
						if (!useSeed) {
							//Fit shifts and find best model
							mmModel = new MultiModeModel(nThreads, shapeModel, ((Data) R1Data.get(0)[0]), isFlank, currFlank, 
									false, isShape, isNSBinding, currKs, nucSym, null);
							mmModel.setLambda(lambda*mmModel.getNCount());
							
							minimizer = newMinimizer(mmModel);
							try {
								currFits = minimizer.shiftPermutation(null, trajectoryFileName, nShifts, R1DataTrainOnly, null);
							} catch (Exception e) {
								e.printStackTrace();
								mmModel.threadPoolShutdown();
								return;
							}
//							printFit(minFit(currFits));
							
//							mmModel.evaluatedDataPoints = 0;
//							minimizer = newMinimizer(mmModel);
//							minimizer.setXStar(minFit(currFits).finalPosition);
//							try {
//								currFits = minimizer.shiftPermutation(null, trajectoryFileName, nShifts, R1DataTrainOnly, null);
//							} catch (Exception e) {
//								e.printStackTrace();
//								mmModel.threadPoolShutdown();
//								return;
//							}
							mmModel.threadPoolShutdown();
							bestFit = minFit(currFits);
							printFit(bestFit);
						}
	
						for (int currMode=1; currMode<nModes; currMode++) {
							//Build seed & update model
							if (useSeed && currMode<2) {
								currPos = generateSeed(k);
							} else {
								currPos = bestFit.finalPosition;
								printFit(bestFit);
							}
							seed = Arrays.copyOfRange(currPos, 0, currPos.length-1);
							seed = Array.cat(seed, new double[4*startKs[currMode]]);
							seed = Array.cat(seed, currPos[currPos.length-1]);
							currKs = Array.cat(currKs, startKs[currMode]);
							if (hasSymmetry && useSymmetry[currMode]) {
								nucSym = Array.cat(nucSym, evenOdd(startKs[currMode]));
							} else {
								nucSym = Array.cat(nucSym, "null");
							}
							//Set active mode to the newest one
							activeModes = new boolean[currMode+1];
							for (int i=0; i<currMode; i++) {
								activeModes[i] = false;
							}
							activeModes[currMode] = true;
							mmModel = new MultiModeModel(nThreads, shapeModel, ((Data) R1Data.get(0)[0]), isFlank, currFlank, 
									false, isShape, isNSBinding, currKs, nucSym, null);
							mmModel.setLambda(lambda*mmModel.getNCount());
							mmModel.setParams(seed);
							mmModel.setActiveModes(activeModes);
							minimizer = newMinimizer(mmModel);
							//Fit shifts and find best model
							try {
								currFits = minimizer.shiftPermutation(seed, trajectoryFileName, nShifts, R1DataTrainOnly, null);
							} catch (Exception e) {
								e.printStackTrace();
								mmModel.threadPoolShutdown();
								return;
							}
							bestFit = minFit(currFits);
							seed = bestFit.finalPosition;
							printFit(bestFit);
							//Activate all modes & refit
							for (int i=0; i<currMode+1; i++) {
								activeModes[i] = true;
							}
							mmModel.setActiveModes(activeModes);
							try {
								if (currMode==nModes-1) {
									currFits = minimizer.minimize(seed, trajectoryFileName, R1Data);
								} else {
									currFits = minimizer.minimize(seed, trajectoryFileName);
								}
							} catch (Exception e) {
								e.printStackTrace();
								mmModel.threadPoolShutdown();
								return;
							}
							mmModel.threadPoolShutdown();
							bestFit = minFit(currFits);
							bestFit.finalPosition = Array.clone(mmModel.normalize());							
						}
						results.addFit(bestFit);
						printFit(bestFit);
						saveR1Results(results);
					}
				}	
			} else {
				for (int flankLength : flankLengths) {
					//Check if lengths are valid
					if (l+2*flankLength-Array.max(startKs)+1<1 || l+2*flankLength > 32) {
						continue;	
					}
					isFlank = true;
					isDinuc = Array.containsTrue(useDinuc);
					System.out.println("Performing sequential fitting with starting flank length "+flankLength);
					//Build base model, ensure isFlank is always true
					startK   = startKs[0];
					currFlank= flankLength;
					currKs   = new int[]{startK};
					nucSym   = null;
					//
					if (hasSymmetry && useSymmetry[0]) {
						nucSym = Array.cat(nucSym, evenOdd(startK));
					} else {
						nucSym = Array.cat(nucSym, "null");
					}
					//Fit shifts and find best model
					mmModel = new MultiModeModel(nThreads, shapeModel, ((Data) R1Data.get(0)[0]), isFlank, currFlank, 
							false, isShape, isNSBinding, currKs, nucSym, null);
					mmModel.setLambda(lambda*mmModel.getNCount());
					minimizer = newMinimizer(mmModel);
					try {
						currFits = minimizer.shiftPermutation(null, trajectoryFileName, nShifts, R1DataTrainOnly, null);
					} catch (Exception e) {
						e.printStackTrace();
						mmModel.threadPoolShutdown();
						return;
					}
					mmModel.threadPoolShutdown();
					bestFit = minFit(currFits);

					for (int currMode=1; currMode<nModes; currMode++) {
						//Build seed & update model
						currPos = bestFit.finalPosition;
						seed = Arrays.copyOfRange(currPos, 0, currPos.length-1);
						seed = Array.cat(seed, new double[4*startKs[currMode]]);
						seed = Array.cat(seed, currPos[currPos.length-1]);
						printFit(bestFit);
						currKs = Array.cat(currKs, startKs[currMode]);
						if (hasSymmetry && useSymmetry[currMode]) {
							nucSym = Array.cat(nucSym, evenOdd(startKs[currMode]));
						} else {
							nucSym = Array.cat(nucSym, "null");
						}
						//Set active mode to the newest one
						activeModes = new boolean[currMode+1];
						for (int i=0; i<currMode; i++) {
							activeModes[i] = false;
						}
						activeModes[currMode] = true;
						mmModel = new MultiModeModel(nThreads, shapeModel, ((Data) R1Data.get(0)[0]), isFlank, currFlank, 
								false, isShape, isNSBinding, currKs, nucSym, null);
						mmModel.setLambda(lambda*mmModel.getNCount());
						mmModel.setParams(seed);
						mmModel.setActiveModes(activeModes);
						minimizer = newMinimizer(mmModel);
						//Fit shifts and find best model
						try {
							currFits = minimizer.shiftPermutation(seed, trajectoryFileName, nShifts, R1DataTrainOnly, null);
						} catch (Exception e) {
							e.printStackTrace();
							mmModel.threadPoolShutdown();
							return;
						}
						bestFit = minFit(currFits);
						seed = bestFit.finalPosition;
						printFit(bestFit);
						//Activate all modes & refit
						for (int i=0; i<currMode+1; i++) {
							activeModes[i] = true;
						}
						mmModel.setActiveModes(activeModes);
						try {
							currFits = minimizer.minimize(seed, trajectoryFileName);
						} catch (Exception e) {
							e.printStackTrace();
							mmModel.threadPoolShutdown();
							return;
						}
						mmModel.threadPoolShutdown();
						bestFit = minFit(currFits);
						bestFit.finalPosition = Array.clone(mmModel.normalize());
					}
					results.addFit(bestFit);
					printFit(bestFit);
					saveR1Results(results);
					System.out.println("Finished performing sequential fitting for starting k.");
					while (true) {
						//Check to see if all modes have hit their largest length. 
						//If not, increase and build seeds simultaneously as well
						hitMax = true;
						nucSym = null;
						seed   = null;
						for (int i=0; i<nModes; i++) {
							if (currKs[i] < maxKs[i]) {
								hitMax = false;
								currKs[i] += 2;
								seed = Array.cat(seed, Array.cat(Array.cat(new double[4], mmModel.getNucBetas(i)), new double[4]));
							} else {
								seed = Array.cat(seed, mmModel.getNucBetas(i));
							}
							//Build symmetry string
							if (hasSymmetry && useSymmetry[i]) {
								nucSym = Array.cat(nucSym, evenOdd(currKs[i]));
							} else {
								nucSym = Array.cat(nucSym, "null");
							}
						}
						seed = Array.cat(seed, mmModel.getNSBeta());
						System.out.println("After growing");
						results.print(true, false, false);
						//Check to see if all modes have hit max length
						if (hitMax) {
							break;
						}
						//Grow flank length and check to see if it is within bounds
						currFlank++;
						if (l+2*currFlank-Array.max(currKs)+1<1 || l+2*currFlank> 32) {
							break;	
						}
						System.out.println("Fitting 2 additional positions from previous motif(s).");
						mmModel = new MultiModeModel(nThreads, shapeModel, ((Data) R1Data.get(0)[0]), isFlank, currFlank, 
								false, isShape, isNSBinding, currKs, nucSym, null);
						mmModel.setLambda(lambda*mmModel.getNCount());
						minimizer = newMinimizer(mmModel);
						//Fit new model (with additional positions)
						try {
							currFits = minimizer.minimize(seed, trajectoryFileName);
						} catch (Exception e) {
							e.printStackTrace();
							mmModel.threadPoolShutdown();
							return;
						}
						mmModel.threadPoolShutdown();
						bestFit = minFit(currFits);
						bestFit.finalPosition = Array.clone(mmModel.normalize());
						results.addFit(bestFit);
						System.out.println("After adding");
						results.print(true, false, false);
						printFit(bestFit);
						saveR1Results(results);
					}
					System.out.println("Finished growing motif(s).");
					
					//Add dinucleotides if requested
					if (isDinuc) {
						//Ensure normalized start point
						mmModel.setParams(mmModel.normalize());
						//Next, create the seed vector by adding dinucleotide 
						//parameters to nuc and shape (if they exist) and define 
						//the symmetry structure 
						seed	= null;
						dinucSym= null;
						for (int currMode=0; currMode<nModes; currMode++) {
							seed = Array.cat(seed, mmModel.getNucBetas(currMode));
							seed = Array.cat(seed, new double[16*(currKs[currMode]-1)]);
							if (isShape) {
								seed = Array.cat(seed, mmModel.getShapeBetas(currMode));
							}
							if (hasSymmetry && useSymmetry[currMode]) {
								dinucSym= Array.cat(dinucSym, evenOdd(currKs[currMode]-1));
							} else {
								dinucSym= Array.cat(dinucSym, "null");
							}
						}
						seed = Array.cat(seed, mmModel.getNSBeta());
						
						//Next, create new model and minimizer
						mmModel = new MultiModeModel(nThreads, shapeModel, ((Data) R1Data.get(0)[0]), isFlank, 
								currFlank, isDinuc, isShape, isNSBinding, currKs, nucSym, dinucSym);
						mmModel.setLambda(lambda*mmModel.getNCount());
						minimizer = newMinimizer(mmModel);
						System.out.println("Perfoming single-shot addition of dinucleotides.");

						try {
							currFits = minimizer.minimize(seed, trajectoryFileName, R1Data);
						} catch (Exception e) {
							e.printStackTrace();
							mmModel.threadPoolShutdown();
							return;
						}
						mmModel.threadPoolShutdown();
						bestFit = minFit(currFits);
						bestFit.finalPosition = mmModel.normalize();
						results.addFit(bestFit);
						printFit(bestFit);
						saveR1Results(results);
					}
				}
			}
		} else {
			//Perform regressions over all available option combinations. First
			//loop over flanks
			for (int flankLength : flankLengths) {						
				System.out.println("Testing flank length "+flankLength);
				//And K's. In the case of the 'growMotif' option, flankLengths and 
				//testKs represent the starting flank length and k. Ensure that isFlank is true
				isFlank = true;
				isDinuc = Array.containsTrue(useDinuc);
				for (int k : r1TestKs) {
					//Test to see if there are any windows that can accomodate said 
					//configuration; reject if k is greater than the total number of
					//windows or if l+f is greater than the capacity of a long
					if (l+2*flankLength-k+1<1 || l+2*flankLength > 32) {
						continue;	
					}
					System.out.println("Testing k "+k);
					//First build base nuc model
					symBuilder(k);
					model= new MultinomialModel(nThreads, shapeModel, ((Data) R1Data.get(0)[0]), k, isFlank, flankLength,
							false, false, isNSBinding, nucSymmetry, dinucSymmetry, false);
					model.setLambda(lambda*model.getNCount());
					minimizer = newMinimizer(model);
					currFits  = runFits(model, minimizer, results, R1Data, testShifts, generateSeed(k));
					//Grow dinucs if dinuc fit is requested
					if (isDinuc) {
						//Find optimal fit
						bestFit = minFit(currFits);
						currPos = bestFit.positionVector();
						seed = seedBuilder(k, Arrays.copyOfRange(currPos, 0, 4*k), null, true, 
								null, isShape, currPos, isNSBinding, shapeModel);
						model = new MultinomialModel(nThreads, shapeModel, ((Data) R1Data.get(0)[0]), k, isFlank, flankLength, 
								true, isShape, isNSBinding, nucSymmetry, dinucSymmetry, false);
						model.setLambda(lambda*model.getNCount());
						minimizer = newMinimizer(model);
						currFits  = runFits(model, minimizer, results, R1Data, testShifts, seed);
					}
					//See if the growMotif option is selected
					if (growMotif) {
						//Ensure that current k is less than max k
						if (k > maxK) {
							continue;
						}
						//Start looping and storage process
						currFlank = flankLength;
						for (int currK=k+2; currK<maxK; currK+=2) {
							//Grow flank length and ensure it is within bounds
							currFlank++;
							if (l+2*currFlank-currK+1<1 || 
									l+2*currFlank> 32 ||
									(isShape && l+2*currFlank+4 > 32)) {
								continue;	
							}
							System.out.println("Fitting "+currK+", 2 additional positions from previous motif");
							//Take optimal fit and create new seed and grow it
							bestFit = minFit(currFits);
							currPos = bestFit.positionVector();
							seed = seedBuilder(currK-2, currPos, isDinuc, isShape, isNSBinding, shapeModel);
							symBuilder(currK);
							model = new MultinomialModel(nThreads, shapeModel, ((Data) R1Data.get(0)[0]), currK, isFlank, 
									currFlank, isDinuc, isShape, isNSBinding, nucSymmetry, dinucSymmetry, false);
							model.setLambda(lambda*model.getNCount());
							minimizer = newMinimizer(model);
							currFits = runFits(model, minimizer, results, R1Data, testShifts, seed);
						}					
					}
				}
			}
		}
		System.out.println("Fitting process completed.");
	}
	
	private static MultinomialFit minFit(Fit[] input) {
		int currIdx = 0;
		double minVal = Double.MAX_VALUE;
		
		for (int i=0; i<input.length; i++) {
			if (input[i]!=null && input[i].trainLikelihood<minVal) {
				minVal = input[i].trainLikelihood;
				currIdx = i;
			}
		}
		return ((MultinomialFit) input[currIdx]);
	}
	
	private static Fit[] runFits(MultinomialModel model, Minimizer min, 
			MultinomialResults results, ArrayList<Object[]> data, 
			boolean testShifts, double[] seed) {
		Fit[] currFits = null;
		
		try {
			if (testShifts) {
				currFits = min.shiftPermutation(seed, trajectoryFileName, 
						nShifts, data, results);
			} else {
				currFits = min.minimize(seed, trajectoryFileName, data);
			}
			if (results!=null) {
				for (int i=0; i<currFits.length; i++) {
					if (currFits[i]!=null) {
						results.addFit(currFits[i]);
					}
				}
				saveR1Results(results);
			}
		} catch (Exception e) {
			e.printStackTrace();
			model.threadPoolShutdown();
		}
		model.threadPoolShutdown();
		return currFits;
	}
	
	private static void printFit(MultinomialFit fit) {
		System.out.println("-----Results for "+fit.betaIdx.length+" mode(s):");
		fit.print(true, true, true);
	}
	
	private static void saveR1Results(MultinomialResults results) {
		results.print(outputDir+"/"+experimentName+"_R1", printPSAM, 
				printSeeds, printRaw);
		results.printList(outputDir+"/"+experimentName+"_R1", false);
		results.printList(outputDir+"/"+experimentName+"_R1", true);
	}
	
	private static double[] generateSeed(int currK) {
		if (!useSeed) {
			return null;
		}		
		boolean isLeft = true;
		double[] output= Array.clone(nucSeed);
		
		while (true) {
			if (output.length == 4*currK) {
				break;
			} else {
				if (output.length > 4*currK) {
					//subtract positions
					if (isLeft) {
						output = Arrays.copyOfRange(output, 4, output.length);
					} else {
						output = Arrays.copyOfRange(output, 0, output.length-4);
					}
				} else {
					//add positions
					if (isLeft) {
						output = Array.cat(new double[4], output);
					} else {
						output = Array.cat(output, new double[4]);
					}
				}
				//switch between left and right
				isLeft = !isLeft;
			}
		}
		output = Array.cat(output, nsbSeed);
		
		return output;
	}
	
	private static double[] seedBuilder(int k, double[] pos, boolean isDinuc, 
			boolean isShape, boolean isNSBinding, Shape shapeModel) {
		double[] nucSeed	= Array.grow(Arrays.copyOfRange(pos, 0, 4*k), 4);
		double[] dinucSeed	= null, shapeSeed = null, nsBindingSeed = null;
		
		if (isDinuc) {
			dinucSeed = Array.grow(Arrays.copyOfRange(pos, 4*k, 4*k+16*(k-1)), 16);
			if (isShape) {
				shapeSeed = Array.grow(Arrays.copyOfRange(pos, 4*k+16*(k-1), 4*k+16*(k-1)+shapeModel.nShapeFeatures()*k), shapeModel.nShapeFeatures());
			}
		} else if (isShape) {
			shapeSeed = Array.grow(Arrays.copyOfRange(pos, 4*k, 4*k+shapeModel.nShapeFeatures()*k), shapeModel.nShapeFeatures());
		}
		if (isNSBinding) {
			nsBindingSeed = new double[]{pos[pos.length-1]};
		}
		return seedBuilder(k+2, nucSeed, dinucSeed, isDinuc, shapeSeed, isShape,
				nsBindingSeed, isNSBinding, shapeModel);
	}
	
	private static double[] seedBuilder(int k, double[] nucSeed, 
			double[] dinucSeed, boolean isDinuc, double[] shapeSeed, 
			boolean isShape, double[] nsBindingSeed, boolean isNSBinding, 
			Shape shapeModel) {
		double[] seed = null;
		if (nucSeed!=null || dinucSeed!=null || shapeSeed!=null || nsBindingSeed!=null) {
			seed = (nucSeed!=null) ? Array.clone(nucSeed) : new double[4*k];
			if (isDinuc) {
				seed = (dinucSeed!=null) ? Array.cat(seed, dinucSeed) : Array.cat(seed, new double[16*(k-1)]);
			}
			if (isShape) {
				seed = (shapeSeed!=null) ? Array.cat(seed, shapeSeed) : Array.cat(seed, new double[shapeModel.nShapeFeatures()*k]);
			}
			if (isNSBinding) {
				seed = (nsBindingSeed!=null) ? Array.cat(seed, nsBindingSeed[nsBindingSeed.length-1]) : Array.cat(seed, 0);
			}
		}
		return seed;
	}
	
	//hasSymmetry, useSymmetry
	private static void readConfigFile(String configFile) throws Exception {
		String currLine;
		String[] temp, temp2;
		ArrayList<String> config= new ArrayList<String>(100);
		BufferedReader br		= new BufferedReader(new FileReader(configFile));
		
		//First load all non-commented lines
		while ((currLine = br.readLine())!=null) {
			if (currLine.startsWith("//")) {
				continue;
			}
			config.add(currLine);
		}
		br.close();
		
		//Read in required values: data file info and output info
		nThreads			= Integer.parseInt(extractValue("nThreads", config, true));
		dataDir				= extractValue("dataDir", config, true);
		outputDir			= extractValue("outputDir", config, true);
		experimentName		= extractValue("experimentName", config, true);
		l					= Integer.parseInt(extractValue("l", config, true));
		trajectoryFileName	= extractValue("trajectoryFileName", config, false);
		if (trajectoryFileName!=null) {
			trajectoryFileName = outputDir + trajectoryFileName;
		}
		printPSAM			= Boolean.parseBoolean(extractValue("printPSAM", config, false));
		printSeeds			= Boolean.parseBoolean(extractValue("printSeeds", config, false));
		printRaw			= Boolean.parseBoolean(extractValue("printRaw", config, false));
		
		//Dataset info
		lFlank				= extractValue("lFlank", config, true);
		rFlank				= extractValue("rFlank", config, true);
		fitR0				= Boolean.parseBoolean(extractValue("fitR0", config, true));
		r0k					= Integer.parseInt(extractValue("r0k", config, true));
		
		//Minimizer Configuration
		isVerbose			= Boolean.parseBoolean(extractValue("isVerbose", config, true));
		lbfgsMem			= Integer.parseInt(extractValue("lbfgsMem", config, true));
		lbfgsMaxIters		= Integer.parseInt(extractValue("lbfgsMaxIters", config, true));
		lbfgsConvergence	= Double.parseDouble(extractValue("lbfgsConvergence", config, true));
		
		//Fit Configuration
		lambda				= Double.parseDouble(extractValue("lambda", config, false));
		r0TestKs			= extractArrayInt("r0TestKs", ",", config, true);
		r1TestKs			= extractArrayInt("r1TestKs", ",", config, true);
		flankLengths		= null;
		temp				= extractArray("useFlank", ",", config, true);
		for (int i=0; i<temp.length; i++) {
			if (Boolean.parseBoolean(temp[i])) {	//consider flanking sequences
				temp2		= extractArray("flankLength", ",", config, true);
				for (int j=0; j<temp2.length; j++) {
					flankLengths	= Array.cat(flankLengths, new int[]{Integer.parseInt(temp2[j])});
				}
			} else {
				flankLengths	= Array.cat(flankLengths, new int[]{0});
			}
		}
		growMotif			= Boolean.parseBoolean(extractValue("growMotif", config, true));
		if (growMotif) {
			maxK			= Integer.parseInt(extractValue("maxK", config, true));
		}
		useDinuc			= extractArrayBoolean("useDinuc", ",", config, true);
		storeHessian		= Boolean.parseBoolean(extractValue("storeHessian", config, true));
		errorBars			= Boolean.parseBoolean(extractValue("errorBars", config, true));
		testShifts			= Boolean.parseBoolean(extractValue("testShifts", config, true));
		if (testShifts)	nShifts = Integer.parseInt(extractValue("nShifts", config, true));
		simpleSymmetry		= Boolean.parseBoolean(extractValue("simpleSymmetry", config, false));
		nucSymmetry			= extractValue("nucSymmetry", config, false);
		dinucSymmetry		= extractValue("dinucSymmetry", config, false);
		multiMode			= Boolean.parseBoolean(extractValue("multiMode", config, false));
		if (multiMode) {
			startKs			= extractArrayInt("startKs", ",", config, true);
			uniformK		= Boolean.parseBoolean(extractValue("uniformK", config, false));
			maxKs			= extractArrayInt("maxKs", ",", config, true);
			nModes			= Integer.parseInt(extractValue("nModes", config, true));
			useSymmetry		= extractArrayBoolean("useSymmetry", ",", config, true);
			hasSymmetry		= Array.containsTrue(useSymmetry);
		}
		round				= Integer.parseInt(extractValue("round", config, true));
		
		//Cross-validation
		crossValidate		= Boolean.parseBoolean(extractValue("crossValidate", config, false));
		if (crossValidate) {
			nFolds			= Integer.parseInt(extractValue("nFolds", config, true));
			testRatio		= Double.parseDouble(extractValue("testRatio", config, true));
		}
		
		//Seeding
		useSeed				= Boolean.parseBoolean(extractValue("useSeed", config, true));
		if (useSeed) {
			nucSeed			= extractArrayDouble("nucSeed", ",", config, true);
			nsbSeed			= Double.parseDouble(extractValue("nsbSeed", config, true));
		}
	}
		
	private static String extractValue(String paramName, 
			ArrayList<String> config, boolean isRequired) {
		String output = null, currLine;
		String[] parsed;
		
		for (int i=0; i<config.size(); i++) {
			//remove all quotes if they exist
			currLine= config.get(i).replaceAll("\"", "");				
			currLine= config.get(i).replaceAll("\'", "");
			//split the first and second halves of the entry
			parsed	= currLine.split("=");
			//Perform case-insensitive matching
			if ((parsed[0].trim()).matches("(?i)"+paramName)) {
				//Remove trailing comments or semicolons
				output = (currLine.split("=")[1]).split("//|;")[0].trim();
			}
		}
		if (output==null && isRequired) {
			throw new IllegalArgumentException("Configuration File Does Not "
					+ "Contain The Required Entry: "+paramName);
		}
		if (!isRequired) {
			if (output!=null && output.matches("(?i)null|(?i)na|(?i)n/a|(?i)false|(?i)no")) {
				output = null;
			}
		}
		
		return output;
	}
		
	private static String[] extractArray(String paramName, String split, 
			ArrayList<String> config, boolean isRequired) {
		String parsed = extractValue(paramName, config, isRequired);
		
		if (parsed==null) {
			return null;
		} else {
			return parsed.replaceAll("\\s+", "").split(split);
		}
	}
	
	private static boolean[] extractArrayBoolean(String paramName, String split,
			ArrayList<String> config, boolean isRequired) {
		boolean[] output;
		String[] parsed = extractArray(paramName, split, config, isRequired);
		
		if (parsed==null) {
			return null;
		} else {
			output = new boolean[parsed.length];
			for (int i=0; i<parsed.length; i++) {
				output[i] = Boolean.parseBoolean(parsed[i]);
			}
			return output;
		}
	}
	
	public static int[] extractArrayInt(String paramName, String split, ArrayList<String> config, boolean isRequired) {
		int[] output;
		String[] parsed = extractArray(paramName, split, config, isRequired);
		
		if (parsed==null) {
			return null;
		} else {
			output = new int[parsed.length];
			for (int i=0; i<parsed.length; i++) {
				output[i] = Integer.parseInt(parsed[i]);
			}
			return output;
		}
	}
	
	public static double[] extractArrayDouble(String paramName, String split, ArrayList<String> config, boolean isRequired) {
		double[] output;
		String[] parsed = extractArray(paramName, split, config, isRequired);
		
		if (parsed==null) {
			return null;
		} else {
			output = new double[parsed.length];
			for (int i=0; i<parsed.length; i++) {
				output[i] = Double.parseDouble(parsed[i]);
			}
			return output;
		}
	}
	
	//TODO
	private static Minimizer newMinimizer(Model model) {
//		return new LBFGS(model, lbfgsMem, lbfgsConvergence, lbfgsMaxIters, true, errorBars, storeHessian, isVerbose);
//		return new sLBFGS(model, 20, 5000, 200, 1000, 20, 400, .005, 1E-7, 0, isVerbose);
//		return new sLBFGS(model, 500, 5000, 50, 1000, 20, 200, .03, 1E-7, 0, isVerbose);
//		return new sLBFGS_PLS(model, 500, 5000, 50, 1000, 20, 200, .02, 1E-7, 0, isVerbose);
		return new sLBFGS_DynStep(model, 500, 5000, 50, 1000, 20, 200, .02, 1E-7, 0, isVerbose);
		//return new kSVRG(model, 10, true, true, 50, 1000, .01, lbfgsConvergence, false, isVerbose);
//		return new sLBFGS_kSVRG(model, 5, 500, 5000, 1000, 20, 50, .03, 1E-7, 0, isVerbose);
//		return new sLBFGS_kSVRG(model, 10, 500, 2500, 1000, 5, 200, .02, lbfgsConvergence, 1E-4, isVerbose);
//		return new sLBFGS_kSVRG_V2(model, 10, 500, 5000, 1000, 10, 80, .01, lbfgsConvergence, 0, isVerbose);
//		return new sLBFGS_kSVRG_V2(model, 5, 5, 10, 1000, 10, 80, .005, lbfgsConvergence, 0, isVerbose);
		//return new sLBFGS_kSVRG(model, 5, 1000, 5000, 1000, 10, 40, .02, 1E-7, 0, isVerbose);
		//return new SGD(model, false, 100, 1000, false, 500, .005, 1E-5, false, true, true);
	}
	
	private static void symBuilder(int k) {
		if (!simpleSymmetry) {
			return;
		}
		nucSymmetry   = evenOdd(k);
		dinucSymmetry = evenOdd(k-1);
	}
	
	private static String evenOdd(int k) {
		String output = "1";
		boolean isEven = (k%2==0);
		
		if (isEven) {
			for (int i=1; i<k/2; i++) {
				output = output.concat(",1");
			}
			for (int i=0; i<k/2; i++) {
				output = output.concat(",-1");
			}
		} else {
			for (int i=1; i<(k-1)/2; i++) {
				output = output.concat(",1");
			}
			output = output.concat(",*");
			for (int i=0; i<(k-1)/2; i++) {
				output = output.concat(",-1");
			}
		}
		return output;
	}
	
	private static Data[] fileParserR0(String filePath) {
		int currCount, nTrainCounts= 0, nTestCounts = 0, nFullCounts = 0;
		String currLine;
		Sequence currSeq;
		BufferedReader br;
		Random generator				= new Random();
		int[] trainC, testC, fullC;
		long[] trainS, testS, fullS;
		Data[] output					= new Data[3];
		String[] temp;
		ArrayList<Integer> trainCounts	= new ArrayList<Integer>();
		ArrayList<Integer> testCounts	= new ArrayList<Integer>();
		ArrayList<Sequence> trainSeq	= new ArrayList<Sequence>();
		ArrayList<Sequence> testSeq		= new ArrayList<Sequence>();
		
		//First build R0 counts 
		try {
			br = new BufferedReader(new FileReader(filePath));
			//Read in all lines
			while ((currLine = br.readLine())!=null) {
				temp = currLine.split("\t");
				//First entry is the sequence
				currSeq = new Sequence(temp[0], 0, l);
				//Get Round 0 Count;
				currCount = Integer.parseInt(temp[1]);
				if (currCount>0) {
					if (generator.nextBoolean()) {
						trainSeq.add(currSeq);
						trainCounts.add(currCount);							
					} else {
						testSeq.add(currSeq);
						testCounts.add(currCount);
					}
				}
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		//Split dataset into test/train pairs and full dataset
		trainC	= new int[trainCounts.size()];
		testC	= new int[testCounts.size()]; 
		fullC	= new int[trainCounts.size()+testCounts.size()];
		trainS	= new long[trainSeq.size()];
		testS	= new long[testSeq.size()];
		fullS	= new long[trainSeq.size()+testSeq.size()];
		
		//Convert to Arrays
		for (int i=0; i<trainCounts.size(); i++) {
			trainC[i]	= trainCounts.get(i);
			nTrainCounts+= trainCounts.get(i);
			trainS[i]	= trainSeq.get(i).getValue();
			fullC[i]	= trainCounts.get(i);
			fullS[i]	= trainSeq.get(i).getValue();
		}
		for (int i=0; i<testCounts.size(); i++) {
			testC[i]	= testCounts.get(i);
			nTestCounts += testCounts.get(i);
			testS[i]	= testSeq.get(i).getValue();
			fullC[i+trainCounts.size()] = testCounts.get(i);
			fullS[i+trainCounts.size()] = testSeq.get(i).getValue();
		}
		nFullCounts = nTrainCounts+nTestCounts;
		
		//Store datasets in following order: Full R0, Train R0, Test R0
		output[0] = new Data(l, nFullCounts, lFlank, rFlank, fullC, fullS, null, null);
		output[1] = new Data(l, nTrainCounts, lFlank, rFlank, trainC, trainS, null, null);
		output[2] = new Data(l, nTestCounts, lFlank, rFlank, testC, testS, null, null);
		return output;
	}
	
	private static ArrayList<Object[]> fileParserR1(String filePath, Round0Model R0) {
		int currCount, nTrainCount = 0, nTestCount = 0, nTotCount = 0;
		String currLine;
		Random generator;
		Sequence currSeq;
		BufferedReader br;
		int[] fCount, counts;
		long[] fProbe, probes;
		double[] R0Prob, probs;
		String[] temp;
		Data[] out						= new Data[1];
		ArrayList<Long> testProbes, trainProbes;
		ArrayList<Integer> testCounts, trainCounts;
		ArrayList<Double> testProb, trainProb;
		trainCounts						= new ArrayList<Integer>();
		trainProbes						= new ArrayList<Long>();
		ArrayList<Object[]> output		= new ArrayList<Object[]>();
		
		//First build R1 counts 
		try {
			br = new BufferedReader(new FileReader(filePath));
			//Read in all lines
			while ((currLine = br.readLine())!=null) {
				temp = currLine.split("\t");
				//First entry is the sequence
				currSeq = new Sequence(temp[0], 0, l);
				//Get Round 1 Count;
				currCount = Integer.parseInt(temp[round+1]);
				if (currCount>0) {
					trainProbes.add(currSeq.getValue());
					trainCounts.add(currCount);							
				}
				
//				//TODO: Terminates file reading
//				if (trainProbes.size()>99) {
//					break;
//				}
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		fCount	= new int[trainCounts.size()];
		fProbe	= new long[trainProbes.size()];
		R0Prob  = new double[trainProbes.size()];
		
		//Convert to Arrays
		for (int i=0; i<trainCounts.size(); i++) {
			fCount[i]	= trainCounts.get(i);
			fProbe[i]	= trainProbes.get(i);
			R0Prob[i]	= R0.getSeqProb(trainProbes.get(i));
			nTotCount	+= trainCounts.get(i);
		}
		
		out[0] = new Data(l, nTotCount, lFlank, rFlank, fCount, fProbe, R0Prob, R0);
		output.add(out);
		
		//Now construct cross-validation set(s)
		if(crossValidate) {
			int trainSet, testSet;
			generator		= new Random();
			
			//Loop over the number of folds
			for (int i=0; i<nFolds; i++) {			
				nTrainCount	= 0;
				nTestCount	= 0;
				trainProbes	= new ArrayList<Long>(10*1000*1000);
				trainCounts	= new ArrayList<Integer>(10*1000*1000);
				trainProb	= new ArrayList<Double>(10*1000*1000);
				testProbes	= new ArrayList<Long>(10*1000*1000);
				testCounts	= new ArrayList<Integer>(10*1000*1000);
				testProb	= new ArrayList<Double>(10*1000*1000);
				out			= new Data[2];

				for (int j=0; j<fCount.length; j++) {
					trainSet	= 0;
					testSet		= 0;
					for (int nReads=0; nReads<fCount[j]; nReads++) {
						//Stash a percentage of reads away as the test set
						if (generator.nextDouble() < testRatio) {
							testSet++;
						} else {
							trainSet++;
						}
					}
					if (testSet!=0) {
						nTestCount += testSet;
						testProbes.add(fProbe[j]);
						testCounts.add(testSet);
						testProb.add(R0Prob[j]);
					}
					if (trainSet!=0) {
						nTrainCount += trainSet;
						trainProbes.add(fProbe[j]);
						trainCounts.add(trainSet);
						trainProb.add(R0Prob[j]);
					}
				}

				//store train data
				probes = new long[trainProbes.size()];
				counts = new int[trainProbes.size()];
				probs  = new double[trainProbes.size()];
				for (int k=0; k<trainProbes.size(); k++) {
					probes[k] = trainProbes.get(k);
					counts[k] = trainCounts.get(k);
					probs[k]  = trainProb.get(k);
				}	
				out[0] = new Data(l, nTrainCount, lFlank, rFlank, counts, probes, probs, R0);
				//store test data
				probes = new long[testProbes.size()];
				counts = new int[testProbes.size()];
				probs  = new double[testProbes.size()];
				for (int k=0; k<testProbes.size(); k++) {
					probes[k] = testProbes.get(k);
					counts[k] = testCounts.get(k);
					probs[k]  = testProb.get(k);
				}	
				out[1] = new Data(l, nTestCount, lFlank, rFlank, counts, probes, probs, R0);
				output.add(out);
			}
		}
		return output;
	}
}
