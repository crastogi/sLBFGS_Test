package model;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;

import base.*;
import minimizers.*;
import utils.FileUtilities;

public class MultiModeRegression{
	static boolean isVerbose, crossValidate, isR0Flank, testShifts, hasSymmetry; 
	static boolean psRandomAxis, lbfgsMCSearch, storeHessian, errorBars, growMotif;
	static boolean printPSAM, printSeeds, printRaw, filterReads, growDinuc, singleShot;
	static int nThreads, l, lbfgsMem, lbfgsMaxIters, r0k, nShifts;
	static int nFolds, maxBaseCount, maxBaseRepeat, nModes;
	static double testRatio, psInitStep, psTheta, psConvergence;
	static double lambda = 0, lbfgsConvergence;
	static String configFile, sequencingRunName, sampleName, samplePath, lFlank;
	static String rFlank, minimizerType, nucSymmetry, dinucSymmetry, shapeDir;
	static String R0ModelPath, outputLocation, outputDataName, trajectoryFileName;
	static boolean[] useDinuc, useNSBinding, evenOdd, useSymmetry;
	static int[] flankLengths, testKs, ks, startK, maxK;
	static double[] nucSeed, dinucSeed, shapeSeed, nsBindingSeed, fullSeed;
	static String[] regexFilter;
	static ArrayList<String[]> shapes;
	static ArrayList<double[]> modeSeeds = null;
	
	public static void main(String[] args) {
		boolean isShape, isFlank, currDinuc, hitMax;
		int mk = 0, currFlank, currK;
		boolean[] activeModes, activeDinucs;
		int[] currKs = null;
		double[] pos, seed;
		String[] nucSym = null, dinucSym = null;
		Shape shapeModel;
		Minimizer minimizer;
		Round0Model R0Model;
		MultiModeModel test;
		FileUtilities reader;
		Fit[] currFits				= null;
		MultinomialResults results	= null;
		ArrayList<Object[]> datasets;
		MultinomialFit latestFit;
		
		//Read and process config file
		if (args.length<1) {	//Has a configuration file been provided?
			throw new IllegalArgumentException("No configuration file defined!");
		}
		configFile = args[0];
		try {
			readConfigFile(configFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		//Define objects
		R0Model	= new Round0Model(R0ModelPath, l, r0k, isR0Flank, lFlank, rFlank);
		results	= new MultinomialResults(l, minimizerType, psInitStep, psTheta, psConvergence, psRandomAxis,
				lbfgsMem, lbfgsMaxIters, lbfgsConvergence, lbfgsMCSearch, configFile);
		reader	= new FileUtilities(l);
		datasets= reader.readSeqFile(samplePath, R0Model, lFlank, rFlank, 
				crossValidate, nFolds, testRatio, filterReads, maxBaseCount,
				maxBaseRepeat, regexFilter);
		
		results.defineDataset(sequencingRunName, sampleName, samplePath, lFlank, 
				rFlank, ((Data) datasets.get(0)[0]).nCount, crossValidate, testRatio, nFolds);
		if (filterReads)	results.defineFilter(maxBaseCount, maxBaseRepeat, reader.removedReads, regexFilter);
		results.defineR0Model(R0Model);
		if (shapeDir!=null)	results.defineShapeModel(shapeDir);
		results.autoSave(outputLocation+"/"+outputDataName);
		
		//Perform regressions over all available option combinations
		for (String[] currShape : shapes) {								//Loop over shape features
			if (currShape[0].equals("0")) {
				isShape		= false;
				shapeModel	= null;
			} else {
				isShape		= true;
				shapeModel	= new Shape(shapeDir, true);
				shapeModel.setFeatures(currShape);
			}
			for (int flankLength : flankLengths) {					//Loop over flanks
				//Test to see if flank length + kmer lengths to be fit are valid
				for (int k : ks) {
					if (k > mk) {
						mk = k;
					}
				}
				if (l+2*flankLength-mk+1<1 || l+2*flankLength > 32 ||
						(isShape && l+2*flankLength+4 > 32)) {
					continue;	
				}	
				isFlank		= (flankLength==0) ? false : true;
				System.out.println("Testing flank length "+flankLength);
				for (boolean isDinuc : useDinuc) {					//Loop over dinucs
					for (boolean isNSBinding : useNSBinding) {		//Loop over NS Binding
						if (growMotif) {
							System.out.println("Performing sequential fitting.");
							//Build base model, ensure isFlank is always true
							isFlank	 = true;
							currK    = startK[0];
							currFlank= flankLength;
							currKs   = new int[]{currK};
							nucSym   = null;
							if (hasSymmetry && useSymmetry[0]) {
								nucSym = Array.cat(nucSym, evenOdd(currK, (currK%2==0) ));
							} else {
								nucSym = Array.cat(nucSym, "null");
							}
							//Fit shifts and find best model
							test = new MultiModeModel(nThreads, shapeModel, ((Data) datasets.get(0)[0]), isFlank, currFlank, 
									false, isShape, isNSBinding, currKs, nucSym, null);
							test.setLambda(lambda*test.getNCount());
							minimizer = new LBFGS(test, lbfgsMem, lbfgsConvergence, lbfgsMaxIters, lbfgsMCSearch, errorBars, storeHessian, isVerbose);
							try {
								currFits = minimizer.shiftPermutation(null, trajectoryFileName, nShifts, datasets, null);
							} catch (Exception e) {
								e.printStackTrace();
								test.threadPoolShutdown();
								return;
							}
							test.threadPoolShutdown();
							latestFit = minFit(currFits);

							for (int currMode=1; currMode<nModes; currMode++) {
								//Build seed & update model
								pos = latestFit.finalPosition;
								seed = (isNSBinding) ? Arrays.copyOfRange(pos, 0, pos.length-1) : Arrays.copyOfRange(pos, 0, pos.length);
								seed = Array.cat(seed, new double[4*startK[currMode]]);
								if (isNSBinding) {
									seed = Array.cat(seed, pos[pos.length-1]);
								}
								fitPrint(latestFit);
								currKs = Array.cat(currKs, startK[currMode]);
								if (hasSymmetry && useSymmetry[currMode]) {
									nucSym = Array.cat(nucSym, evenOdd(startK[currMode], (startK[currMode]%2==0) ));
								} else {
									nucSym = Array.cat(nucSym, "null");
								}
								//Set active mode to the newest one
								activeModes = new boolean[currMode+1];
								for (int i=0; i<currMode; i++) {
									activeModes[i] = false;
								}
								activeModes[currMode] = true;
								test = new MultiModeModel(nThreads, shapeModel, ((Data) datasets.get(0)[0]), isFlank, currFlank, 
										false, isShape, isNSBinding, currKs, nucSym, null);
								test.setLambda(lambda*test.getNCount());
								test.setParams(seed);
								test.setActiveModes(activeModes);
								minimizer = new LBFGS(test, lbfgsMem, lbfgsConvergence, lbfgsMaxIters, lbfgsMCSearch, errorBars, storeHessian, isVerbose);
								//Fit shifts and find best model
								try {
									currFits = minimizer.shiftPermutation(seed, trajectoryFileName, nShifts, datasets, null);
								} catch (Exception e) {
									e.printStackTrace();
									test.threadPoolShutdown();
									return;
								}
								latestFit = minFit(currFits);
								seed = latestFit.finalPosition;
								fitPrint(latestFit);
								//Activate all modes & refit
								for (int i=0; i<currMode+1; i++) {
									activeModes[i] = true;
								}
								test.setActiveModes(activeModes);
								try {
									currFits = minimizer.minimize(seed, trajectoryFileName);
								} catch (Exception e) {
									e.printStackTrace();
									test.threadPoolShutdown();
									return;
								}
								test.threadPoolShutdown();
								latestFit = minFit(currFits);
								latestFit.finalPosition = test.normalize();
							}
							fitPrint(latestFit);
							results.addFit(latestFit);
							System.out.println("Finished performing sequential fitting for starting k.");
							while (true) {
								//Check to see if all modes have hit their largest length. If not, increase
								//Build seeds simultaneously as well
								hitMax = true;
								nucSym = null;
								seed   = null;
								for (int i=0; i<nModes; i++) {
									if (currKs[i] < maxK[i]) {
										hitMax = false;
										currKs[i] += 2;
										seed = Array.cat(seed, Array.cat(Array.cat(new double[4], test.getNucBetas(i)), new double[4]));
									} else {
										seed = Array.cat(seed, test.getNucBetas(i));
									}
									//Build symmetry string
									if (hasSymmetry && useSymmetry[i]) {
										nucSym = Array.cat(nucSym, evenOdd(currKs[i], (currKs[i]%2==0) ));
									} else {
										nucSym = Array.cat(nucSym, "null");
									}
								}
								if (isNSBinding) {
									seed = Array.cat(seed, test.getNSBeta());
								}
								//Check to see if all modes have hit max length
								if (hitMax) {
									break;
								}
								//Grow flank length and check to see if it is within bounds
								currFlank++;
								if (l+2*currFlank-currK+1<1 || l+2*currFlank> 32 || (isShape && l+2*currFlank+4 > 32)) {
									break;	
								}
								System.out.println("Fitting 2 additional positions from previous motif(s).");
								test = new MultiModeModel(nThreads, shapeModel, ((Data) datasets.get(0)[0]), isFlank, currFlank, 
										false, isShape, isNSBinding, currKs, nucSym, null);
								test.setLambda(lambda*test.getNCount());
								minimizer = new LBFGS(test, lbfgsMem, lbfgsConvergence, lbfgsMaxIters, lbfgsMCSearch, errorBars, storeHessian, isVerbose);
								//Fit new model (with additional positions)
								try {
									currFits = minimizer.minimize(seed, trajectoryFileName);
								} catch (Exception e) {
									e.printStackTrace();
									test.threadPoolShutdown();
									return;
								}
								latestFit = minFit(currFits);
								latestFit.finalPosition = test.normalize();
								fitPrint(latestFit);
								results.addFit(latestFit);
								test.threadPoolShutdown();
							}
							results.print(outputLocation+"/"+outputDataName, printPSAM, printSeeds, printRaw);
							results.printList(outputLocation+"/"+outputDataName, false);
							results.printList(outputLocation+"/"+outputDataName, true);
							System.out.println("Finished growing motif(s).");
							
							//Add dinucleotides if requested
							if (isDinuc) {
								//Ensure normalized start point
								test.setParams(test.normalize());
								
								//Next, create the seed vector by adding 
								//dinucleotide parameters to nuc and shape 
								//(if they exist) and define the symmetry 
								//structure 
								seed	= null;
								dinucSym= null;
								for (int currMode=0; currMode<nModes; currMode++) {
									seed = Array.cat(seed, test.getNucBetas(currMode));
									seed = Array.cat(seed, new double[16*(currKs[currMode]-1)]);
									if (isShape) {
										seed = Array.cat(seed, test.getShapeBetas(currMode));
									}
									if (hasSymmetry && useSymmetry[currMode]) {
										dinucSym= Array.cat(dinucSym, evenOdd(currKs[currMode]-1, !(currKs[currMode]%2==0)));
									} else {
										dinucSym= Array.cat(dinucSym, "null");
									}
								}
								if (isNSBinding) {
									seed = Array.cat(seed, test.getNSBeta());
								}
								
								//Next, create new model and minimizer
								test = new MultiModeModel(nThreads, shapeModel, ((Data) datasets.get(0)[0]), isFlank, 
										currFlank, isDinuc, isShape, isNSBinding, currKs, nucSym, dinucSym);
								test.setLambda(lambda*test.getNCount());
								minimizer = new LBFGS(test, lbfgsMem, lbfgsConvergence, lbfgsMaxIters, lbfgsMCSearch, 
										errorBars, storeHessian, isVerbose);
								System.out.println("Perfoming single-shot addition of dinucleotides.");

								try {
									currFits = minimizer.minimize(seed, trajectoryFileName);
								} catch (Exception e) {
									e.printStackTrace();
									test.threadPoolShutdown();
									return;
								}
								test.threadPoolShutdown();
								latestFit = minFit(currFits);
								latestFit.finalPosition = test.normalize();
								fitPrint(latestFit);
								results.addFit(latestFit);
								results.print(outputLocation+"/"+outputDataName, printPSAM, printSeeds, printRaw);
								results.printList(outputLocation+"/"+outputDataName, false);
								results.printList(outputLocation+"/"+outputDataName, true);
							}
						} 						
						//Single-shot fitting
						else if (singleShot && !growDinuc) {
							System.out.println("Performing Single-Shot fitting.");
							//Build base model and seed
							nucSym	= null;
							dinucSym= null;
							if (hasSymmetry) {
								for (int i=0; i<ks.length; i++) {
									if (useSymmetry[i]) {
										nucSym = Array.cat(nucSym, evenOdd(ks[i], evenOdd[i]));
										if (isDinuc) {
											dinucSym = Array.cat(dinucSym, evenOdd(ks[i]-1, !evenOdd[i]));
										}
									} else {
										nucSym = Array.cat(nucSym, "null");
										if (isDinuc)	dinucSym = Array.cat(dinucSym, "null");
									}
								}								
							}
							test = new MultiModeModel(nThreads, shapeModel, ((Data) datasets.get(0)[0]), isFlank, flankLength,
									isDinuc, isShape, isNSBinding, ks, nucSym, dinucSym);
							test.setLambda(lambda*test.getNCount());
							
							//First, fit without any shifts
							minimizer = new LBFGS(test, lbfgsMem, lbfgsConvergence, lbfgsMaxIters, lbfgsMCSearch, errorBars, storeHessian, isVerbose);
							try {
								if (fullSeed!=null) {
									currFits = minimizer.minimize(fullSeed, trajectoryFileName);
								} else {
									currFits = minimizer.minimize(null, trajectoryFileName);
								}
								latestFit= minFit(currFits);
								latestFit.finalPosition = test.normalize();
							} catch (Exception e) {
								e.printStackTrace();
								test.threadPoolShutdown();
								return;
							}
							fitPrint(latestFit);
							test.threadPoolShutdown();
							results.addFit(latestFit);
							results.print(outputLocation+"/"+outputDataName, printPSAM, printSeeds, printRaw);
							results.printList(outputLocation+"/"+outputDataName, false);
							results.printList(outputLocation+"/"+outputDataName, true);					
						} //Modes have been seeded; need to perform mode-regression first
						else if (modeSeeds!=null) {
							System.out.println("Performing Mode-Mode regression to resolve relative affinities of the seeds.");
							//Build base model and seed
							seed	= null;
							nucSym	= null;
							dinucSym= null;
							for (int i=0; i<ks.length; i++) {
								if (hasSymmetry) {
									if (useSymmetry[i]) {
										nucSym = Array.cat(nucSym, evenOdd(ks[i], evenOdd[i]));
										if (isDinuc) {
											dinucSym = Array.cat(dinucSym, evenOdd(ks[i]-1, !evenOdd[i]));
										}
									} else {
										nucSym = Array.cat(nucSym, "null");
										if (isDinuc)	dinucSym = Array.cat(dinucSym, "null");
									}
								}
								seed = Array.cat(seed, modeSeeds.get(i));
							}
							if (isNSBinding) {
								seed = Array.cat(seed, nsBindingSeed);
							}
							//First, perform mode-mode regression
							test = new MultiModeModel(nThreads, shapeModel, ((Data) datasets.get(0)[0]), isFlank, flankLength,
										isDinuc, isShape, isNSBinding, ks, nucSym, dinucSym);
							test.setLambda(lambda*test.getNCount());
							test.setParams(seed);
							test.setModeRegression(true);
							minimizer = new LBFGS(test, lbfgsMem, lbfgsConvergence, lbfgsMaxIters, lbfgsMCSearch, errorBars, storeHessian, isVerbose);
							try {
								currFits = minimizer.minimize(null, trajectoryFileName);
								latestFit= minFit(currFits);
								latestFit.finalPosition = test.getMergedPositionVector();
							} catch (Exception e) {
								e.printStackTrace();
								test.threadPoolShutdown();
								return;
							}
							System.out.println("Switching to regular regression.");
							test.setModeRegression(false);
							seed = latestFit.finalPosition;
							minimizer = new LBFGS(test, lbfgsMem, lbfgsConvergence, lbfgsMaxIters, lbfgsMCSearch, errorBars, storeHessian, isVerbose);
							try {
								test.setParams(seed);
								currFits = minimizer.minimize(seed, trajectoryFileName);
							} catch (Exception e) {
								e.printStackTrace();
								test.threadPoolShutdown();
								return;
							}
							test.threadPoolShutdown();
							latestFit= minFit(currFits);
							latestFit.finalPosition = test.normalize();
						} else {
							System.out.println("Performing sequential fitting.");
							//Check to see if dinucs will be used and grown
							currDinuc = isDinuc && !growDinuc;
							//Build base model
							currKs   = null;
							nucSym   = null;
							dinucSym = null;
							currKs   = Array.cat(currKs, ks[0]);
							if (hasSymmetry) {
								if (useSymmetry[0]) {
									nucSym = Array.cat(nucSym, evenOdd(ks[0], evenOdd[0]));
									if (currDinuc) {
										dinucSym = Array.cat(dinucSym, evenOdd(ks[0]-1, !evenOdd[0]));
									}
								} else {
									nucSym = Array.cat(nucSym, "null");
									if (currDinuc)	dinucSym = Array.cat(dinucSym, "null");
								}
							}
							//Fit shifts and find best model
							test = new MultiModeModel(nThreads, shapeModel, ((Data) datasets.get(0)[0]), isFlank, flankLength, 
									currDinuc, isShape, isNSBinding, currKs, nucSym, dinucSym);
							test.setLambda(lambda*test.getNCount());
							minimizer = new LBFGS(test, lbfgsMem, lbfgsConvergence, lbfgsMaxIters, lbfgsMCSearch, errorBars, storeHessian, isVerbose);
							try {
								currFits = minimizer.shiftPermutation(null, trajectoryFileName, nShifts, datasets, null);
							} catch (Exception e) {
								e.printStackTrace();
								test.threadPoolShutdown();
								return;
							}
							test.threadPoolShutdown();
							latestFit = minFit(currFits);

							for (int currMode=1; currMode<ks.length; currMode++) {
								//Build seed & update model
								pos = latestFit.finalPosition;
								seed = (isNSBinding) ? Arrays.copyOfRange(pos, 0, pos.length-1) : Arrays.copyOfRange(pos, 0, pos.length);
								seed = Array.cat(seed, new double[4*ks[currMode]]);
								if (currDinuc) {
									seed = Array.cat(seed, new double[16*(ks[currMode]-1)]);
								}
								if (isNSBinding) {
									seed = Array.cat(seed, pos[pos.length-1]);
								}
								fitPrint(latestFit);
								currKs = Array.cat(currKs, ks[currMode]);
								if (hasSymmetry) {
									if (useSymmetry[currMode]) {
										nucSym = Array.cat(nucSym, evenOdd(ks[currMode], evenOdd[currMode]));
										if (currDinuc) {
											dinucSym = Array.cat(dinucSym, evenOdd(ks[currMode]-1, !evenOdd[currMode]));
										}
									} else {
										nucSym = Array.cat(nucSym, "null");
										if (currDinuc)	dinucSym = Array.cat(dinucSym, "null");
									}
								}
								//Set active mode to the newest one
								activeModes = new boolean[currMode+1];
								for (int i=0; i<currMode; i++) {
									activeModes[i] = false;
								}
								activeModes[currMode] = true;
								test = new MultiModeModel(nThreads, shapeModel, ((Data) datasets.get(0)[0]), isFlank, flankLength, 
										currDinuc, isShape, isNSBinding, currKs, nucSym, dinucSym);
								test.setLambda(lambda*test.getNCount());
								test.setParams(seed);
								test.setActiveModes(activeModes);
								minimizer = new LBFGS(test, lbfgsMem, lbfgsConvergence, lbfgsMaxIters, lbfgsMCSearch, errorBars, storeHessian, isVerbose);
								//Fit shifts and find best model
								try {
									currFits = minimizer.shiftPermutation(seed, trajectoryFileName, nShifts, datasets, null);
								} catch (Exception e) {
									e.printStackTrace();
									test.threadPoolShutdown();
									return;
								}
								latestFit = minFit(currFits);
								seed = latestFit.finalPosition;
								fitPrint(latestFit);
								//Activate all modes & refit
								for (int i=0; i<currMode+1; i++) {
									activeModes[i] = true;
								}
								test.setActiveModes(activeModes);
								try {
									currFits = minimizer.minimize(seed, trajectoryFileName);
								} catch (Exception e) {
									e.printStackTrace();
									test.threadPoolShutdown();
									return;
								}
								test.threadPoolShutdown();
								latestFit = minFit(currFits);
								latestFit.finalPosition = test.normalize();
							}
							
							//Now try to 'grow' dinucleotides, if requested
							if (isDinuc && growDinuc) {
								//First save the current fit status
								fitPrint(latestFit);
								results.addFit(latestFit);
								results.print(outputLocation+"/"+outputDataName, printPSAM, printSeeds, printRaw);
								results.printList(outputLocation+"/"+outputDataName, false);
								results.printList(outputLocation+"/"+outputDataName, true);
								
								//Ensure normalized start point
								test.setParams(test.normalize());
								
								//Next, create the seed vector by adding 
								//dinucleotide parameters to nuc and shape 
								//(if they exist) and define the symmetry 
								//structure (will be fixed since there are no 
								//new modes being added)
								seed	= null;
								nucSym	= null;
								dinucSym= null;
								for (int currMode=0; currMode<ks.length; currMode++) {
									seed = Array.cat(seed, test.getNucBetas(currMode));
									seed = Array.cat(seed, new double[16*(ks[currMode]-1)]);
									if (isShape) {
										seed = Array.cat(seed, test.getShapeBetas(currMode));
									}
									if (hasSymmetry) {
										if (useSymmetry[currMode]) {
											nucSym	= Array.cat(nucSym, evenOdd(ks[currMode], evenOdd[currMode]));
											dinucSym= Array.cat(dinucSym, evenOdd(ks[currMode]-1, !evenOdd[currMode]));
										} else {
											nucSym	= Array.cat(nucSym, "null");
											dinucSym= Array.cat(dinucSym, "null");
										}
									}
								}
								if (isNSBinding) {
									seed = Array.cat(seed, test.getNSBeta());
								}
								
								//Next, create new model and minimizer (these 
								//will be the only ones used)
								test = new MultiModeModel(nThreads, shapeModel, ((Data) datasets.get(0)[0]), isFlank, 
										flankLength, isDinuc, isShape, isNSBinding, currKs, nucSym, dinucSym);
								test.setLambda(lambda*test.getNCount());
								minimizer = new LBFGS(test, lbfgsMem, lbfgsConvergence, lbfgsMaxIters, lbfgsMCSearch, 
										errorBars, storeHessian, isVerbose);
								
								if (singleShot) {
									System.out.println("Single-shot addition of dinucleotides.");

									try {
										currFits = minimizer.minimize(seed, trajectoryFileName);
									} catch (Exception e) {
										e.printStackTrace();
										test.threadPoolShutdown();
										return;
									}
									latestFit = minFit(currFits);
									latestFit.finalPosition = test.normalize();
								} else {
									System.out.println("Sequentially adding dinucleotides.");
									activeModes = new boolean[ks.length];
									activeDinucs= new boolean[ks.length];								
																		
									//Next, loop over all modes, adding dinucs at each step
									for (int currMode=0; currMode<ks.length; currMode++) {								
										//First define active modes and active dinucs									
										//Inital fit is with nuc+dinuc freedom for single mode
										for (int i=0; i<ks.length; i++) {
											activeModes[i] = false;
											activeDinucs[i]= false;
										}
										activeModes[currMode] = true;
										activeDinucs[currMode]= true;
										
										//Set parameters and actives
										test.setParams(seed);
										test.setActiveModes(activeModes);
										test.setActiveDinucs(activeDinucs);
										//Fit shifts and find best model
										try {
											currFits = minimizer.shiftPermutation(seed, trajectoryFileName, nShifts, datasets, null);
										} catch (Exception e) {
											e.printStackTrace();
											test.threadPoolShutdown();
											return;
										}
										latestFit = minFit(currFits);
										seed = latestFit.finalPosition;
										fitPrint(latestFit);
										
										//Activate all modes but NOT all dinucs and refit
										for (int i=0; i<ks.length; i++) {
											activeModes[i] = true;
											if (i<=currMode) {
												activeDinucs[i] = true;
											}
										}
										test.setActiveModes(activeModes);
										test.setActiveDinucs(activeDinucs);
										try {
											currFits = minimizer.minimize(seed, trajectoryFileName);
										} catch (Exception e) {
											e.printStackTrace();
											test.threadPoolShutdown();
											return;
										}
										latestFit = minFit(currFits);
										latestFit.finalPosition = test.normalize();
										seed = latestFit.finalPosition;
										fitPrint(latestFit);
									}
								}
								//Process completed, shut down threads.
								test.threadPoolShutdown();
							}
						}
						//Print and store latest fit
//						fitPrint();
						results.addFit(latestFit);
						results.print(outputLocation+"/"+outputDataName, printPSAM, printSeeds, printRaw);
						results.printList(outputLocation+"/"+outputDataName, false);
						results.printList(outputLocation+"/"+outputDataName, true);
					}
				}
			}
		}
		System.out.println("Process complete");
	}
	
	private static String evenOdd(int k, boolean isEven) {
		String output = "1";
		
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
	
	private static void fitPrint(MultinomialFit fit) {
		System.out.println("-----Results for "+fit.betaIdx.length+" mode(s):");
		fit.print(true, true, true);
	}
	
	public static void readConfigFile(String configFile) throws Exception {
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
		sequencingRunName	= extractValue("sequencingRunName", config, true);
		sampleName			= extractValue("sampleName", config, true);
		R0ModelPath			= extractValue("R0ModelDir", config, true)+"/"+extractValue("R0ModelName", config, true);
		l					= Integer.parseInt(extractValue("l", config, true));
		samplePath			= extractValue("selexWorkingDir", config, true)+"/"+sequencingRunName+"."+sampleName+"."+l+".dat";
		outputLocation		= extractValue("outputLocation", config, true);
		outputDataName		= extractValue("outputFileName", config, false);
		if (outputDataName==null) {
			outputDataName	= sequencingRunName+"_"+sampleName;
		}
		trajectoryFileName	= extractValue("trajectoryFileName", config, false);
		if (trajectoryFileName!=null) {
			trajectoryFileName = outputLocation + trajectoryFileName;
		}
		printPSAM			= Boolean.parseBoolean(extractValue("printPSAM", config, false));
		printSeeds			= Boolean.parseBoolean(extractValue("printSeeds", config, false));
		printRaw			= Boolean.parseBoolean(extractValue("printRaw", config, false));
		
		//Dataset info
		lFlank				= extractValue("lFlank", config, true);
		rFlank				= extractValue("rFlank", config, true);
		r0k					= Integer.parseInt(extractValue("r0k", config, true));
		isR0Flank			= Boolean.parseBoolean(extractValue("isR0Flank", config, true));
		filterReads			= Boolean.parseBoolean(extractValue("filterReads", config, false));
		if (filterReads) {
			maxBaseRepeat	= Integer.parseInt(extractValue("maxBaseRepeat", config, true));
			maxBaseCount	= Integer.parseInt(extractValue("maxBaseCount", config, true));
			regexFilter		= extractArray("regexFilter", ",", config, true);
		}
		
		//Minimizer Configuration
		isVerbose			= Boolean.parseBoolean(extractValue("isVerbose", config, true));
		minimizerType		= extractValue("minimizerType", config, true);
		if (minimizerType.matches("(?i)lbfgs")) {			//case insensitive matching
			minimizerType	= "LBFGS";
			lbfgsMem		= Integer.parseInt(extractValue("lbfgsMem", config, true));
			lbfgsMaxIters	= Integer.parseInt(extractValue("lbfgsMaxIters", config, true));
			lbfgsConvergence= Double.parseDouble(extractValue("lbfgsConvergence", config, true));
			lbfgsMCSearch	= Boolean.parseBoolean(extractValue("lbfgsMCSearch", config, true));
		}
		if (minimizerType.matches("(?i)ps|(?i)patternsearch|(?i)pattern search")) {
			minimizerType	= "Pattern Search";
			psInitStep		= Double.parseDouble(extractValue("psInitStep", config, true));
			psTheta			= Double.parseDouble(extractValue("psTheta", config, true));
			psConvergence	= Double.parseDouble(extractValue("psConvergence", config, true));
			psRandomAxis	= Boolean.parseBoolean(extractValue("psRandomAxis", config, true));
		}
		
		//Fit Configuration
		crossValidate		= Boolean.parseBoolean(extractValue("crossValidate", config, true));
		if (crossValidate)	{
			testRatio		= Double.parseDouble(extractValue("testRatio", config, true));
			currLine		= extractValue("nFolds", config, false);
			if (currLine==null) {
				nFolds = 1;
			} else {
				nFolds		= Integer.parseInt(currLine);
				if (nFolds<1)	nFolds = 1;
			}
		} 
		temp				= extractArray("testKs", ",", config, true);
		lambda				= Double.parseDouble(extractValue("lambda", config, false));
		testKs				= new int[temp.length];
		for (int i=0; i<temp.length; i++) {
			testKs[i]		= Integer.parseInt(temp[i]);
		}
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
		useDinuc			= extractArrayBoolean("useDinuc", ",", config, true);
		if (Array.containsTrue(useDinuc)) {
			growDinuc		= Boolean.parseBoolean(extractValue("growDinuc", config, true));
		}
		shapeDir			= extractValue("shapeDir", config, false);
		shapes				= new ArrayList<String[]>();
		temp				= extractArray("useShape", ",", config, true);
		for (int i=0; i<temp.length; i++) {
			if (Boolean.parseBoolean(temp[i])) {	//consider shape features
				shapeDir	= extractValue("shapeDir", config, true);
				temp2		= extractArray("shapes", "\\|", config, true);
				for (int j=0; j<temp2.length; j++) {
					shapes.add(temp2[j].split(","));
				}				
			} else {								//do NOT consider shape features
				shapes.add(new String[]{"0"});
			}
		}
		useNSBinding		= extractArrayBoolean("useNSBinding", ",", config, true);
		storeHessian		= Boolean.parseBoolean(extractValue("storeHessian", config, true));
		errorBars			= Boolean.parseBoolean(extractValue("errorBars", config, true));
		testShifts			= Boolean.parseBoolean(extractValue("testShifts", config, true));
		if (testShifts)	nShifts = Integer.parseInt(extractValue("nShifts", config, true));
		nucSymmetry			= extractValue("nucSymmetry", config, false);
		dinucSymmetry		= extractValue("dinucSymmetry", config, false);
		nucSeed				= extractArrayDouble("nucSeed", ",", config, false);
		dinucSeed			= extractArrayDouble("dinucSeed", ",", config, false);
		shapeSeed			= extractArrayDouble("shapeSeed", ",", config, false);
		nsBindingSeed		= extractArrayDouble("nsBindingSeed", ",", config, false);
		if (nsBindingSeed!=null) {
			nsBindingSeed	= new double[]{nsBindingSeed[0]};
		}
		
		ks					= extractArrayInt("ks", ",", config, true);
		useSymmetry			= extractArrayBoolean("useSymmetry", ",", config, true);
		hasSymmetry			= false;
		for (boolean currVal : useSymmetry) {
			if (currVal) {
				hasSymmetry = true;
				break;
			}
		}
		if(hasSymmetry)		evenOdd	= extractArrayBoolean("evenOdd", ",", config, true);
		if(extractArrayDouble("modeSeed1", ",", config, false)!=null) {
			modeSeeds		= new ArrayList<double[]>();
			for (int i=1; i<ks.length+1; i++) {
				modeSeeds.add(extractArrayDouble("modeSeed"+i,",", config, true));
			}
		}
		singleShot			= Boolean.parseBoolean(extractValue("singleShot", config, true));
		growMotif			= Boolean.parseBoolean(extractValue("growMotif", config, true));
		if (growMotif) {
			startK			= extractArrayInt("startK", ",", config, true);
			maxK			= extractArrayInt("maxK", ",", config, true);
			nModes			= Integer.parseInt(extractValue("nModes", config, true));
		}
		fullSeed			= extractArrayDouble("fullSeed", ",", config, false);
	}
	
	public static String extractValue(String paramName, ArrayList<String> config, boolean isRequired) {
		String output = null, currLine;
		String[] parsed;
		
		for (int i=0; i<config.size(); i++) {
			currLine= config.get(i).replaceAll("\"", "");					//remove all quotes if they exist
			currLine= config.get(i).replaceAll("\'", "");
			parsed	= currLine.split("=");									//split the first and second halves of the entry
			if ((parsed[0].trim()).matches("(?i)"+paramName)) {				//Perform case-insensitive matching
				output = (currLine.split("=")[1]).split("//|;")[0].trim();	//Remove trailing comments or semicolons
			}
		}
		if (output==null && isRequired) {
			throw new IllegalArgumentException("Configuration File Does Not Contain The Required Entry: "+paramName);
		}
		if (!isRequired) {
			if (output!=null && output.matches("(?i)null|(?i)na|(?i)n/a|(?i)false|(?i)no")) {
				output = null;
			}
		}
		
		return output;
	}
	
	public static String[] extractArray(String paramName, String split, ArrayList<String> config, boolean isRequired) {
		String parsed = extractValue(paramName, config, isRequired);
		
		if (parsed==null) {
			return null;
		} else {
			return parsed.replaceAll("\\s+", "").split(split);
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
	
	public static boolean[] extractArrayBoolean(String paramName, String split, ArrayList<String> config, boolean isRequired) {
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
}
