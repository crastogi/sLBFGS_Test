package model;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;

import base.*;
import minimizers.*;
import utils.FileUtilities;

public class MultinomialRegression {
	static boolean isVerbose, crossValidate, isR0Flank, testShifts; 
	static boolean psRandomAxis, lbfgsMCSearch, storeHessian, errorBars;
	static boolean printPSAM, printSeeds, printRaw, filterReads;
	static boolean growDinuc, growMotif, simpleSymmetry;
	static int nThreads, l, lbfgsMem, lbfgsMaxIters, r0k, nShifts;
	static int nFolds, maxBaseCount, maxBaseRepeat, maxK;
	static double testRatio, psInitStep, psTheta, psConvergence;
	static double lambda = 0, lbfgsConvergence;
	static String configFile, sequencingRunName, sampleName, samplePath, lFlank;
	static String rFlank, R0ModelPath, nucSymmetry, shapeDir, trajectoryFileName;
	static String minimizerType, outputLocation, outputDataName, dinucSymmetry;
	static boolean[] useDinuc, useNSBinding;
	static int[] flankLengths, testKs;
	static double[] nucSeed, dinucSeed, shapeSeed, nsBindingSeed;
	static String[] regexFilter;
	static ArrayList<String[]> shapes;
	
	public static void main(String[] args) {
		boolean isShape, isFlank;
		int currIdx = 0, currFlank;
		double minVal;
		double[] seed, currPos;
		Round0Model R0Model;
		Shape shapeModel;
		FileUtilities reader;
		ArrayList<Object[]> datasets;
		MultinomialModel model;
		Fit[] currFits				= null;
		MultinomialResults results	= null;
		Minimizer minimizer;
		
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
		results	= new MultinomialResults(l, minimizerType, psInitStep, psTheta, 
				psConvergence, psRandomAxis, lbfgsMem, lbfgsMaxIters, 
				lbfgsConvergence, lbfgsMCSearch, configFile);
		reader	= new FileUtilities(l);
		datasets= reader.readSeqFile(samplePath, R0Model, lFlank, rFlank, 
				crossValidate, nFolds, testRatio, filterReads, maxBaseCount,
				maxBaseRepeat, regexFilter);
		
		results.defineDataset(sequencingRunName, sampleName, samplePath, lFlank, 
				rFlank, ((Data) datasets.get(0)[0]).nCount, crossValidate, 
				testRatio, nFolds);
		if (filterReads)	results.defineFilter(maxBaseCount, maxBaseRepeat, 
				reader.removedReads, regexFilter);
		results.defineR0Model(R0Model);
		if (shapeDir!=null)	results.defineShapeModel(shapeDir);
		results.autoSave(outputLocation+"/"+outputDataName);
		
		//Perform regressions over all available option combinations. First
		//loop over shape features
		for (String[] currShape : shapes) {
			if (currShape[0].equals("0")) {
				isShape		= false;
				shapeModel	= null;
			} else {
				isShape		= true;
				shapeModel	= new Shape(shapeDir, true);
				shapeModel.setFeatures(currShape);
			}
			//Next, loop over NSBinding
			for (boolean isNSBinding : useNSBinding) {
				//Now loop over flanks
				for (int flankLength : flankLengths) {						
					isFlank		= (flankLength==0) ? false : true;
					System.out.println("Testing flank length "+flankLength);
					//And K's. In the case of the 'growMotif' option, 
					//flankLengths and testKs represent the starting flank 
					//length and k.
					for (int k : testKs) {
						//Test to see if there are any windows that can 
						//accomodate said configuration; reject if k is greater
						//than the total number of windows or if l+f is greater
						//than the capacity of a long
						if (l+2*flankLength-k+1<1 || l+2*flankLength > 32 ||
								(isShape && l+2*flankLength+4 > 32)) {
							continue;	
						}
						System.out.println("Testing k "+k);
						//Now see if the growMotif option is selected
						if (growMotif) {
							//Ensure that current k is less than max k
							if (k > maxK) {
								continue;
							}
							//Ensure that isFlank is true
							isFlank = true;
							//Loop over dinucs
							for (boolean isDinuc : useDinuc) {
								//Check if growDinuc is used for the initial fit
								if (isDinuc && growDinuc) {
									//growDinuc used; build nuc seed
									seed = seedBuilder(k, false, false, isNSBinding, null);
									symBuilder(k);
									model= new MultinomialModel(nThreads, shapeModel, ((Data) datasets.get(0)[0]), k, isFlank, flankLength,
											false, false, isNSBinding, nucSymmetry, dinucSymmetry, false);
									model.setLambda(lambda*model.getNCount());
									minimizer = newMinimizer(model);
									currFits  = runFits(model, minimizer, null, datasets, testShifts, seed);
									//Find optimal fit
									minVal = Double.MAX_VALUE;
									for (int i=0; i<currFits.length; i++) {
										if (currFits[i]!=null && currFits[i].trainLikelihood<minVal) {
											minVal	= currFits[i].trainLikelihood;
											currIdx = i;
										}
									}
									currPos = currFits[currIdx].positionVector();
									seed = seedBuilder(k, Arrays.copyOfRange(currPos, 0, 4*k), null, true, 
											shapeSeed, isShape, currPos, isNSBinding, shapeModel);
									model = new MultinomialModel(nThreads, shapeModel, ((Data) datasets.get(0)[0]), k, isFlank, flankLength, 
											true, isShape, isNSBinding, nucSymmetry, dinucSymmetry, false);
									model.setLambda(lambda*model.getNCount());
									minimizer = newMinimizer(model);
									currFits  = runFits(model, minimizer, results, datasets, testShifts, seed);
								} else {
									//growDinuc not used; grow as per usual
									//First build seed for base k
									seed = seedBuilder(k, isDinuc, isShape, isNSBinding, shapeModel);
									symBuilder(k);
									//Fit on the whole dataset, and if cross validation is requested, perform that as well
									model= new MultinomialModel(nThreads, shapeModel, ((Data) datasets.get(0)[0]), k, isFlank, flankLength, 
											isDinuc, isShape, isNSBinding, nucSymmetry, dinucSymmetry, false);
									model.setLambda(lambda*model.getNCount());
									minimizer = newMinimizer(model);
									currFits = runFits(model, minimizer, results, datasets,testShifts, seed);
								}
								//Start looping and storage process
								currFlank = flankLength;
								for (int currK=k+2; currK<maxK; currK+=2) {
									//Grow flank length and check to see if it
									//is within bounds
									currFlank++;
									if (l+2*currFlank-currK+1<1 || 
											l+2*currFlank> 32 ||
											(isShape && l+2*currFlank+4 > 32)) {
										continue;	
									}
									System.out.println("Fitting "+currK+", 2 additional positions from previous motif");
									//Find optimal fit from the previous results
									minVal = Double.MAX_VALUE;
									for (int i=0; i<currFits.length; i++) {
										if (currFits[i]!=null && currFits[i].trainLikelihood<minVal) {
											minVal	= currFits[i].trainLikelihood;
											currIdx = i;
										}
									}
									//Take optimal fit and create new seed and grow it
									currPos = currFits[currIdx].positionVector();
									seed = seedBuilder(currK-2, currPos, isDinuc, isShape, isNSBinding, shapeModel);
									symBuilder(currK);
									model = new MultinomialModel(nThreads, shapeModel, ((Data) datasets.get(0)[0]), currK, isFlank, 
											currFlank, isDinuc, isShape, isNSBinding, nucSymmetry, dinucSymmetry, false);
									model.setLambda(lambda*model.getNCount());
									minimizer = newMinimizer(model);
									currFits = runFits(model, minimizer, results, datasets, testShifts, seed);
								}
							}
						} else {
							//growMotif option not used. Check to see if the 
							//growDinuc option is selected and if the useDinuc 
							//vector contains a true element
							if (Array.containsTrue(useDinuc) && growDinuc) {
								seed = seedBuilder(k, false, false, isNSBinding, null);
								symBuilder(k);
								model= new MultinomialModel(nThreads, shapeModel, ((Data) datasets.get(0)[0]), k, isFlank, flankLength,
										false, false, isNSBinding, nucSymmetry, dinucSymmetry, false);
								model.setLambda(lambda*model.getNCount());
								minimizer = newMinimizer(model);
								currFits  = runFits(model, minimizer, results, datasets, testShifts, seed);
								//Find optimal fit
								minVal = Double.MAX_VALUE;
								for (int i=0; i<currFits.length; i++) {
									if (currFits[i]!=null && currFits[i].trainLikelihood<minVal) {
										minVal	= currFits[i].trainLikelihood;
										currIdx = i;
									}
								}
								currPos = currFits[currIdx].positionVector();
								seed = seedBuilder(k, Arrays.copyOfRange(currPos, 0, 4*k),null, true, shapeSeed, 
										isShape, currPos, isNSBinding, shapeModel);
								symBuilder(k);
								model = new MultinomialModel(nThreads, shapeModel, ((Data) datasets.get(0)[0]), k, isFlank, flankLength, 
										true, isShape, isNSBinding, nucSymmetry, dinucSymmetry, false);
								model.setLambda(lambda*model.getNCount());
								minimizer = newMinimizer(model);
								runFits(model, minimizer, results, datasets, testShifts, seed);
							} else {
								//growDinuc option not used. Loop over dinucs
								//and build seeds as per usual
								for (boolean isDinuc : useDinuc) {
									seed = seedBuilder(k, isDinuc, isShape, isNSBinding, shapeModel);
									symBuilder(k);
									//Fit on the whole dataset, and if cross validation is requested, perform that as well
									model= new MultinomialModel(nThreads, shapeModel, ((Data) datasets.get(0)[0]), k, isFlank, flankLength, 
											isDinuc, isShape, isNSBinding, nucSymmetry, dinucSymmetry, false);
									model.setLambda(lambda*model.getNCount());
									minimizer = newMinimizer(model);
									runFits(model, minimizer, results, datasets,testShifts, seed);
								}
							}
						}
					}
				}
			}
		}
		System.out.println("Fitting process completed.");
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
				results.print(outputLocation+"/"+outputDataName, printPSAM, 
						printSeeds, printRaw);
				results.printList(outputLocation+"/"+outputDataName, false);
				results.printList(outputLocation+"/"+outputDataName, true);	
			}
		} catch (Exception e) {
			e.printStackTrace();
			model.threadPoolShutdown();
		}
		model.threadPoolShutdown();
		return currFits;
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
	
	private static double[] seedBuilder(int k, boolean isDinuc, boolean isShape, 
			boolean isNSBinding, Shape shapeModel) {
		return seedBuilder(k, nucSeed, dinucSeed, isDinuc, shapeSeed, isShape,
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
		growMotif			= Boolean.parseBoolean(extractValue("growMotif", config, true));
		if (growMotif) {
			maxK			= Integer.parseInt(extractValue("maxK", config, true));
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
		currLine			= extractValue("simpleSymmetry", config, false);
		simpleSymmetry		= (currLine==null) ? false : Boolean.parseBoolean(currLine);
		nucSymmetry			= extractValue("nucSymmetry", config, false);
		dinucSymmetry		= extractValue("dinucSymmetry", config, false);
		nucSeed				= extractArrayDouble("nucSeed", ",", config, false);
		dinucSeed			= extractArrayDouble("dinucSeed", ",", config, false);
		shapeSeed			= extractArrayDouble("shapeSeed", ",", config, false);
		nsBindingSeed		= extractArrayDouble("nsBindingSeed", ",", config, false);
		if (nsBindingSeed!=null) {
			nsBindingSeed	= new double[]{nsBindingSeed[0]};
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
	
	
	private static double[] extractArrayDouble(String paramName, String split, 
			ArrayList<String> config, boolean isRequired) {
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
	
	
	private static Minimizer newMinimizer(Model model) {
		if (minimizerType.equals("LBFGS")) {
			return new LBFGS(model, lbfgsMem, lbfgsConvergence, lbfgsMaxIters, 
					lbfgsMCSearch, errorBars, storeHessian, isVerbose);
		} else {
			return new PatternSearch(model, psInitStep, psConvergence, psTheta, 
					psRandomAxis, errorBars, storeHessian, isVerbose);
		}
	}
	
	private static void symBuilder(int k) {
		if (!simpleSymmetry) {
			return;
		}
		String output = "1";
		
		if (k%2==0) {
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
		nucSymmetry = output;
		output = "1";
		k--;
		if (k%2==0) {
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
		dinucSymmetry = output;
	}
}
