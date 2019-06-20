package utils;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;

import base.Sequence;
import base.Shape;
import model.MultinomialResults;
import model.Round0Model;

public class MaxSeq {
	public static String dir		= "/Users/chaitanya/Documents/Research/SELEX/Biological Analysis/";
	public static String R0Dir		= "./Round0/";
	public static String shapeDir	= "./../Shape/";
	public static String results	= dir+"HM Dataset/Fits/Antp.dat";
	public static int modelIndex	= 23;
	public static int modeIndex		= -1;
	public static int l				= 16;
	public static double maxValue	= 1;
	public static double divSize	= .1;
	public static String bestMotif	= "AAAAACACCTAAATCTTT";
	public static int nMotifs		= 5;
	public static boolean save		= false;
	public static String savePath	= "";
	
	@SuppressWarnings("unchecked")
	public static void main(String[] args) {
		boolean unfinished			= true;
		boolean unconverged			= true;
		int maxBinCount				= 1000;
		int k, nDivs, currBin, nextOffset, counter=0, currLim=4;
		long currSeq, bestSeq, mask, addMask, topSeq;
		Shape shapeModel;
		Random generator			= new Random();
		Round0Model R0;
		ScoringObject so;
		MultinomialResults mr		= new MultinomialResults(results);
		ArrayList<Long> tempSeq;
		int[] binCounts;
		long[] testSeqs				= null;
		double[] weight, outWeight	= null;
		String[] parsedPath;
		HashSet<Long>[] probeSet;
		ArrayList<Long> probes 		= new ArrayList<Long>();
		PrintStream outputFile;
		
//		mr.l = 16;
//		mr.lFlank = "TGGGCCTGG";
//		mr.rFlank = "CCAGGGAGGTGGAGTAGG";
//		mr.getFit(modelIndex).isDinuc = true;
//		mr.getFit(modelIndex).k = 11;
//		mr.getFit(modelIndex).finalPosition = new double[]{-0.218678632, -0.274194455, -0.133516615, -0.041342229, -0.292684844, -0.1462958, -0.458435976, 0.229684689, 0.123128894, -0.039655521, 0.138763202, -0.889968507, 0.227621764, 0.017077599, 0.107340996, -1.019772291, 0.162714728, -0.2319694, -0.368079965, -0.230397294, -0.09146499, -0.257700797, -0.257700797, -0.09146499, -0.230397294, -0.368079965, -0.2319694, 0.162714728, -1.019772291, 0.107340996, 0.017077599, 0.227621764, -0.889968507, 0.138763202, -0.039655521, 0.123128894, 0.229684689, -0.458435976, -0.1462958, -0.292684844, -0.041342229, -0.133516615, -0.274194455, -0.218678632, -0.066793968, 0.10690226, -0.017766486, 0.222211095, -0.241921829, -0.266221268, -0.434527707, 0.023034984, 0.342877081, 0.371016862, 0.232581613, 0.359324141, -0.453406131, -0.316883913, -0.281847081, 0.060288083, 0.048259861, -0.011854843, 0.120899855, -0.57654972, -0.081110555, -0.268447345, 0.217343477, 0.027028364, 0.292814235, 0.046503092, 0.233004404, -1.073881391, -2.89E-05, 0.156767764, -0.094301254, 0.602420669, -0.195744144, 0.057412456, -0.231587642, 0.629853995, -0.147821264, -0.171796111, -0.269693305, 0.512279348, -0.21371515, -0.050921993, 0.065492688, 0.676090937, -0.457517792, -0.483835964, -0.527244019, 0.447615696, -0.440596656, -0.23480326, -0.159081385, -0.180317048, -0.353646162, -0.068436919, -0.154327532, -0.072730999, -0.434373788, -0.162669398, -0.27371141, -0.092277682, 0.396254771, 0.721156202, 0.592891174, 0.55553783, -0.141626433, -0.217842687, -0.16280972, -0.310082994, 0.143987465, 0.216598944, -0.05894595, -0.046393834, 0.159333912, 0.062592087, -0.287646263, 0.07149111, 0.094404492, 0.140295828, 0.040169482, -0.064657702, -0.066786656, -0.013170211, 0.02596393, 0.310092372, 0.136427299, -0.114319752, -0.063262167, 0.242798793, -0.418151992, -0.04919822, -0.356480765, 0.354598527, -0.279309184, -0.187181161, -0.078180505, 0.19502743, -1.421158841, 0.097255167, 0.400389357, 0.295693784, -2.083819611, 0.843308832, 0.362121135, 0.5145203, -1.127726305, 0.266635274, 0.023328241, 0.365803283, 1.077430868, 0.00329286, -0.051178207, 0.072971601, -1.535017708, -0.814100653, -0.585948653, -0.620206876, 0.676140195, 0.096900506, 0.249855859, 0.187595572, 0.433572946, 0.171006143, 0.140256041, -0.010174605, 0.78550148, 0.155312918, 0.179217178, 0.128957392, 0.139558211, 0.139856649, -0.004697018, 0.08547907, -0.055409795, -0.342243336, -0.028481957, 0.035254003, -0.117392591, 0.079296069, 0.00997056, 0.011506388, -0.142713711, 0.051578193, -0.108971797, -0.113721202, 0.11115481, -0.38617747, 0.101222394, -0.022198643, 0.257978386, -0.31417156, 0.153512084, -0.131573816, 0.190376038, -0.379324236, 0.192744901, -0.056713817, 0.294702949, -0.269520768, 0.270269105, -0.171136925};
		
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
		
		so			= new ScoringObject(mr, modelIndex, modeIndex, true, false, 
							shapeModel, R0);
		k			= so.getK();
		nDivs		= (int) Math.ceil(maxValue/divSize)+1;

		//Create unique list of binned sequences
		bestSeq		= (new Sequence(bestMotif, 0, k)).getValue();
		topSeq		= bestSeq;
		currSeq		= bestSeq;
		
		while(unconverged) {
			currLim		= 4;
			mask		= (long) Math.pow(4, k) - 1;
			addMask 	= currLim - 1;
			binCounts	= new int[nDivs];
			probeSet	= (HashSet<Long>[]) new HashSet[nDivs];
			probes 		= new ArrayList<Long>();
			for (int i=0; i<nDivs; i++) {
				probeSet[i] = new HashSet<Long>();
			}

			//Generate random sequences until bin conditions are satisfied
			testSeqs	= mutator(bestSeq, k);
			//Dispatch & Bin
			System.out.println(so.seqEval(bestSeq));
			weight		= so.seqEval(testSeqs);
			for (int i=0; i<testSeqs.length; i++) {
				currBin	= (int) (Math.floor((weight[i])/divSize));
				probeSet[currBin].add(testSeqs[i]);	
			}
			for (int i=0; i<nDivs; i++) {
				binCounts[i] = probeSet[i].size();
			}
			while (unfinished) {
				//Find a random sequence within bin and mutate it
				for (int i=0; i<nDivs; i++) {
					if (binCounts[i]>maxBinCount || binCounts[i]==0) continue;
					tempSeq	= new ArrayList<Long>(probeSet[i]);
					testSeqs= mutator(tempSeq.get(generator.nextInt(tempSeq.size())), k);
					//Dispatch & Bin
					weight	= so.seqEval(testSeqs);
					for (int j=0; j<testSeqs.length; j++) {
						currBin	= (int) (Math.floor((weight[j])/divSize));
						probeSet[currBin].add(testSeqs[j]);	
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
				for (int i=0; i<100000; i++) {
					currBin	= (int) (Math.floor((weight[i])/divSize));
					probeSet[currBin].add(testSeqs[i]);					
				}

				//Test counts for termination
				unfinished = false;
				for (int i=0; i<nDivs; i++) {
					binCounts[i] = probeSet[i].size();
					if (binCounts[i]>1E6) {
						binCounts[i]=0;
						probeSet[i] = new HashSet<Long>();
					}
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
			
			//Find at least the first nMotifs
			for (int i=nDivs-1; i>=0; i--) {
				if (binCounts[i]==0) continue;
				probes.addAll(probeSet[i]);
				if (probes.size() > nMotifs) {
					break;
				}
			}
			
			//Rescore and re-order
			weight = new double[probes.size()];
			for (int i=0; i<probes.size(); i++) {
				weight[i] = so.seqEval(probes.get(i));
			}
			
			//Build top nMotif list
			testSeqs = new long[nMotifs];
			outWeight= new double[nMotifs];
			for (int i=0; i<nMotifs; i++) {
				currBin			= maxValue(weight);
				testSeqs[i]		= probes.get(currBin);
				outWeight[i]	= weight[currBin];
				weight[currBin]	= Double.MIN_VALUE;
			}
			if (testSeqs[0]==topSeq) {
				unconverged = false;
			} else {
				topSeq	= testSeqs[0];
				bestSeq	= testSeqs[0];
			}
		}
		
		//Print
		Sequence seq, rcSeq;
		if (save) {
			try {
				outputFile = new PrintStream(new FileOutputStream(savePath));
				System.setOut(outputFile);	
				System.out.println("Forward,Reverse,Weight");
				for (int i=0; i<nMotifs; i++) {
					seq		= new Sequence(testSeqs[i], k);
					rcSeq	= new Sequence(seq.reverseComplement(), k);
					System.out.println(seq.getString()+","+rcSeq.getString()+","
							+ outWeight[i]);
				}
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		} else {
			System.out.println("Forward\tReverse\tWeight");
			for (int i=0; i<nMotifs; i++) {
				seq		= new Sequence(testSeqs[i], k);
				rcSeq	= new Sequence(seq.reverseComplement(), k);
				System.out.println(seq.getString()+"\t"+rcSeq.getString()+"\t"+
						outWeight[i]);
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
	
	private static int maxValue(double[] inputArray) {
		int bestIdx		= 0;
		double currMax	= inputArray[0];
		
		for (int i=0; i<inputArray.length; i++) {
			if (inputArray[i]>currMax) {
				currMax = inputArray[i];
				bestIdx = i;
			}
		}
		return bestIdx;
	}
}
