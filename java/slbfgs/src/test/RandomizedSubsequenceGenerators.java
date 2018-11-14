package test;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

import base.MersenneTwisterFast;

public class RandomizedSubsequenceGenerators {
	static Random generator;
	static MersenneTwisterFast mtfast;
	
	public static void main(String[] args) {
		int testCases = 1*1000;
		int N = 1000*1000;
		int k = 1000;
		int shuffleSamples;
		double tStart, tEnd;
		generator = new Random();
		mtfast = new MersenneTwisterFast();
		int[] testOutput, avgFreq;
		
		// Test Reservoir sampling
		tStart = System.nanoTime();	
		for (int i=0; i<testCases; i++) {
			ReservoirSample(N, k);
		}
		tEnd = System.nanoTime();
		avgFreq = new int[N];
		for (int i=0; i<testCases; i++) {
			testOutput = ReservoirSample(N, k);
			for (int j=0; j<k; j++) {
				avgFreq[testOutput[j]]++;
			}
		}
		System.out.println("Reservoir Sampling");
		System.out.println((tEnd-tStart)/(1E9*testCases)+"\t"+sd(avgFreq, k));
//		Array.print(avgFreq);
		
		// Test Simple sampling
		tStart = System.nanoTime();	
		for (int i=0; i<testCases; i++) {
			Simple(N, k);
		}
		tEnd = System.nanoTime();
		avgFreq = new int[N];
		for (int i=0; i<testCases; i++) {
			testOutput = Simple(N, k);
			for (int j=0; j<k; j++) {
				avgFreq[testOutput[j]]++;
			}
		}
		System.out.println("Simple Sampling");
		System.out.println((tEnd-tStart)/(1E9*testCases)+"\t"+sd(avgFreq, k));
//		Array.print(avgFreq);
		
		// Test Fisher Yates. Note that this produces a full shuffle, so we have to give it the benefit of doing that
		shuffleSamples = testCases*k/N;
		tStart = System.nanoTime();
		for (int i=0; i<shuffleSamples; i++) {
			FisherYates(N);
		}
		tEnd = System.nanoTime();
		avgFreq = new int[N];
		for (int i=0; i<shuffleSamples; i++) {
			testOutput = FisherYates(N);
			for (int j=0; j<N; j++) {
				avgFreq[testOutput[j]]++;
			}
		}
		System.out.println("Fisher-Yates Shuffle");
		System.out.println((tEnd-tStart)/(1E9*testCases)+"\t"+sd(avgFreq, k));
	}
	
	// Reservoir sampling
	public static int[] ReservoirSample(int N, int k) {
		int randIdx;
		int[] R = new int[k];
		
		// Fill array
		for (int i=0; i<k; i++) {
			R[i] = i;
		}
		
		// Replace elements with decreasing probability
		for (int i=k; i<N; i++) {
			randIdx = mtfast.nextInt(i+1); //tester.choose(0, i+1); //generator.nextInt(i+1);
			if (randIdx<k) {
				R[randIdx] = i;
			}
		}
		
		return R;
	}
	
	// Simple duplicate removal
	public static int[] Simple(int N, int k) {
		int idx = 0;
		int[] output = new int[k];
		HashSet<Integer> h = new HashSet<Integer>();
		
		while (true) {
			h.add(mtfast.nextInt(N)); //tester.choose(0, N)); //generator.nextInt(N));
			if (h.size() >= k) {
				break;
			}
		}
		Iterator<Integer> i = h.iterator(); 
        while (i.hasNext()) {
        	output[idx] = i.next();
        	idx++;
        }
		
		return output;
	}
	
	// Fisher-Yates shuffle
	public static int[] FisherYates(int N) {
		int idx, hold;
		int[] output = new int[N];
		
		// Populate array
		for (int i=0; i<N; i++) {
			output[i] = i;
		}
		
		// Shuffle
		for (int i=N-1; i>0; i--) {
			idx = mtfast.nextInt(i+1);
			// Swap at index
			hold = output[idx];
			output[idx] = output[i];
			output[i] = hold;
		}
		
		return output;
	}
	
	// Compute sample standard deviation
	public static double sd(int[] input, int N) {
		double mean = 0, sse = 0;
		
		// First compute mean
		for (int i=0; i<input.length; i++) {
			mean += input[i];
		}
		mean /= input.length;		
		// Next compute sum squared error
		for (int i=0; i<input.length; i++) {
			sse += (input[i]-mean)*(input[i]-mean);
		}
		
		return Math.sqrt(sse/(input.length-1));
	}
}
