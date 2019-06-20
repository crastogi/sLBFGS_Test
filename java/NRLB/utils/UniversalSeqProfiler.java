package utils;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Arrays;

import base.Array;
import base.Sequence;
import base.Shape;
import model.MultinomialFit;
import model.MultinomialResults;
import model.Round0Model;

public class UniversalSeqProfiler{
	public static String dir		= "C:/Users/Chaitanya/Documents/Columbia Work/Research/SELEX/Biological Analysis/";
	public static String R0Dir		= "./Round0/";
	public static String shapeDir	= "./../Shape/";
	public static String results	= dir+"HM Dataset/Fits/UbxIVa.dat";
	public static int modelIndex	= 23;
	public static int modeIndex		= -1;
	public static String outFile	= dir+"Ubx_Score.txt";
	public static int lBases		= 4;
	public static int testK			= 12;
	public static double filterLevel= 5E-5;
	
	/*
	 * Does NOT evaluate the NSBinding component of kappa_i OR in the error 
	 * propagation 
	 */
	
	public static void main(String[] args) {
		double weight;
		Shape shapeModel;
		Sequence currSeq;
		Round0Model R0;
		ScoringObject so;
		MultinomialResults mr	= new MultinomialResults(results);
		MultinomialFit fit		= mr.getFit(modelIndex);
		double[] nuc, dinuc;
		String[] parsedPath;
		PrintStream outputFile	= null;
		
		//Load Round0 model from results file
		parsedPath	= mr.r0ModelPath.split("/");
		R0			= new Round0Model(R0Dir+parsedPath[parsedPath.length-1], mr.l, 
				mr.r0k, mr.r0Flank, mr.lFlank, mr.rFlank);

		//Create new sliding window system
		fit.isFlank		= false;
		fit.flankLength	= 0;
		fit.isNSBinding	= false;
		//Fix the offset and betas
		nuc = Arrays.copyOfRange(fit.finalPosition, lBases*4, lBases*4+testK*4);
		if (fit.isDinuc) {
			dinuc = Arrays.copyOfRange(fit.finalPosition, fit.k*4+lBases*16, fit.k*4+lBases*16+(testK-1)*16);
		} else {
			dinuc = null;
		}
		if (fit.isShape) {
			shapeModel = new Shape(shapeDir, true);
			shapeModel.setFeatures(fit.shapes);
		} else {
			shapeModel = null;
		}
		fit.k			= testK;
		mr.l			= testK;
		fit.finalPosition = Array.cat(Array.cat(nuc, dinuc), 0);
		so = new ScoringObject(mr, modelIndex, modeIndex, false, false, shapeModel, R0);
		
		//Generate all possible MITOMI sequences and print to file
		try {
			outputFile = new PrintStream(new FileOutputStream(outFile));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		System.setOut(outputFile);
		for (long i=0; i<Math.pow(4, testK); i++) {
			currSeq = new Sequence(i, testK);
			weight	= so.seqEval(i);
			if (weight>filterLevel) {
				System.out.println(currSeq.getString()+","+weight);
			}
		}
	}
}
