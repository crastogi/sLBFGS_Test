package utils;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

import base.Array;
import base.Sequence;
import base.Shape;
import model.MultinomialFit;
import model.MultinomialResults;
import model.Round0Model;

public class MITOMIScorer {
	public static String dir		= "/Users/chaitanya/Documents/Research/SELEX/Biological Analysis/";
	public static String R0Dir		= "./Round0/";
	public static String shapeDir	= "./../Shape/";
	public static String results	= dir+"Max-MITOMI/Max.dat";
	public static int modelIndex	= 83;
	public static int modeIndex		= -1;
	public static String outFile	= dir+"MITOMI_83.txt";

	/*
	 * Does NOT evaluate the NSBinding component of kappa_i OR in the error 
	 * propagation 
	 */
	
	public static void main(String[] args) {
		boolean isError;
		double weight, totError;
		long lSeq, rSeq, scoreSeq;
		String lFlank			= "AATTG";
		String rFlank			= "GTGGGTGTCCGGCGGTATGAC";
		Shape shapeModel;
		Sequence currSeq;
		Round0Model R0;
		ScoringObject so;
		MultinomialResults mr	= new MultinomialResults(results);
		MultinomialFit fit		= mr.getFit(modelIndex);
		double[] errors			= null;
		double[] gradient;
		String[] parsedPath;
		PrintStream outputFile	= null;
		
		//Load Round0 model from results file
		parsedPath	= mr.r0ModelPath.split("/");
		R0			= new Round0Model(R0Dir+parsedPath[parsedPath.length-1], mr.l, 
				mr.r0k, mr.r0Flank, mr.lFlank, mr.rFlank);

		//Create new sliding window system
		mr.l			= 30;
		fit.isFlank		= false;
		fit.flankLength	= 0;
		isError			= fit.isError;
		if (isError) {
			errors		= Array.clone(fit.errors);
		}
		fit.isNSBinding	= false;
		if (fit.isShape) {
			shapeModel = new Shape(shapeDir, true);
			shapeModel.setFeatures(fit.shapes);
		} else {
			shapeModel = null;
		}
		so = new ScoringObject(mr, modelIndex, modeIndex, false, false, shapeModel, R0);
		
		//Generate all possible MITOMI sequences and print to file
		lSeq = (new Sequence(lFlank, 0, 5)).getValue();
		rSeq = (new Sequence(rFlank, 0, 21)).getValue();
		try {
			outputFile = new PrintStream(new FileOutputStream(outFile));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		System.setOut(outputFile);
		for (long l=0; l<(long) Math.pow(4, 4); l++) {
			currSeq = new Sequence(l, 4);
			scoreSeq= (((lSeq << 4*2) | currSeq.getValue()) << 21*2) | rSeq;
			weight	= so.seqEval(scoreSeq);
			if (isError) {
				gradient = so.gradEval(scoreSeq);
				totError = 0;
				for (int i=0; i<gradient.length; i++) {
					totError = gradient[i]*errors[i]*gradient[i]*errors[i];
				}
				totError = Math.sqrt(totError);
				System.out.println(currSeq.getString()+","+weight+","+totError);
			} else {
				System.out.println(currSeq.getString()+","+weight);
			}
		}
	}
}
