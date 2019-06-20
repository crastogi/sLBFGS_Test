package utils;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;

import base.Shape;
import model.*;

public class GenomeScorer {
	public static String dir		= "/Users/chaitanya/Documents/Research/SELEX/Biological Analysis/";
	public static String R0Dir		= "./Round0/";
	public static String shapeDir	= "./../Shape/";
	public static String results	= dir+"/HM Dataset/Fits/AbdB_Monomer.dat";
	public static int modelIndex	= 23;
	public static int modeIndex		= 0;
	public static String chrDir		= dir+"/SVB Enhancer/Genomes/dm3/";
	public static String outDir		= dir+"/SVB Enhancer/scoring/dm3/AbdBMonomer_Di";
	
	public static void main(String[] args) {
		int totBases, k;
		long currSeq;
		FA chr;
		Shape shapeModel;
		Round0Model R0;
		ScoringObject so;
		DataOutputStream osF, osR, fs;
		MultinomialResults mr	= new MultinomialResults(results);
		int[] status			= new int[1];
		String[] parsedPath;
		File[] chrPaths			= (new File(chrDir)).listFiles();
		
		//Load Round0 model from results file
		parsedPath	= mr.r0ModelPath.split("/");
		R0			= new Round0Model(R0Dir+parsedPath[parsedPath.length-1], mr.l, 
						mr.r0k, mr.r0Flank, mr.lFlank, mr.rFlank);
		//Create new sliding window system
		if (mr.getFit(modelIndex).isShape) {
			shapeModel = new Shape(shapeDir, true);
			shapeModel.setFeatures(mr.getFit(modelIndex).shapes);
		} else {
			shapeModel = null;
		}
		so = new ScoringObject(mr, modelIndex, modeIndex, true, false, shapeModel, R0);
		if (mr.getFit(modelIndex).isMulti) {
			k = mr.getFit(modelIndex).ks[modeIndex];
		} else {
			k = mr.getFit(modelIndex).k;
		}
		
		double tStart = System.currentTimeMillis();
		for (File currChr : chrPaths) {
			if (currChr.isFile()) {
				chr = new FA(chrDir+currChr.getName());
				if (!chr.name.equals("chrX")) {
					continue;
				}
				System.out.println("Currently Processing "+chr.name);
				//Create output file
				try {
					osF = new DataOutputStream(new FileOutputStream(outDir+"_"+chr.name+"_F.dat"));
					osR = new DataOutputStream(new FileOutputStream(outDir+"_"+chr.name+"_R.dat"));
					fs = new DataOutputStream(new FileOutputStream(outDir+"_"+chr.name+"_Flag.dat"));
					//Evaluate chromosome
					totBases = chr.nBases;
					for (int start=0; start<totBases-k+1; start++) {
						currSeq = chr.getLong(start, start+k-1, status);
						if (status[0]==-1) {
							osF.writeDouble(0.0);
							osR.writeDouble(0.0);
							fs.writeByte(-1);
						} else if (status[0]==1) {
							osF.writeDouble(so.seqEval(currSeq));
							osR.writeDouble(so.seqEval(reverseComplement(currSeq, k)));
							fs.writeByte(1);
						} else {
							osF.writeDouble(so.seqEval(currSeq));
							osR.writeDouble(so.seqEval(reverseComplement(currSeq, k)));
							fs.writeByte(0);
						}
					}
					osF.close();
					fs.close();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		System.out.println(System.currentTimeMillis()-tStart);
	}
	
	static long reverseComplement(long input, int length) {
		long output = 0;
		input = ~input;
		output = output | (input & 3);
		for (int i=1; i<length; i++) {
			input = input >> 2;
			output = output << 2;
			output = output | (input & 3);
		}
		return output;
	}
}
