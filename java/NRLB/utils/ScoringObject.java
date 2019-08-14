package utils;

import java.util.Arrays;

import base.*;
import model.*;

public class ScoringObject {
	private boolean singleScore, sumScore;
	private int nModes, totFeatures, totFrames, type, nShapeClasses, shapeWindows;
	private double nsBindingValue;
	private SingleModeSW sEval;
	private int[] ks;
	private SingleModeSW[] sw;
	
	public ScoringObject(MultinomialResults results, int index, 
			Shape shapeModel, boolean sumScore, Round0Model R0) {
		this(results, index, -1, false, sumScore, shapeModel, R0);
	}
	
	public ScoringObject(MultinomialResults results, int index, 
			boolean singleScore, boolean sumScore, Shape shapeModel, Round0Model R0) {
		this(results, index, -1, singleScore, sumScore, shapeModel, R0);
	}
	
	public ScoringObject(MultinomialResults results, int index, int mode, 
			boolean singleScore, boolean sumScore, Shape shapeModel, Round0Model R0) {
		boolean isFlank;
		int l, flankLength;
		double[] betas;
		int[][][] betaIdx;
		MultinomialFit fit	= results.getFit(index);
		Shape currShapeModel	= null;
	
		this.singleScore	= singleScore;
		this.sumScore		= (singleScore) ? sumScore : false;
		if (fit.isShape) {
			currShapeModel = shapeModel;
			currShapeModel.setFeatures(fit.shapes);
		}
		if (fit.isMulti) {				//Multi-mode model
			if (singleScore) {			//Score a single window
				if (mode==-1) {			//No specific mode has been selected
					throw new IllegalArgumentException("Must provide a mode "
							+ "index for a MultiModeModel!");
				}
				l			= fit.ks[mode];
				ks			= new int[]{fit.ks[mode]};
				nModes		= 1;
				isFlank		= false;
				flankLength	= 0;
			} else {					//Score a probe
				l			= results.l;
				if (mode==-1) {			//No Specific Mode has been selected
					ks		= fit.ks;
					nModes	= ks.length;
				} else {				//Specific mode has been selected
					ks		= new int[]{fit.ks[mode]};
					nModes	= 1;
				}
				isFlank		= fit.isFlank;
				flankLength	= fit.flankLength;
			}
			sw			= new SingleModeSW[nModes];
			betas		= fit.finalPosition;
			betaIdx		= fit.betaIdx;
			totFrames	= 0;
			totFeatures	= 0;			
			for (int i=0; i<nModes; i++) {
				sw[i] = new SingleModeSW(currShapeModel, R0, l, ks[i], 
						isFlank, flankLength, results.lFlank, 
						results.rFlank, fit.isDinuc, fit.isShape);
				sw[i].setOffsets(totFeatures, totFrames);
				totFrames	+= 2*sw[i].maxFrames;
				totFeatures	+= sw[i].nFeatures;
				if (mode!=-1) {
					sw[i].setParams(Arrays.copyOfRange(betas, betaIdx[mode][0][0], 
							betaIdx[mode][0][1]+1));					
				} else {
					sw[i].setParams(Arrays.copyOfRange(betas, betaIdx[i][0][0], 
							betaIdx[i][0][1]+1));										
				}
				sw[i].reverseBetas();	//Keep betas reversed
			}		
		} else {						//Multinomial model
			l			= (singleScore) ? fit.k : results.l;
			ks			= new int[]{fit.k};
			isFlank		= (singleScore) ? false : fit.isFlank;
			flankLength	= (singleScore) ? 0 : fit.flankLength;
			nModes		= 1;
			sw			= new SingleModeSW[nModes];
			betas		= fit.finalPosition;
			sw[0] 		= new SingleModeSW(currShapeModel, R0, l, fit.k, 
						isFlank, flankLength, results.lFlank, 
						results.rFlank, fit.isDinuc, fit.isShape);
			sw[0].setOffsets(0, 0);
			totFrames	= 2*sw[0].maxFrames;
			totFeatures	= sw[0].nFeatures;
			sw[0].setParams(betas);
			sw[0].reverseBetas();			//Keep betas reversed			
		}
		if (fit.isNSBinding && !singleScore) {
			nsBindingValue = Math.exp(betas[betas.length-1]);
			totFeatures++;
		}
		//What type of sliding window should be used?
		if (fit.isFlank) {
			if (fit.isDinuc) {
				if (fit.isShape) {
					type = 3;		//Nuc+Dinuc+Shape
				} else {
					type = 1;		//Nuc+Dinuc
				}
			} else {
				if (fit.isShape) {
					type = 2;		//Nuc+Shape
				} else {
					type = 0;		//Nuc
				}
			}			
		} else {
			if (fit.isDinuc) {
				if (fit.isShape) {
					type = 7;		//Nuc+Dinuc+Shape
				} else {
					type = 5;		//Nuc+Dinuc
				}
			} else {
				if (fit.isShape) {
					type = 6;		//Nuc+Shape
				} else {
					type = 4;		//Nuc
				}
			}
		}
	}
	
	public ScoringObject(MultinomialResults results, int index, int mode, 
			Shape shapeModel, Round0Model R0, Shape evalShapeModel) {
		boolean isFlank;
		int l, flankLength;
		double[] betas;
		int[][][] betaIdx;
		MultinomialFit fit	= results.getFit(index);
				
		if (fit.isMulti) {				//Multi-mode model
			l			= results.l;
			if (mode==-1) {				//No Specific Mode has been selected
				ks		= fit.ks;
				nModes	= ks.length;
			} else {					//Specific mode has been selected
				ks		= new int[]{fit.ks[mode]};
				nModes	= 1;
			}
			isFlank		= fit.isFlank;
			flankLength	= fit.flankLength;
			sw			= new SingleModeSW[nModes];
			betas		= fit.finalPosition;
			betaIdx		= fit.betaIdx;
			totFrames	= 0;
			totFeatures	= 0;			
			for (int i=0; i<nModes; i++) {
				sw[i]		= new SingleModeSW(shapeModel, R0, l, ks[i], 
							isFlank, flankLength, results.lFlank, 
							results.rFlank, fit.isDinuc, fit.isShape);
				sw[i].setOffsets(totFeatures, totFrames);
				totFrames	+= 2*sw[i].maxFrames;
				totFeatures	+= sw[i].nFeatures;
				sw[i].setParams(Arrays.copyOfRange(betas, betaIdx[i][0][0], 
						betaIdx[i][0][1]+1));
				sw[i].reverseBetas();	//Keep betas reversed
			}
			//Set up the shape evaluator
			sEval		= new SingleModeSW(evalShapeModel, R0, l, fit.k, 
						isFlank, flankLength, results.lFlank, results.rFlank, 
						fit.isDinuc, fit.isShape);
			shapeWindows= (isFlank) ? l+2*flankLength : l;
		} else {						//Multinomial model
			l			= results.l;
			ks			= new int[]{fit.k};
			isFlank		= fit.isFlank;
			flankLength	= fit.flankLength;
			nModes		= 1;
			sw			= new SingleModeSW[nModes];
			betas		= fit.finalPosition;
			sw[0] 		= new SingleModeSW(shapeModel, R0, l, fit.k, 
						isFlank, flankLength, results.lFlank, 
						results.rFlank, fit.isDinuc, fit.isShape);
			sw[0].setOffsets(0, 0);
			totFrames	= 2*sw[0].maxFrames;
			totFeatures	= sw[0].nFeatures;
			sw[0].setParams(betas);
			sw[0].reverseBetas();			//Keep betas reversed		
			//Set up the shape evaluator
			sEval		= new SingleModeSW(evalShapeModel, R0, l, fit.k, 
						isFlank, flankLength, results.lFlank, results.rFlank, 
						fit.isDinuc, fit.isShape);
			shapeWindows= (isFlank) ? l+2*flankLength-4 : l-4;
		}
		if (fit.isNSBinding) {
			nsBindingValue = Math.exp(betas[betas.length-1]);
			totFeatures++;
		}
		nShapeClasses	= evalShapeModel.nShapeFeatures();
		//What type of sliding window should be used?
		if (fit.isFlank) {
			if (fit.isDinuc) {
				if (fit.isShape) {
					type = 3;		//Nuc+Dinuc+Shape
				} else {
					type = 1;		//Nuc+Dinuc
				}
			} else {
				if (fit.isShape) {
					type = 2;		//Nuc+Shape
				} else {
					type = 0;		//Nuc
				}
			}			
		} else {
			if (fit.isDinuc) {
				if (fit.isShape) {
					type = 7;		//Nuc+Dinuc+Shape
				} else {
					type = 5;		//Nuc+Dinuc
				}
			} else {
				if (fit.isShape) {
					type = 6;		//Nuc+Shape
				} else {
					type = 4;		//Nuc
				}
			}
		}
	}
	
	public int getK() {
		if (sw.length > 1) {
			throw new IllegalArgumentException("No MultiMode mode selected!");
		}
		return sw[0].k;
	}
	
	public double getZ() {
		double totalZ = 0;
		
		for (int i=0; i<nModes; i++) {
			totalZ += sw[i].getZ();
		}
		totalZ /= sw[0].R0Model.getZ();
		return (totalZ + nsBindingValue);
	}
	
	public double[] seqEval(long[] input) {
		double sum;
		double[] kappas = new double[totFrames];
		double[] output = new double[input.length];
		
		if (singleScore) {
			if (type==0) {
				for (int k=0; k<input.length; k++) {
					for (int j=0; j<nModes; j++) {
						sw[j].swNucleotide(input[k], kappas);					
					}
					output[k] = (sumScore) ? kappas[0]+kappas[1] : kappas[0];
				}
			} else if(type==1) {
				for (int k=0; k<input.length; k++) {
					for (int j=0; j<nModes; j++) {
						sw[j].swNucleotideDinucleotide(input[k], kappas);					
					}
					output[k] = (sumScore) ? kappas[0]+kappas[1] : kappas[0];
				}
			} else if(type==2) {
				for (int k=0; k<input.length; k++) {
					for (int j=0; j<nModes; j++) {
						sw[j].swNucleotideShape(input[k], kappas);					
					}
					output[k] = (sumScore) ? kappas[0]+kappas[1] : kappas[0];
				}
			} else if(type==3){
				for (int k=0; k<input.length; k++) {
					for (int j=0; j<nModes; j++) {
						sw[j].swNucleotideDinucleotideShape(input[k], kappas);					
					}
					output[k] = (sumScore) ? kappas[0]+kappas[1] : kappas[0];
				}
			} else if(type==4){
				for (int k=0; k<input.length; k++) {
					for (int j=0; j<nModes; j++) {
						sw[j].swNucleotideNoFlank(input[k], kappas);					
					}
					output[k] = (sumScore) ? kappas[0]+kappas[1] : kappas[0];
				}
			} else if(type==5) {
				for (int k=0; k<input.length; k++) {
					for (int j=0; j<nModes; j++) {
						sw[j].swNucleotideDinucleotideNoFlank(input[k], kappas);					
					}
					output[k] = (sumScore) ? kappas[0]+kappas[1] : kappas[0];
				}
			} else if(type==6) {
				for (int k=0; k<input.length; k++) {
					for (int j=0; j<nModes; j++) {
						sw[j].swNucleotideShapeNoFlank(input[k], kappas);					
					}
					output[k] = (sumScore) ? kappas[0]+kappas[1] : kappas[0];
				}
			} else if(type==7) {
				for (int k=0; k<input.length; k++) {
					for (int j=0; j<nModes; j++) {
						sw[j].swNucleotideDinucleotideShapeNoFlank(input[k], kappas);					
					}
					output[k] = (sumScore) ? kappas[0]+kappas[1] : kappas[0];
				}
			}
		} else {
			if (type==0) {
				for (int k=0; k<input.length; k++) {
					for (int j=0; j<nModes; j++) {
						sw[j].swNucleotide(input[k], kappas);					
					}
					sum = nsBindingValue;
					for (int j=0; j<totFrames; j++) {
						sum += kappas[j];
					}
					output[k] = sum;
				}
			} else if(type==1) {
				for (int k=0; k<input.length; k++) {
					for (int j=0; j<nModes; j++) {
						sw[j].swNucleotideDinucleotide(input[k], kappas);					
					}
					sum = nsBindingValue;
					for (int j=0; j<totFrames; j++) {
						sum += kappas[j];
					}
					output[k] = sum;
				}
			} else if(type==2) {
				for (int k=0; k<input.length; k++) {
					for (int j=0; j<nModes; j++) {
						sw[j].swNucleotideShape(input[k], kappas);					
					}
					sum = nsBindingValue;
					for (int j=0; j<totFrames; j++) {
						sum += kappas[j];
					}
					output[k] = sum;
				}
			} else if(type==3){
				for (int k=0; k<input.length; k++) {
					for (int j=0; j<nModes; j++) {
						sw[j].swNucleotideDinucleotideShape(input[k], kappas);					
					}
					sum = nsBindingValue;
					for (int j=0; j<totFrames; j++) {
						sum += kappas[j];
					}
					output[k] = sum;
				}
			} else if(type==4){
				for (int k=0; k<input.length; k++) {
					for (int j=0; j<nModes; j++) {
						sw[j].swNucleotideNoFlank(input[k], kappas);					
					}
					sum = nsBindingValue;
					for (int j=0; j<totFrames; j++) {
						sum += kappas[j];
					}
					output[k] = sum;
				}
			} else if(type==5) {
				for (int k=0; k<input.length; k++) {
					for (int j=0; j<nModes; j++) {
						sw[j].swNucleotideDinucleotideNoFlank(input[k], kappas);					
					}
					sum = nsBindingValue;
					for (int j=0; j<totFrames; j++) {
						sum += kappas[j];
					}
					output[k] = sum;
				}
			} else if(type==6) {
				for (int k=0; k<input.length; k++) {
					for (int j=0; j<nModes; j++) {
						sw[j].swNucleotideShapeNoFlank(input[k], kappas);					
					}
					sum = nsBindingValue;
					for (int j=0; j<totFrames; j++) {
						sum += kappas[j];
					}
					output[k] = sum;
				}
			} else if(type==7) {
				for (int k=0; k<input.length; k++) {
					for (int j=0; j<nModes; j++) {
						sw[j].swNucleotideDinucleotideShapeNoFlank(input[k], kappas);					
					}
					sum = nsBindingValue;
					for (int j=0; j<totFrames; j++) {
						sum += kappas[j];
					}
					output[k] = sum;
				}
			}
		}
		return output;
	}
	
	public double seqEval(long input) {
		double sum = nsBindingValue;
		double[] kappas = new double[totFrames];
		
		if (singleScore) {
			if (type==0) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotide(input, kappas);
				}
				return ((sumScore) ? kappas[0]+kappas[1] : kappas[0]);
			} else if(type==1) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideDinucleotide(input, kappas);
				}
				return ((sumScore) ? kappas[0]+kappas[1] : kappas[0]);
			} else if(type==2) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideShape(input, kappas);
				}
				return ((sumScore) ? kappas[0]+kappas[1] : kappas[0]);
			} else if(type==3){
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideDinucleotideShape(input, kappas);
				}
				return ((sumScore) ? kappas[0]+kappas[1] : kappas[0]);
			} else if(type==4){
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideNoFlank(input, kappas);
				}
				return ((sumScore) ? kappas[0]+kappas[1] : kappas[0]);
			} else if(type==5) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideDinucleotideNoFlank(input, kappas);
				}
				return ((sumScore) ? kappas[0]+kappas[1] : kappas[0]);
			} else if(type==6) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideShapeNoFlank(input, kappas);
				}
				return ((sumScore) ? kappas[0]+kappas[1] : kappas[0]);
			} else if(type==7) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideDinucleotideShapeNoFlank(input, kappas);
				}
				return ((sumScore) ? kappas[0]+kappas[1] : kappas[0]);
			}
		} else {
			if (type==0) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotide(input, kappas);
				}
				for (int j=0; j<totFrames; j++) {
					sum += kappas[j];
				}
			} else if(type==1) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideDinucleotide(input, kappas);
				}
				for (int j=0; j<totFrames; j++) {
					sum += kappas[j];
				}
			} else if(type==2) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideShape(input, kappas);
				}
				for (int j=0; j<totFrames; j++) {
					sum += kappas[j];
				}
			} else if(type==3){
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideDinucleotideShape(input, kappas);
				}
				for (int j=0; j<totFrames; j++) {
					sum += kappas[j];
				}
			} else if(type==4){
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideNoFlank(input, kappas);
				}
				for (int j=0; j<totFrames; j++) {
					sum += kappas[j];
				}
			} else if(type==5) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideDinucleotideNoFlank(input, kappas);
				}
				for (int j=0; j<totFrames; j++) {
					sum += kappas[j];
				}
			} else if(type==6) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideShapeNoFlank(input, kappas);
				}
				for (int j=0; j<totFrames; j++) {
					sum += kappas[j];
				}
			} else if(type==7) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideDinucleotideShapeNoFlank(input, kappas);
				}
				for (int j=0; j<totFrames; j++) {
					sum += kappas[j];
				}
			}
		}
		return sum;
	}
	
	public double[] perWindowSeqEval(long input) {
		double[] kappas = new double[totFrames];
		
		if (type==0) {
			for (int j=0; j<nModes; j++) {
				sw[j].swNucleotide(input, kappas);
			}
		} else if(type==1) {
			for (int j=0; j<nModes; j++) {
				sw[j].swNucleotideDinucleotide(input, kappas);
			}
		} else if(type==2) {
			for (int j=0; j<nModes; j++) {
				sw[j].swNucleotideShape(input, kappas);
			}
		} else if(type==3){
			for (int j=0; j<nModes; j++) {
				sw[j].swNucleotideDinucleotideShape(input, kappas);
			}
		} else if(type==4){
			for (int j=0; j<nModes; j++) {
				sw[j].swNucleotideNoFlank(input, kappas);
			}
		} else if(type==5) {
			for (int j=0; j<nModes; j++) {
				sw[j].swNucleotideDinucleotideNoFlank(input, kappas);
			}
		} else if(type==6) {
			for (int j=0; j<nModes; j++) {
				sw[j].swNucleotideShapeNoFlank(input, kappas);
			}
		} else if(type==7) {
			for (int j=0; j<nModes; j++) {
				sw[j].swNucleotideDinucleotideShapeNoFlank(input, kappas);
			}
		}
		return kappas;
	}
	
	public double[] gradEval(long input) {
		double[] kappas		= new double[totFrames];
		double[] gradients	= new double[totFeatures];
		double[] gradVar	= new double[totFeatures];
		
		if (type==0) {
			for (int j=0; j<nModes; j++) {
				sw[j].swNucleotide(input, kappas);
			}
			for (int j=0; j<nModes; j++) {						
				sw[j].swGradNucleotide(input, 1, 1, kappas, gradients, gradVar);
			}
		} else if(type==1) {
			for (int j=0; j<nModes; j++) {
				sw[j].swNucleotideDinucleotide(input, kappas);
			}
			for (int j=0; j<nModes; j++) {						
				sw[j].swGradNucleotideDinucleotide(input, 1, 1, kappas, gradients, gradVar);
			}
		} else if(type==2) {
			for (int j=0; j<nModes; j++) {
				sw[j].swNucleotideShape(input, kappas);
			}
			for (int j=0; j<nModes; j++) {						
				sw[j].swGradNucleotideShape(input, 1, 1, kappas, gradients, gradVar);
			}
		} else if(type==3){
			for (int j=0; j<nModes; j++) {
				sw[j].swNucleotideDinucleotideShape(input, kappas);
			}
			for (int j=0; j<nModes; j++) {						
				sw[j].swGradNucleotideDinucleotideShape(input, 1, 1, kappas, gradients, gradVar);
			}
		} else if(type==4){
			for (int j=0; j<nModes; j++) {
				sw[j].swNucleotideNoFlank(input, kappas);
			}
			for (int j=0; j<nModes; j++) {						
				sw[j].swGradNucleotideNoFlank(input, 1, 1, kappas, gradients, gradVar);
			}
		} else if(type==5) {
			for (int j=0; j<nModes; j++) {
				sw[j].swNucleotideDinucleotideNoFlank(input, kappas);
			}
			for (int j=0; j<nModes; j++) {						
				sw[j].swGradNucleotideDinucleotideNoFlank(input, 1, 1, kappas, gradients, gradVar);
			}
		} else if(type==6) {
			for (int j=0; j<nModes; j++) {
				sw[j].swNucleotideShapeNoFlank(input, kappas);
			}
			for (int j=0; j<nModes; j++) {						
				sw[j].swGradNucleotideShapeNoFlank(input, 1, 1, kappas, gradients, gradVar);
			}
		} else if(type==7) {
			for (int j=0; j<nModes; j++) {
				sw[j].swNucleotideDinucleotideShapeNoFlank(input, kappas);
			}
			for (int j=0; j<nModes; j++) {						
				sw[j].swGradNucleotideDinucleotideShapeNoFlank(input, 1, 1, kappas, gradients, gradVar);
			}
		}
		return gradients;
	}
	
	public Object[] shapeEval(long[] input) {
		int maxIdx;
		int[] position		= new int[input.length];
		double[] kappas 	= new double[totFrames];
		double[] weight 	= new double[input.length];
		double[][][] profile= new double[input.length][nShapeClasses][shapeWindows];
		
		if (type==0) {
			for (int k=0; k<input.length; k++) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotide(input[k], kappas);					
				}
				maxIdx		= Array.whichMax(kappas);
				weight[k]	= kappas[maxIdx];
				position[k]	= (int) maxIdx/2;
				profile[k]	= sEval.swShapeProfile(input[k]);
			}
		} else if(type==1) {
			for (int k=0; k<input.length; k++) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideDinucleotide(input[k], kappas);					
				}
				maxIdx		= Array.whichMax(kappas);
				weight[k]	= kappas[maxIdx];
				position[k]	= (int) maxIdx/2;
				profile[k]	= sEval.swShapeProfile(input[k]);
			}
		} else if(type==2) {
			for (int k=0; k<input.length; k++) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideShape(input[k], kappas);					
				}
				maxIdx		= Array.whichMax(kappas);
				weight[k]	= kappas[maxIdx];
				position[k]	= (int) maxIdx/2;
				profile[k]	= sEval.swShapeProfile(input[k]);
			}
		} else if(type==3){
			for (int k=0; k<input.length; k++) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideDinucleotideShape(input[k], kappas);					
				}
				maxIdx		= Array.whichMax(kappas);
				weight[k]	= kappas[maxIdx];
				position[k]	= (int) maxIdx/2;
				profile[k]	= sEval.swShapeProfile(input[k]);
			}
		} else if(type==4){
			for (int k=0; k<input.length; k++) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideNoFlank(input[k], kappas);					
				}
				maxIdx		= Array.whichMax(kappas);
				weight[k]	= kappas[maxIdx];
				position[k]	= (int) maxIdx/2;
				profile[k]	= sEval.swShapeProfile(input[k]);
			}
		} else if(type==5) {
			for (int k=0; k<input.length; k++) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideDinucleotideNoFlank(input[k], kappas);					
				}
				maxIdx		= Array.whichMax(kappas);
				weight[k]	= kappas[maxIdx];
				position[k]	= (int) maxIdx/2;
				profile[k]	= sEval.swShapeProfile(input[k]);
			}
		} else if(type==6) {
			for (int k=0; k<input.length; k++) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideShapeNoFlank(input[k], kappas);					
				}
				maxIdx		= Array.whichMax(kappas);
				weight[k]	= kappas[maxIdx];
				position[k]	= (int) maxIdx/2;
				profile[k]	= sEval.swShapeProfile(input[k]);
			}
		} else if(type==7) {
			for (int k=0; k<input.length; k++) {
				for (int j=0; j<nModes; j++) {
					sw[j].swNucleotideDinucleotideShapeNoFlank(input[k], kappas);					
				}
				maxIdx		= Array.whichMax(kappas);
				weight[k]	= kappas[maxIdx];
				position[k]	= (int) maxIdx/2;
				profile[k]	= sEval.swShapeProfile(input[k]);
			}
		}
		return new Object[]{weight, position, profile};
	}
}
