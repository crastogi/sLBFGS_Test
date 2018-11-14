package base;

import java.util.Formatter;

public class Minimizer {
	
	public double GradientTime = 0; 
	
	protected double epsilon, randomSeedScale = .05;
	protected Model model;
	
	public Fit doMinimize(double[] seed, String trajectoryFile) 
			throws Exception {
		return null;
	}
	
	protected Model.CompactGradientOutput evaluate(double[] input) throws Exception {
		model.setParams(input);
		return model.gradientEval();	
	}
	
	protected Model.CompactGradientOutput stochasticEvaluate(double[] input) throws Exception {
		double tStart = System.nanoTime();
		
		model.setParams(input);
		
		GradientTime += (System.nanoTime()-tStart)/1E9;
		
		return model.stochasticGradientEval();
	}
	
	protected void printStep(int iterations, int calls, double likelihood, 
			double distance, double ... params) {
		Formatter fmt = new Formatter();
		
		System.out.printf("   %7d      %7d   ", iterations, calls);
		fmt.format("%18.18s   %18.18s", String.format("%10.15f", likelihood), 
				String.format("%10.15f", distance));
		System.out.print(fmt);
		fmt = new Formatter();
		for (int i=0; i<params.length; i++) {
			fmt.format("   %18.18s", String.format("%10.15f", params[i]));
		}
		System.out.print(fmt+"\n");
		fmt.close();
	}
}
