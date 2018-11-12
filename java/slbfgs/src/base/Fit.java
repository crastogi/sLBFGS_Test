package base;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

public abstract class Fit {
	public int fitSteps, functionCalls;
	public double fitTime, likelihood;
	public double[] finalPosition				= null;
	public ArrayList<double[]> trajectory		= new ArrayList<double[]>();

	public abstract void recordFit(int fitSteps, int functionCalls, 
			double fitTime, double likelihood, Model input); 
	
	public double[] positionVector() {
		return Array.clone(finalPosition);
	}

	public void addStep(double[] positionVector) {
		trajectory.add(positionVector);
	}
	
	public void printTrajectories() {
		for (int i=0; i<trajectory.size(); i++) {
			Array.print(trajectory.get(i));
		}
	}
	
	public void printTrajectories(String filePath) {
		try {

			PrintStream outputFile	= new PrintStream(
					new FileOutputStream(filePath, false));
			PrintStream original	= System.out;
			
			System.setOut(outputFile);
			printTrajectories();
			System.setOut(original);
		} catch (FileNotFoundException e) {
			System.out.println("Cannot create trajectory file at this "
					+ "location: "+filePath);
			e.printStackTrace();
		}
	}
	
	public void printTrajectories(String filePath, boolean append) {
		if(append == false || trajectory.size() <= 10 ) //Rewrites the full file the first couple of iterations.
			printTrajectories(filePath);
		else {
			try {
				PrintStream outputFile	= new PrintStream(
						new FileOutputStream(filePath, true));
				PrintStream original	= System.out;
				
				System.setOut(outputFile);
				Array.print(trajectory.get(trajectory.size()-1));
				System.setOut(original);
				
			} catch (FileNotFoundException e) {
				System.out.println("Cannot create trajectory file at this "
						+ "location: "+filePath);
				e.printStackTrace();
			}
		}
	}
	
}