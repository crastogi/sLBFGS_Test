package test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

import base.Array;

public class GenericTestFile {
	
	public static void main(String[] args) {
		int maxCols;
        String csvFile = "./src/featurized_df.csv";
        String line;
        double[] temp;
        String[] split;
        double[][] data;
        ArrayList<double[]> parsed = new ArrayList<double[]>();

        try (BufferedReader br = new BufferedReader(new FileReader(csvFile))) {
        	// Ignore first line
        	line = br.readLine();
        	maxCols = line.split(",").length;
            // Read, split, and ensure that the remaining lines have the same number of columns
        	// Store the remaining values in a matrix
        	while ((line = br.readLine()) != null) {
            	split = line.split(",");
            	if (split.length != maxCols) {
            		throw new Exception("Variable number of columns!");
            	}
            	temp = new double[maxCols];
            	for (int i=0; i<maxCols; i++) {
            		temp[i] = Double.parseDouble(split[i]);
            	}
            	parsed.add(temp);
            }
        	// Store values in data
        	data = new double[parsed.size()][maxCols];
        	for (int i=0; i<parsed.size(); i++) {
        		for (int j=0; j<maxCols; j++) {
        			data[i][j] = parsed.get(i)[j];
        		}
        	}
        	Array.print(data);
        } catch (Exception e) {
            e.printStackTrace();
        }
	}
}
