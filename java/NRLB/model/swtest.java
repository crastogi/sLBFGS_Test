package model;

import java.util.UUID;

import base.Array;
import base.Shape;

public class swtest {
	public static int k = 10;
	public static int l = 24;
	public static int flankLength = 3;
	public static int maxFrames = l-k+1+flankLength*2;
	public static long frameMask = (long) Math.pow(4, k) - 1;
	public static double[] nucBetas = Array.randomDouble(4*k);
	public static long rFlankingSequence = 0;
	public static long fFlankingSequence = 0;
	
	public static void main(String[] args) {
		if (true) {
			String test = "blah2 = 'poop'";
			String currLine, output = null;
			String[] parsed;
			
			System.out.println(test);
			currLine= test.replaceAll("\"", "");					//remove all quotes if they exist
			currLine= currLine.replaceAll("\'", "");
			System.out.println(currLine);
			parsed	= currLine.split("=");									//split the first and second halves of the entry
			if ((parsed[0].trim()).matches("(?i)"+"blah2")) {				//Perform case-insensitive matching
				output = (currLine.split("=")[1]).split("#|;")[0].trim();	//Remove trailing comments or semicolons
			}
			
			System.out.println(output);
			
			return;
		}
		
		
		MultinomialResults r = new MultinomialResults("/Users/chaitanya/Documents/Research/SELEX/Multinomial Paper/Figures/Fig3/ATF4_CEBPb/CEBPb.dat");
		
		r.printList(false);
		
		for (int i=159; i>156; i--) {
			r.fits.remove(i);
		}
		
		r.fits.remove(155);
		r.fits.remove(153);
		r.fits.remove(151);
		r.fits.remove(149);
		r.fits.remove(147);
		r.fits.remove(145);
		r.fits.remove(144);
		r.fits.remove(142);
		r.fits.remove(140);
		r.fits.remove(138);
		r.fits.remove(137);
		r.fits.remove(135);
		r.fits.remove(133);
		r.fits.remove(131);
		r.fits.remove(130);
		r.fits.remove(128);
		r.fits.remove(125);
		
		r.printList(false);
		
		r.printList("/Users/chaitanya/Desktop/CEBPb", true);
		
//		double tStart = System.currentTimeMillis();
//		for (long i=0; i<1E6; i++) {
//			swNucleotide(i);
//		}
//		System.out.println(System.currentTimeMillis()-tStart);
	}
	
	private static double swNucleotide(long input) {
		double totalSum = 0;
		double fSubSum	= 0;
		double rSubSum	= 0;
		long forwardSubString;
		long reverseSubString;
		long forwardStrand = fFlankingSequence | (input << 2*flankLength);
		long reverseStrand = rFlankingSequence | (reverseComplement(input, l) << 2*flankLength);
		
		for (int j=0; j<maxFrames; j++) {
			forwardSubString = forwardStrand & frameMask;
			reverseSubString = reverseStrand & frameMask;
			fSubSum = 0;
			rSubSum = 0;
			for (int loc=0; loc<k; loc++) {
				fSubSum += nucBetas[loc*4+((int) (forwardSubString&3))];
				rSubSum += nucBetas[loc*4+((int) (reverseSubString&3))];
				forwardSubString >>= 2;
				reverseSubString >>= 2;
			}
			totalSum += Math.exp(fSubSum) + Math.exp(rSubSum);
			forwardStrand >>= 2;
			reverseStrand >>= 2;
		}
		return totalSum;
	}
	
	protected static long reverseComplement(long input, int length) {
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
