package main;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Properties;
import java.util.zip.GZIPInputStream;

import base.DebugLog;
import base.Sequence;

public class FastQFileScanner
{
	/**
	 * Calculates PSFM for the entire FASTQ file
	 * @param path
	 * @param props
	 * @return
	 * @throws Exception
	 */
	public Object[]  calculatePSFM(String path, Properties props) throws Exception
	{
		return calculatePSFM( path, 0, null, null, props);
	}
	
	/**
	 * Calculates PSFM for the lines match left/right barcodes
	 * 
	 * @param path
	 * @param variableLength
	 * @param leftBarcode
	 * @param rightBarcode
	 * @param props
	 * @return
	 * @throws Exception
	 */
	public Object[]  calculatePSFM(String path, int variableLength,
			String leftBarcode, String rightBarcode,
			Properties props) throws Exception
	{
		boolean leftBarcodeAtIndex0 = true;		//default option, to be extended if needed	
		boolean checkBarcode = (leftBarcode!=null && rightBarcode!=null);
		DebugLog.log("calculatePSFM() : checkBarcode = "+ checkBarcode);
		if(checkBarcode)
		{
			DebugLog.log("leftBarcode    = " + leftBarcode);
			DebugLog.log("rightBarcode   = " + rightBarcode);
			DebugLog.log("variableLength = " + variableLength);
		}
		long[][] freqs = null;
		int k=0;
		
		BufferedReader br = null;
		if (path.endsWith(".gz"))
		{
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(
					new FileInputStream(path))));
		} else
		{
			br = new BufferedReader(new InputStreamReader(new FileInputStream(
					path)));
		}

		long lineCount = 0;
		long totalReads = 0;
		String line = null;
		long numberOfValidSeq = 0;

		DebugLog.log("Reading...");
		while ((line = br.readLine()) != null)
		{
			lineCount++;

			if (lineCount % 4 != 2)
				continue;
			boolean valid = true;

			totalReads++;

			int len = line.length();
			
			if(checkBarcode)
			{
				int idx1=-1;
				//Barcodes checking here
				if (leftBarcodeAtIndex0)
				{
					if (line.startsWith(leftBarcode))
						idx1 = 0;
				} else
				{
					idx1 = line.indexOf(leftBarcode);
				}
				
				int rightIdx = idx1 + leftBarcode.length() + variableLength;

				if (idx1 != -1 && line.regionMatches(rightIdx, rightBarcode, 0, rightBarcode.length()))
				{
					//pass
				}
				else
				{
					continue;
				}
			}
			
			//
			for (int i = 0; i < len; i++)
			{
				if (!Sequence.isValidSymbol(line.charAt(i)))
				{
					valid = false;
					break;
				}
			}
			//
			
			if(totalReads% (1000*1000*2) ==0)
			{
				DebugLog.log("Processed reads:"+totalReads);
			}
			if (valid)
			{
				numberOfValidSeq++;
				
				if( freqs == null)
				{
					k=line.length();
					freqs = new long[4][k];
				}
				
				for(int i=0;i<k;i++)
				{
					char c=line.charAt(i);
					int idx=(int)Sequence.getCharCode(c);
					freqs[idx][i]++;
				}
			}
			
		}

		DebugLog.log("Processed reads:"+totalReads);
		double[][] ratios = new double[4][k];
		for(int j=0;j<k;j++)
		{
			long total =0;
			for(int i=0;i<4;i++)
			{
				total += freqs[i][j];
			}

			for(int i=0;i<4;i++)
			{
				ratios[i][j] = freqs[i][j]*1.0/total;
				//System.out.print(ratios[i][j]+"\t");
			}
			//System.out.println();
		}

		Object[] results=new Object[4];
		for(int i=0;i<4;i++)
		{
			results[i] = ratios[i];
		}
		
		props.put("fastqTotalLine",lineCount);
		props.put("fastqTotalRead",totalReads);
		props.put("fastqValidRead",numberOfValidSeq);
		
		DebugLog.log(props);

		return results;
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception
	{
		Properties ps=new Properties();
		FastQFileScanner scanner=new FastQFileScanner();
//		Object[] psfm = scanner.calculatePSFM("./sample_data/exdUbx.exdScr.0.fastq.gz",ps);
//		printPFSM(psfm);
		Object[] psfm2 = scanner.calculatePSFM("./sample_data/exdUbx.exdScr.0.fastq.gz",
				16, "TGG", "CCAGCTG", ps);
		printPFSM(psfm2);
	}
	
	public static void printPFSM(Object[] psfm)
	{
		DebugLog.log("PSFM : ");
		int len=Array.getLength(psfm[0]);
		for(int i=0;i<len;i++)
		{
			StringBuffer sb=new StringBuffer();
			for(Object ds:psfm)
			{
				sb.append(Array.get(ds, i)+"\t");
			}
			DebugLog.log(sb.toString());
		}
	}

}
