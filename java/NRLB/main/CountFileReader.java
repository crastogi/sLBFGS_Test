package main;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Random;

import base.ArraysMergerHeap;
import base.CountObject;
import base.DebugLog;

/**
 * This is a tool to display the content of a count file.
 * It takes in the path on a count file, and outputs the 
 * list of the sequences and their counts.
 *
 */
public class CountFileReader
{

	/**
	 * @param args
	 * @throws FileNotFoundException 
	 */
	public static void main(String[] args) throws FileNotFoundException
	{
		//"tmp/kmax-4.dat"
		//String filePath = "tmp/traning.4.dat";
		if(args.length==0)
		{
			System.out.println("<count_file_path> [output_file_path]");
			return;
		}
		String filePath = args[0];
		PrintStream writer = null ; 
		if(args.length ==2)
		{
			writer = new PrintStream( args[1] );
		}
		else
		{
			writer = System.out;
		}
		ArraysMergerHeap.MinHeapNode node=new ArraysMergerHeap.MinHeapNode(filePath);
		CountObject obj=null;
		while((obj=node.peek())!=null)
		{
			writer.println(obj.toString2());
			node.pop();
		}
		writer.close();
	}

	private static Random randomObj = new Random();
	
	private static int getIdx(Long counts[], Long avgNum, int k)
	{
		int idx= randomObj.nextInt(k);
		while(counts[idx] > avgNum)
		{
			idx--;
			if(idx<0)
			{
				idx=k-1;
			}
		}
		return idx;
	}
	
	public static void split(String filePath, Long totalReads,Integer k,String folder) throws Exception
	{
		Long numPerFile = (long) Math.ceil(totalReads*1.0/k);
		DebugLog.log("Total Size:"+ totalReads);
		DebugLog.log("Avg Size:"+ numPerFile);
		ArraysMergerHeap.MinHeapNode node=new ArraysMergerHeap.MinHeapNode(filePath);
		
		DataOutputStream os[]=new DataOutputStream[k];
		Long counts[] =new Long[k];
		for(int i=0;i<k;i++)
		{
			counts[i] =0L;
			String path = folder+ "/"+(i)+".part" ;
			os[i] = new DataOutputStream(
					new BufferedOutputStream(new FileOutputStream(path), 10*1000*1000));
			DebugLog.log(path);
		}

		CountObject obj=null;
		while((obj=node.peek())!=null)
		{
			//System.out.println(obj.toString2());
			int idx= getIdx(counts , numPerFile, k);
			CountObject.steamOutSequence(os[idx],obj);
			counts[idx]++;
			//
			node.pop();
		}

		for(int i=0;i<k;i++)
		{
			CountObject.steamOutSequenceClose(os[i]);
			os[i].flush();
			DebugLog.log("Bucket["+(i)+"] Size:"+ counts[i]);
		}
		
	}
}
