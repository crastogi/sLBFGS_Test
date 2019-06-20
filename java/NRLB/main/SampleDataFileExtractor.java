package main;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Scanner;

import base.DebugLog;

public class SampleDataFileExtractor
{
	public static String baseFolder = "/data/selex/";
	public static String fileList = baseFolder + "filelist";
	
	public static ArrayList<String> dataFiles = new  ArrayList<String>();
		
	public static void extractSampleDataFiles(String outptuFolder) throws IOException
	{
		InputStream in = SampleDataFileExtractor.class.getResourceAsStream(fileList);
		Scanner s=new Scanner(in);
		while(s.hasNext())
		{
			String line = s.nextLine().trim();
			if(line.length()>0)
			{
				dataFiles.add(line);
			}
		}
		s.close();
		
		for(String name:dataFiles)
		{
			String path = baseFolder + name;
			String outputPath  = outptuFolder +"/" + name;
			DebugLog.log("Extracting the file : "+ path + "\t to "+outputPath);
			extractFile(path, outputPath);
		}
		DebugLog.log("Done");
	}
	
	public static void extractFile(String inPath, String outPath) throws IOException
	{
		InputStream in = SampleDataFileExtractor.class.getResourceAsStream(inPath);
		FileOutputStream os =new FileOutputStream(outPath);
		byte[] buffer=new byte[1024*1024*10];
		int len=-1;
		while((len=in.read(buffer))!=-1)
		{
			os.write(buffer, 0, len);
		}
		os.close();
		in.close();
	}
	
	public static void main(String[] args) throws IOException
	{
		extractSampleDataFiles(args[0]);
	}

}
