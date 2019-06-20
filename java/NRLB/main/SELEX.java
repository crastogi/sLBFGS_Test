package main;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.Properties;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;

import config.AdditionalDataSet;
import config.CrossValidationSetting;
import config.ExperimentReference;
import config.InformationGainStats;
import config.InputDataSetStats;
import config.KmerCountStats;
import config.MarkovModelStats;
import config.Model;
import config.Round0;
import config.SELEXSequencingConfig;
import config.SELEXStats;
import config.Sample;
import config.SequencingRunInfo;

import base.ArraysMergerHeap;
import base.CountObject;
import base.DebugLog;
import base.MarkovModelInfo;
import base.MarkovModelOption;
import base.RegexOption;
import base.SELEXConfigReader;
import base.Sequence;
import base.Util;


public class SELEX
{
	private static long start ;
	private static SELEXConfigReader reader;
	
	/** Measure time elapsed in SELEX.main()  
	 * 
	 */
	private static void logStart()
	{
		start =  System.currentTimeMillis();
	}

	/** Measure time elapsed in SELEX.main()  
	 * 
	 */
	private static void logEnd()
	{
		long end =  System.currentTimeMillis();
		long seconds = (end-start)/1000;
		DebugLog.log("Time Elapsed: " + seconds +" seconds.");
	}

	/** Returns the singleton of SELEXConfigReader
	 * 
	 */
	public static SELEXConfigReader getConfigReader()
	{
		if(reader==null)
		{
			reader =new SELEXConfigReader();
		}
		return reader;
	}
	
	/**
	 * Reads in the XML configuration.
	 * @param settingPath
	 * @param dataFolder If specified, it is used as the root folder of the relative paths in the XML config.
	 * @return a SELEXConfigReader object
	 */
	public static SELEXConfigReader loadProperties(String settingPath, String dataFolder)
	{
		DebugLog.log("Temp Folder:"+TEMP_FOLDER);
		
		try
		{
			SELEX.getConfigReader().readConfig(settingPath, dataFolder);
		} catch (Exception e)
		{
			DebugLog.log(e); throw new RuntimeException(e);
		}
		
		return reader;
	}

	/**
	 * Loads a SequencingRunInfo object
	 *  
	 * @param seqName
	 * @param seqFile
	 * @param sampleName
	 * @param round
	 * @param length
	 * @param left
	 * @param right
	 * @param leftFlank
	 * @param rightFlank
	 */
	public static void addSequenceInfo(String seqName,String seqFile,String sampleName,
			Integer round,Integer length,
			String left,String right, String leftFlank, String rightFlank)
	{
		addSequenceInfo( seqName, seqFile, sampleName,
				 round, length,
				 left, right,  leftFlank,  rightFlank,
				 null,  null);
	}
	
	/**
	 * Loads a SequencingRunInfo object
	 * 
	 * Called in R.
	 * 
	 * @param seqName
	 * @param seqFile
	 * @param sampleName
	 * @param round
	 * @param length
	 * @param left
	 * @param right
	 * @param leftFlank
	 * @param rightFlank
	 * @param round0SeqName
	 * @param round0Sample
	 */
	public static void addSequenceInfo(String seqName,String seqFile,String sampleName,
			Integer round,Integer length,
			String left,String right, String leftFlank, String rightFlank,
			String round0SeqName, String round0Sample)
	{
		if(reader==null)
		{
			reader =new SELEXConfigReader();
		}

		if(leftFlank==null)
		{
			leftFlank=left;
		}
		if(rightFlank==null)
		{
			rightFlank=right;
		}
		
		SequencingRunInfo info=new SequencingRunInfo();
		info.setName(seqName);
		info.setDataFile(seqFile);
		Sample sample=new Sample();
		sample.setName(sampleName);
		sample.setLeftBarcode(left);
		sample.setRightBarcode(right);
		sample.setLeftFlank(leftFlank);
		sample.setRightFlank(rightFlank);
		sample.setRound(round);
		sample.setVariableRegionLength(length);
		if(!Util.isEmpty(round0SeqName) && !Util.isEmpty(round0Sample))
		{
			Round0 r0 = new Round0();
			r0.setSequencingName(round0SeqName);
			r0.setSampleName(round0Sample);
			sample.setRound0(r0);
		}
		info.getSample().add(sample);
		reader.addSequencingConfig(info);
	}

	/** 
	 * Saves the current status to a XML file.
	 * 
	 * Called in R.
	 * 
	 */
	public static void saveStats()
	{	
		//xml
		if(SELEX.stats!=null)
		{
			try
			{
				String stats_path = getStatsXMLOutputPath();
				DebugLog.log("Saving stats to:"+stats_path);
				SELEXConfigReader.writeConfig(stats_path,SELEX.stats);
			}
			catch(Exception ex)
			{
				DebugLog.log(ex);
			}
		}
	}
	
	/**
	 * Loads a SELEXProcessConfig file and run the SELEX process.
	 * @param args
	 */
	public static void main(String[] args)
	{
		if(args.length!=2 && args.length!=3)
		{
			throw new RuntimeException("SELEX <a configuration file> <a working directory> [data folder]");
		}
		if(args.length == 3)
		{
			loadConfigFile(args[0], args[2]);
		}
		else
		{
			loadConfigFile(args[0], null);
		}
		setWorkingDirectory(args[1]);
		runSELEX();
		
	}
	
	private static String TEMP_FOLDER = "./SELEXworkspace/";

	/**
	 * 
	 * Loads the XML configuration when the workspace is set.
	 * 
	 * Called in R.
	 * 
	 * @param path
	 */
	public static void setWorkingDirectory(String path)
	{
		//
		if(!(path.endsWith("/")||!path.endsWith("\\")))
		{
			path = path + "/" ;
		}
		//
		createFolders(path);
		//
		String logFolder = path +"/log/";
		createFolders(logFolder);
		DebugLog.setDefaultLogFolder(logFolder);
		
		File f = new File(path);
		if(!f.isDirectory())
		{
			throw new RuntimeException("Path ["+path+"] should be a folder instead of a file.");
		}
		TEMP_FOLDER = f.getAbsolutePath()+"/";
		DebugLog.log("TempFolder:"+TEMP_FOLDER);
		loadStats();
	}
	
	private static void createFolders(String path)
	{
		File f = new File(path);
		if(!f.exists())
		{
			f.mkdirs();
		}
	}

	/**
	 * Get a formatted path string.
	 * 
	 * @param dataSet
	 * @param id
	 * @param offset
	 * @param asText
	 * @return
	 */
	private static String getOutputPath(ExperimentReference dataSet, String id,Integer offset, Boolean asText)
	{
		String dsID = getSampleID(dataSet);
		if(offset==null || offset<0)
			offset=-1;
		if(asText)
		{
			return TEMP_FOLDER+"/"+dsID+"."+id+".o"+offset+".txt";
		}
		else
		{
			return TEMP_FOLDER+"/"+dsID+"."+id+".o"+offset+".dat";
		}
	}
	
	/**
	 * Get a formatted path string.
	 * @param dataSet
	 * @param id
	 * @param asText
	 * @return
	 */
	public static String getOutputPath(ExperimentReference dataSet, String id, Boolean asText)
	{
		String dsID = getSampleID(dataSet);
		if(asText)
		{
			return TEMP_FOLDER+"/"+dsID+"."+id+".txt";
		}
		else
		{
			return TEMP_FOLDER+"/"+dsID+"."+id+".dat";
		}
	}
	
	public static String getOutputFolder(String dataSet)
	{
		return TEMP_FOLDER+"/"+dataSet+"/";
	}
	
	public static String getStatsOutputPath()
	{
		return TEMP_FOLDER+"/stats.txt";
	}

	public static String getStatsXMLOutputPath()
	{
		return TEMP_FOLDER+"/stats.xml";
	}
	
	/**
	 * Returns a list of newly created data set names and their counts.
	 * 
	 * K-mer counting of the variable length is performed to generate the binary cache file.
	 * The cache file is the binary representation of the sample in the original FASTQ file.
	 * The cache file is then split into several pieces.
	 * 
	 */
	public static Object[] splitDataSet( ExperimentReference ds, Integer variableLength,
			SimpleKmerCount obj, Double[] porportions, String name )
	{
		String dataSet  = getSampleID(ds);
		doMinimalCounting(ds, variableLength, obj, null, false,-1, false, -1, null);

		InputDataSetStats inputStats = getInputDataSetStats(ds, null);
		
		Long totalReads = inputStats.getValidReadCount();
		DebugLog.log("Total Reads of ["+dataSet+"]["+variableLength+"]: "+ totalReads);
		
		String outputFolder=getOutputFolder(dataSet+"_"+name);
		if(!new File(outputFolder).exists())
		{
			new File(outputFolder).mkdirs();
		}
		
		Properties ps=new Properties();
		List<InputDataSetStats> newDS=new ArrayList<InputDataSetStats>();
		Object[] results = obj.splitCacheFile(reader,  ds, ps, name ,totalReads, porportions, outputFolder,newDS);
		removeAdditionalDataSet(inputStats.getOriginalDataFilePath(), name);
		for(Entry entry:ps.entrySet())
		{
			AdditionalDataSet ads=new AdditionalDataSet();
			ads.setConfigurationPath((String)entry.getValue());
			ads.setExperimentReference(ds);
			ads.setOriginalDataFilePath(inputStats.getOriginalDataFilePath());
			ads.setTag(name);
			addAdditionalDataSet(ads);
		}
		
		for(InputDataSetStats stats:newDS)
		{
			SELEX.addInputDataSetStats(stats);
		}	
		
		saveStats();

		reader.reloadAdditionalSequenceConfig(SELEX.stats.getAdditionalDataSet());
		
		return results;
	}
	
	/**
	 * Queries a target InputDataSetStats object
	 * 
	 * Called in R
	 */
	public static InputDataSetStats getInputDataSetStats(ExperimentReference ref, String filter)
	{
		if(SELEX.stats==null)
			return null;
		if(filter==null)
			filter=new RegexOption().getVariableRegionRegexFormattedString();
		String id = getSampleID(ref);
		for(InputDataSetStats dataSet: SELEX.stats.getInputDataSetStats())
		{
			if(dataSet.getId().equals(id) &&  isStringEqual(filter, dataSet.getFilter()))
			{
				return dataSet;
			}
		}
		return null;
	}

	private static void addInputDataSetStats(InputDataSetStats seqName)
	{
		if(Util.isEmpty(seqName.getFilter()))
		{
			seqName.setFilter(new RegexOption().getVariableRegionRegexFormattedString());
		}
		List<InputDataSetStats> list=SELEX.stats.getInputDataSetStats();
		for(int i=0;i<list.size();i++)
		{
			InputDataSetStats stats = list.get(i);
			if(stats.getId().equals(seqName.getId()) && isStringEqual(seqName.getFilter(), stats.getFilter()))
			{
				DebugLog.log("Replacing <InputDataSetStats>="+seqName.getId());
				list.set(i,seqName);
				return;
			}
		}
		DebugLog.log("Adding <InputDataSetStats>="+seqName.getId());
		SELEX.stats.getInputDataSetStats().add(seqName);
	}

	private static void addKmerCountStats(KmerCountStats seqName)
	{
		List<KmerCountStats> list=SELEX.stats.getKmerCountStats();
		for(int i=0;i<list.size();i++)
		{
			KmerCountStats stats = list.get(i);
			if(stats.getId().equals(seqName.getId()) && stats.getLength() == seqName.getLength() 
					&& stats.getOffset() == seqName.getOffset() && isStringEqual(stats.getFilter(), seqName.getFilter()) )
			{
				DebugLog.log("Replacing <KmerCountStats>="+seqName.getId()+" Length:"+ stats.getLength()+ " Offset:"+stats.getOffset() );
				list.set(i,seqName);
				return;
			}
		}
		SELEX.stats.getKmerCountStats().add(seqName);
	}
	
	private static boolean isStringEqual(String a,String b)
	{
		if(Util.isEmpty(a) && Util.isEmpty(b)) return true;
		return a!=null && a.equals(b);
	}
	
	private static void addMarkovModelStats(MarkovModelStats seqName)
	{
		List<MarkovModelStats> list=SELEX.stats.getMarkovModelStats();
		for(int i=0;i<list.size();i++)
		{
			MarkovModelStats stats = list.get(i);
			if(stats.getId().equals(seqName.getId()) && stats.getModelLength() == seqName.getModelLength()
					&& isStringEqual(seqName.getFilter(), stats.getFilter())
					&& isStringEqual(seqName.getMarkovModelMethod(), stats.getMarkovModelMethod()))
			{
				DebugLog.log("Replacing <MarkovModelStats>="+seqName.getId());
				list.set(i,seqName);
				return;
			}
		}
		SELEX.stats.getMarkovModelStats().add(seqName);
	}


	private static void removeAdditionalDataSet(String orignalPath, String tag)
	{
		List<AdditionalDataSet> list=SELEX.stats.getAdditionalDataSet();
		for(int i=0;i<list.size();)
		{
			AdditionalDataSet stats = list.get(i);
			if(stats.getOriginalDataFilePath()!=null && 
					stats.getOriginalDataFilePath().equals(orignalPath) && 
					stats.getTag() !=null &&
					stats.getTag().equals(tag))
			{
				DebugLog.log("Removing <AdditionalDataSet>="+orignalPath+" TAG:"+tag);
				list.remove(i);
				continue;
			}
			i++;
		}
	}
	
	public static void addAdditionalDataSet(AdditionalDataSet seqName)
	{
		SELEX.stats.getAdditionalDataSet().add(seqName);
	}

	public static void addInformationGainStats(InformationGainStats seqName)
	{
		List<InformationGainStats> list=SELEX.stats.getInformationGainStats();
		for(int i=0;i<list.size();i++)
		{
			InformationGainStats stats = list.get(i);
			if(stats.getId().equals(seqName.getId()) && stats.getLength() == seqName.getLength()
					&& isStringEqual(stats.getFilter(), seqName.getFilter())
					&& isStringEqual(stats.getMarkovModelMethod(), seqName.getMarkovModelMethod()))
			{
				DebugLog.log("Replacing <InformationGainStats>="+seqName.getId());
				list.set(i,seqName);
				return;
			}
		}

		SELEX.stats.getInformationGainStats().add(seqName);
	}
	

	
	public static MarkovModelStats getMarkovStats(ExperimentReference ds, Integer length, RegexOption regexOption,
			String markovModelMethod)
	{
		String filter =Util.getKmerFilterString(regexOption);
		for(MarkovModelStats stats : SELEX.stats.getMarkovModelStats())
		{	
			if((length==stats.getModelLength())
					&&compareExperimentReferences(stats.getExperimentReference(),ds)
					&& isStringEqual(stats.getFilter(), filter)
					&& isStringEqual(stats.getMarkovModelMethod(), markovModelMethod))
			{
				return stats;
			}
		}
		return null;
	}
	
	/**
	 * Returns the K-mer count statistics
	 * 
	 * Called in R.
	 * 
	 * @param ds
	 * @param i
	 * @param offSet
	 * @param regexOption
	 * @return
	 */
	public static KmerCountStats getKMerCountStats(ExperimentReference ds, Integer i, Integer offSet,
			RegexOption regexOption)
	{
		String filter = Util.getKmerFilterString(regexOption);
		if(offSet==null || offSet <0)
			offSet=-1;
		for(KmerCountStats stats : SELEX.stats.getKmerCountStats())
		{	
			if((i==stats.getLength()) &&  offSet == stats.getOffset()
					&& compareExperimentReferences(stats.getExperimentReference(),ds))
			{
				if(isStringEqual(stats.getFilter(), filter))
					return stats;
			}
		}
		return null;
	}
	

	public static InformationGainStats getInformationGainStats(ExperimentReference ds, Integer i,
			String filter, String mmMethod)
	{
//		String filter = Util.getKmerFilterString(regexOption);
		for(InformationGainStats stats : SELEX.stats.getInformationGainStats())
		{	
			if((i==stats.getLength()) &&compareExperimentReferences(stats.getExperimentReference(),ds)
					&& isStringEqual(stats.getFilter(), filter)
					&& isStringEqual(stats.getMarkovModelMethod(), mmMethod) )
			{
				return stats;
			}
		}
		return null;
	}
	
	public static boolean compareExperimentReferences(ExperimentReference a,ExperimentReference b)
	{
		return a.getSequencingName().equals(b.getSequencingName()) &&
				 a.getSampleName().equals(b.getSampleName()) &&
				 a.getSampleRound() == b.getSampleRound();
	}

	public static boolean compareBarcodes(Sample a,Sample b)
	{
		return a.getLeftBarcode().equals(b.getLeftBarcode()) &&
		 		a.getRightBarcode().equals(b.getRightBarcode()) ;
	}

	/**
	 * Performs K-mer counting.
	 * 
	 * Calls doMinimalCounting() method, with 
	 * 		SimpleKmerCount is set to null.
	 * 		left offset is set to null.
	 * 
	 * Called in R.
	 */
	public static int doMinimalCounting(ExperimentReference ds, Integer i, 
			MarkovModelOption mmOption, Boolean saveOutputAsText, Integer threshold, RegexOption regexOption)
	{
		return doMinimalCounting( ds,  i,  null,  mmOption,  saveOutputAsText, null, false, threshold, regexOption);
	}

	/**
	 * Performs K-mer counting.
	 * 
	 * Calls doMinimalCounting() method, with 
	 * 		SimpleKmerCount is set to null.
	 * 		threshold is set to 100.
	 * 
	 * Called in R.
	 */
	public static int doMinimalCounting(ExperimentReference ds, Integer i,  
			MarkovModelOption mmOption, Boolean saveOutputAsText, Integer leftOffset, Boolean recalculate, RegexOption regexOption)
	{
		return doMinimalCounting( ds,  i,  null, mmOption,  saveOutputAsText,  leftOffset,  recalculate, 100, regexOption);
	}

	/**
	 * Does K-mer counting. Skips if the K-mer table is already calculated. 
	 * 
	 * @param ds Sample data set
	 * @param i	Length of K-mer
	 * @param obj The object that keeps k-mer counting results. If not provided, a local SimpleKmerCount object will be created.
	 * @param mmOption Indicates whether to build a Markov model and the Markov model method
	 * @param saveOutputAsText Saves the result in txt format. This mighe be useful for testing purpose.
	 * @param leftOffset At the specific offset (to the right of the left barcode) the K-mer counting will perform.
	 * 				If not provided or less then 0, the sliding window method is used. 
	 * @param recalculate Force recalculation.
	 * @param threshold Set the threshold for the minimal count. 
	 * 		The method returns 1 if the minimal count >= threshold, otherwise 0.
	 * @param regexOption Regular expression to filter lines and K-mers
	 * @return The method returns 1 if the minimal count >= threshold, otherwise 0.
	 */
	public static int doMinimalCounting(ExperimentReference ds, Integer i, SimpleKmerCount obj, 
			MarkovModelOption mmOption, Boolean saveOutputAsText, Integer leftOffset, Boolean recalculate, 
			Integer threshold, RegexOption regexOption)
	{
		if(mmOption==null)
		{
			mmOption = new MarkovModelOption();
		}
		
		String filter = Util.getKmerFilterString(regexOption);
		
		String dataSet = getSampleID(ds);
		String seqName = getSequenceID(ds);
		String outputPath=getOutputPath(ds, i+"", saveOutputAsText);
		int minimalCount  = threshold; 
		boolean overThreshold = false;
	
		KmerCountStats countStats = getKMerCountStats( ds,  i, leftOffset, regexOption);
		MarkovModelStats markovStats = getMarkovStats( ds, i, regexOption, mmOption.getMethod());

		if(obj==null)
		{
			obj=new SimpleKmerCount();
		}
		obj.setTempFolder(TEMP_FOLDER);
		obj.initTraining(SELEX.getConfigReader(), ds);
		obj.initKMerLength(i);
		obj.setRegexOption(regexOption);
		if(leftOffset!=null && leftOffset>=0) // disable the sliding window for a specific offset
		{
			obj.setLeftOffset(leftOffset);
			obj.setUseSlidingWindow(false);
			outputPath=getOutputPath(ds, i+"", leftOffset, saveOutputAsText);
		}
		outputPath = obj.setOutputPath(outputPath);
		obj.setSaveAsText(saveOutputAsText);
		
		boolean isBinary = obj.getIsBinaryFile();
		boolean isFASTQ = obj.getIsFASTQFile();
		
		//a quick check whether the K-mer count table is there. skip calculation if so
		if( (isBinary || !isFASTQ || new File(obj.getCacheFileName()).exists())
				&& countStats!=null 
				&& new File(countStats.getKmerCountStatsTablePath()).exists()
				&& !recalculate)
		{
			//
			DebugLog.log("KMerCount for ["+dataSet+"]["+i+"] has already been done. Skipping...");
			long minCount =countStats.getLowestCount();
			DebugLog.log("KMerCount.minCount for ["+dataSet+"]["+i+"] = "+ minCount);
			overThreshold =minCount> minimalCount;
			
			if(!mmOption.getBuildMarkovModel() || 
					(markovStats!=null && new File(markovStats.getTransitionTablePath()).exists()))
			{
				return overThreshold?1:0; 
			}
		}

		KmerCountStats subCountStats = getKMerCountStats( ds,  i-1, leftOffset, regexOption);

		if(mmOption.getBuildMarkovModel())//get sub k-mer for Markov model building
		{
			if(i!=1)
			{
				if(subCountStats==null || !new File(subCountStats.getKmerCountStatsTablePath()).exists())
				{
					DebugLog.log("SubKmerCountTable["+subCountStats.getKmerCountStatsTablePath()+"] not found.");
					DebugLog.log("Calculating for SubKmerCountTable["+subCountStats.getKmerCountStatsTablePath()+"]");
					doMinimalCounting( ds, i-1, obj, mmOption, false, -1, false, threshold, null);
				}
				subCountStats = getKMerCountStats( ds,  i-1, leftOffset, regexOption);
			}
		}
		
		DebugLog.log("Scanning ["+dataSet+"][ K = " + i +" ]");

		//perform k-mer counting
		try
		{
			obj.process();
		} catch (Exception e)
		{
			DebugLog.log(e); throw new RuntimeException(e);
		}

		//process k-mer counting results and stats
		Properties runStats= obj.getStats();
		Long lowestCount =(Long) runStats.get("lowestCount");
		Boolean fastqFileProcessed = (Boolean) runStats.get("fastqFileProcessed");
		DebugLog.log("runStats="+runStats);
		
		overThreshold= lowestCount > minimalCount;

		// save k-mer count stats
		KmerCountStats kmerCountStats =new KmerCountStats();
		kmerCountStats.setId(dataSet);
		kmerCountStats.setKmerCountStatsTablePath(outputPath);
		kmerCountStats.setLength(i);
		kmerCountStats.setLowestCount( (Long)runStats.get("lowestCount"));
		kmerCountStats.setHighestCount( (Long)runStats.get("highestCount"));
		if(leftOffset!=null && leftOffset>=0)
			kmerCountStats.setOffset(leftOffset);
		else
			kmerCountStats.setOffset(-1);
		kmerCountStats.setValidReadCount( (Long)runStats.get("noOfValidRead"));
		kmerCountStats.setExperimentReference(ds);
		kmerCountStats.setFilter(filter);
		addKmerCountStats(kmerCountStats);
		
		saveStats();
		DebugLog.log("fastqFileProcessed = "+ fastqFileProcessed);
		if(fastqFileProcessed ) 
				// no filters : processing the original fastq file , saving the stats
		{
			Long fastqTotalLineCount = (Long) runStats.get("fastqTotalLineCount");
			Long fastqTotalReadCount = (Long) runStats.get("fastqTotalReadCount");
			Long fastqBarcodeMatchedReadCount = (Long) runStats.get("fastqBarcodeMatchedReadCount");
			Long fastqValidReadCount = (Long) runStats.get("fastqValidReadCount");
			
			String fastqCacheFileName = (String) runStats.get("fastqCacheFileName");

			//xml stats
			InputDataSetStats dataSetStats= getInputDataSetStats( ds, Util.getVariableRegionFilterString(regexOption));

			if(dataSetStats==null)
			{
				dataSetStats=new InputDataSetStats();
				dataSetStats.setId(dataSet);
				dataSetStats.setOriginalDataFilePath(reader.getSequencingRunInfo(ds).getDataFile());
				dataSetStats.setCachedDataFilePath(fastqCacheFileName);
				dataSetStats.setExperimentReference(ds);
				Sample s=reader.getSample(ds);
				dataSetStats.setLeftBarcode(s.getLeftBarcode());
				dataSetStats.setRightBarcode(s.getRightBarcode());
				dataSetStats.setTotalLineCount(fastqTotalLineCount);
				dataSetStats.setTotalReadCount(fastqTotalReadCount);
				dataSetStats.setTotalBarcodeMatchedReadCount(fastqBarcodeMatchedReadCount);
				dataSetStats.setValidReadCount(fastqValidReadCount);
				dataSetStats.setFilter(Util.getVariableRegionFilterString(regexOption));
				DebugLog.log("Adding InputDataSetStats = "+ dataSetStats.getId() );
				addInputDataSetStats(dataSetStats);
			}
			else
			{
				DebugLog.log("Found InputDataSetStats = "+ dataSetStats.getId() );
			}
			
		}
		
		// build the Markov model if needed
		// this computes the transition probability table for given k-mer 
		// for instance, in the table P(ACGT) means P(T|ACG) 
		if(mmOption.getBuildMarkovModel())
		{
			if(mmOption.getMethod().contains(MarkovModelOption.TRANSITION) )
			{
				obj.buildMarkovModelEfficient(mmOption);
			}
			else
			{
				if(subCountStats==null && i==1)
				{
					obj.buildMarkovModelDivisionMethod(null,0,mmOption);
				}
				else
				{
					obj.buildMarkovModelDivisionMethod(
							subCountStats.getKmerCountStatsTablePath(),
							subCountStats.getValidReadCount(),mmOption);
				}
			}
			
			//xml
			MarkovModelStats modelStats =new MarkovModelStats();
			
			modelStats.setId(dataSet);
			modelStats.setExperimentReference(ds);
			modelStats.setModelLength(i); 
			modelStats.setModelOrder(i-1);
			modelStats.setFilter(filter);
			modelStats.setModelLogPath(mmOption.getModelLogFile());
			modelStats.setTransitionTablePath(mmOption.getModelOutputFile());
			modelStats.setMarkovModelMethod(mmOption.getMethod());
			addMarkovModelStats(modelStats);
		}

		saveStats();
		return overThreshold?1:0;
		
	}
	
	/**
	 * Build a Markov model without doing cross validation
	 * 
	 * Called in R.
	 * 
	 * @param ds
	 * @param modelLength
	 * @param obj
	 * @return
	 */
	public static MarkovModelInfo trainMarkovModel(ExperimentReference ds, Integer modelLength, RegexOption regexOption,
			String markovModelMethod)
	{
		return trainMarkovModel( ds, null, modelLength , 0, regexOption,markovModelMethod );
	}
	
	/**
	 * Builds a Markov model from a sample data set and performs cross validation if validationDS provided.
	 * The output of this method is a probability table for all k-mers, 
	 * 		the probability of the last nucleotide given first k-1 nucleotides  
	 * 
	 * Called in R.
	 * @param ds Sample data set.
	 * @param validationDS Validation data set.
	 * @param modelLength Markov order + 1
	 * @param kmax Kmax is the length of k-mers in the cross validation data set
	 * 				 for which the expected counts are calculated.
	 * @param regexOption
	 * @param markovModelMethod Possible options: TRANSITION, DIVISION, TRANSITION|WITH_LEFT_FLANK, DIVISION|WITH_LEFT_FLANK
	 * @return
	 */
	public static MarkovModelInfo trainMarkovModel(ExperimentReference ds, ExperimentReference validationDS, 
			Integer modelLength, Integer kmax, RegexOption regexOption, String markovModelMethod)
	{
		try
		{	
			if(markovModelMethod==null)
			{
				markovModelMethod  = MarkovModelOption.DIVISION;
			}

			String filter = Util.getKmerFilterString(regexOption);		
			String dataSet =  getSampleID(ds);

			KmerCountStats countStats = getKMerCountStats( ds,  modelLength, -1, regexOption);
			MarkovModelStats markovStats = getMarkovStats( ds,  modelLength, regexOption, markovModelMethod);
			
			//a rough check to see if the Markov model has been calculated
			if(countStats!=null && new File(countStats.getKmerCountStatsTablePath()).exists()
					&& markovStats!=null &&  new File( markovStats.getTransitionTablePath()).exists())
			{
				boolean calculationNeeded=false;
				//calculation needed if kmax doesn't match
				if(validationDS!=null)
				{
					CrossValidationSetting cv = markovStats.getCrossValidationSetting();
					if(cv==null || cv.getExperimentReference()==null ) //cross validation not done
					{
						calculationNeeded=true;
					}
					else
					{
						DebugLog.log("kmax:"+kmax+ " cv length:"+cv.getLength());
						
						if(!compareExperimentReferences(cv.getExperimentReference(), validationDS) 
								|| cv.getLength()!=kmax) // validation doesn't match
						{
							calculationNeeded=true;
						}
					}
					
				}
				
				if(!calculationNeeded)
				{
					DebugLog.log("MarkovModel for ["+dataSet+"]["+modelLength+
							"] has already been done. CalculationNeeded["+calculationNeeded+"]. Skipping...");
					MarkovModelInfo mm = getMarkModelInfo( ds ,validationDS, modelLength, regexOption, markovModelMethod);
					return mm;
				}
				
			}
			
			// calculates (k-1)-mer table for DIVISION method
			// this can be skipped for TRANSITION method
			if(modelLength!=1)
			{
				//to verify
				//
				KmerCountStats subCountStats = getKMerCountStats( ds,  modelLength-1, -1, regexOption);
				
				if(subCountStats==null || !new File(subCountStats.getKmerCountStatsTablePath()).exists() )
				{
					DebugLog.log("Sub Kmer ["+(modelLength-1)+"] is not calculated. Calculating ...");
					doMinimalCounting(ds, modelLength-1, null, null, false, -1, false, -1, regexOption );
				}
			}

			// pass an SimpleKmerCount obj in doMinimalCounting. 
			// it keeps the k-mer count table in memory, which is useful for cross validation
			SimpleKmerCount obj=new SimpleKmerCount();
			
			//
			MarkovModelOption mmOption = new MarkovModelOption();
			mmOption.setBuildMarkovModel(true);
			if(markovModelMethod!=null)
			{
				mmOption.setMethod(markovModelMethod);
			}
			
			//build markov model
			doMinimalCounting(ds, modelLength, obj, mmOption, false, -1, true, -1, regexOption);
			
			if(validationDS==null)//skip cross validation
			{
				MarkovModelInfo mm = getMarkModelInfo( ds ,validationDS, modelLength, regexOption, markovModelMethod);
				return mm;
			}

			//make sure kmax count table is there
			{
				SimpleKmerCount.includeReversed = false;
				doMinimalCounting(validationDS, kmax, null, null, false, -1, false, -1, null);
			}
			KmerCountStats validationCountStats = getKMerCountStats( validationDS,  kmax, -1, null);
			long noValidCounts=  validationCountStats.getValidReadCount();
			DebugLog.log("KMax["+kmax+"] Total Read:"+noValidCounts);
			
			try
			{
				String kmaxCountTablePath = validationCountStats.getKmerCountStatsTablePath();
				
				boolean useLeftFlank = false;
				if(markovModelMethod.contains(MarkovModelOption.WITH_LEFT_FLANK) )
				{
					useLeftFlank = true;
				}
				double R=obj.crossValidatePR(kmaxCountTablePath, noValidCounts, useLeftFlank);
				
				//xml config
				MarkovModelStats modelStats = getMarkovStats(ds,  modelLength, regexOption, markovModelMethod); 
				CrossValidationSetting cv=new CrossValidationSetting();
				cv.setExperimentReference(validationDS);
				cv.setLength(kmax);
				modelStats.setCrossValidationSetting(cv);
				modelStats.setPearsonCoefficient(R);
				modelStats.setFilter(filter);
				
				MarkovModelInfo mm = getMarkModelInfo( ds ,validationDS, modelLength, regexOption, markovModelMethod);
				
				saveStats();
				return mm;
			} catch (Exception e)
			{
				DebugLog.log(e);
				throw new RuntimeException(e);
			}
		}
		catch (Exception e)
		{
			DebugLog.log(e);
		}
		
		saveStats();
		return null;
	}
	
	/**
	 * Calculates information gain
	 */
	public static double calculateInformationGain(ExperimentReference ds, Integer i,
			  MarkovModelInfo mm, RegexOption regexOption)
	{
		return calculateInformationGain( ds,  i,  mm,  regexOption, false);
	}
	
	/**
	 * Calculates information gain
	 * 
	 * Called in R.
	 * 
	 * @param ds Sample to calculate
	 * @param i Length of the K-mer
	 * @param mm Markov model to predict expected counts
	 * @param regexOption 
	 * @param checkBarcode Verifies the barcodes of the sample and the barcodes of Markov model are the same if specified true.
	 * @return
	 */
	public static double calculateInformationGain(ExperimentReference ds, Integer i,
			  MarkovModelInfo mm, RegexOption regexOption, Boolean checkBarcode)
	{	
		String filter = "InfoGainFilter:"+ Util.getKmerFilterString(regexOption)
						+ "\nMarkovModelFilter:"+mm.getFilter();
		
		String dataSet = getSampleID(ds);
		
		if(checkBarcode)
		{
			Sample b = reader.getSample(ds) ; 
			Sample a = reader.getSample(mm.getSample());
			if(!compareBarcodes(a,b))
			{
				throw new RuntimeException("Barcodes in the infomration gain sample don't match those of Markov model's.");
			}
		}

		KmerCountStats countStats = getKMerCountStats( ds,  i,  -1, regexOption );
		InformationGainStats infoGainStats = getInformationGainStats( ds,  i, filter, mm.getMarkovModelMethod());
		if(infoGainStats!=null && new File(countStats.getKmerCountStatsTablePath()).exists() )
		{
			DebugLog.log("InfomationGain for ["+dataSet+"]["+i+"] has already been done. Skipping...");
			double info = infoGainStats.getInformationGainValue();
			return info;
		}
		
		SimpleKmerCount obj=new SimpleKmerCount();
		doMinimalCounting(ds, i, obj, null, false, -1, true, -1, regexOption);
		
		try
		{
			double info= obj.calculateInformationGain(mm);
			
			//xml config
			InformationGainStats infoStats =new InformationGainStats();
			infoStats.setExperimentReference(ds);
			infoStats.setId(dataSet);
			infoStats.setInformationGainValue(info);
			infoStats.setLength(i);
			infoStats.setFilter(filter);
			infoStats.setModelID(getSampleID(mm.getSample())+".Order"+(mm.getMarkovLength()-1));
			infoStats.setMarkovModelMethod(mm.getMarkovModelMethod());
			addInformationGainStats(infoStats);
			
			saveStats();
			return info;
		}
		catch(Exception ex)
		{
			DebugLog.log(ex);
		}
		return -1;
	}
	
	/**
	 * Returns a formatted string ID of the sample
	 * 
	 * Called in R.
	 * 
	 * @param dataSet
	 * @return
	 */
	public static String getSampleID(ExperimentReference dataSet)
	{
		String dsID = dataSet.getSequencingName() +"."+ dataSet.getSampleName()+"."+dataSet.getSampleRound();
		return dsID;
	}
	
	public static String getSequenceID(ExperimentReference dataSet)
	{
		String dsID = dataSet.getSequencingName();
		return dsID;
	}
	
	public static MarkovModelInfo getMarkModelInfo(ExperimentReference dataSet, ExperimentReference cvDataSet, 
			Integer markovLength, RegexOption regexOption, String mmMethod)
	{
		//old way
		KmerCountStats countStats = getKMerCountStats( dataSet,  markovLength, -1, regexOption);
		//alternative way 
		//KmerCountStats countStats = getKMerCountStats( dataSet,  markovLength, 0, regexOption);
		MarkovModelStats markovStats = getMarkovStats( dataSet, markovLength, regexOption, mmMethod);
		String markovCountsPath =countStats.getKmerCountStatsTablePath();
		String markovObjPath = markovStats.getTransitionTablePath();
		long markovLengthTotalCount = countStats.getValidReadCount();
		Double markovR2 =markovStats.getPearsonCoefficient();
		
		MarkovModelInfo mm = new MarkovModelInfo();
		mm.setFilter(markovStats.getFilter());
		mm.setMarkovLength( markovLength);
		mm.setMarkovLengthTotalCount( markovLengthTotalCount);
		mm.setMarkovCountsPath( markovCountsPath);
		mm.setMarkovObjPath( markovObjPath );
		mm.setMarkovR2( markovR2);
		mm.setSample(dataSet);
		mm.setCrossValidationSample(cvDataSet);
		mm.setMarkovModelMethod(mmMethod);
		
		return mm;
	}
	
	/**
	 * Calculates the k-max of a sample
	 * @param testingDataSet
	 * @return
	 */
	public static int kmax(ExperimentReference testingDataSet)
	{
		//calculate kmax
		int kmax=0;
		int threshold = 100;
		
		MarkovModelOption options=new MarkovModelOption();
		options.setBuildMarkovModel(true);
		for(int i=1;i<16;i++)
		{
			kmax = i-1;
			
			//returning 0/1 instead of boolean to 
			int result = doMinimalCounting(testingDataSet, i, null, options, false, -1, false, threshold, null);
			
			if(result==0)
				break;
			
			saveStats();
			
		}

		InputDataSetStats stats = getInputDataSetStats(testingDataSet, null); 
		stats.setKmax(kmax);
		return kmax;
	}
	
	/**
	 * Runs the standard SELEX process.
	 * Calculates K-max of the testing samples.
	 * Builds the optimal Markov model from the training samples.
	 * Calculates the highest information gain for the later round samples.
	 */
	public static void runSELEX()
	{
		//checkConfigFilePath();

		String mmMethod = MarkovModelOption.DIVISION;
		Model model =null;
		try{
			model = reader.getProcessConfig().getProcess().getModel(); 
		}
		catch(Exception ex)
		{
			throw new RuntimeException("Please provide a SELEX process configuration file.",ex);
		}

		ExperimentReference trainingDataSet = model.getDataSource().getExperimentReference();
		ExperimentReference testingDataSet  = model.getCrossValidation().getExperimentReference();
		ExperimentReference infoGainDataSet = model.getInformationGain().getExperimentReference();
		
		logStart();
		try
		{
			//calculate kmax
			int kmax=0;
			int markovLength=0;
			int infoGainRound=0;
			int threshold = 100;

			InputDataSetStats stats = getInputDataSetStats(testingDataSet, null); 
			if(stats!=null)
			{
				Integer km= stats.getKmax();
				if(km!=null)
				{
					kmax = km;
				}
			}
			
			//Calculate K-max
			if(kmax==0)
			{	
				for(int i=1;i<16;i++)
				{
					kmax = i-1;
					
					//returning 0/1 instead of boolean to 
					int result = doMinimalCounting(testingDataSet, i, null, null, false, -1, false, threshold, null);
					
					if(result==0)
						break;
					
				}
			}

			stats = getInputDataSetStats(testingDataSet,null); 
			stats.setKmax(kmax);
			saveStats();
			
			DebugLog.log("Kmax = " + kmax);
			
			if(markovLength==0)
			{	
				//training Markov model
				
				double maxR=0;
				for(int i=1;i<=kmax;i++)
				{
					double R=trainMarkovModel(trainingDataSet, testingDataSet, i ,  kmax ,null, mmMethod).getMarkovR2();
					if(maxR<R)
					{
						maxR=R;
						markovLength=i;
					}

					saveStats();
				}
				
			}
			
			DebugLog.log("markov.optimal.length = " + markovLength);
			//System.exit(0);

			if(infoGainRound==0)
			{
				infoGainRound = markovLength ; 
			}
			else
			{
				infoGainRound++;
			}

			DebugLog.log("information gain starting from k = " + infoGainRound);
			
			//calculate information gain

			MarkovModelInfo mm = getMarkModelInfo( trainingDataSet,testingDataSet,  markovLength, null, mmMethod);
			
			for(int i=kmax;i<=16;i++)
			{
				calculateInformationGain(infoGainDataSet, i ,  mm, null);
				saveStats();
			}

			System.out.println("===========INFOGAIN========");
			Object[] obj1=getInformationGain();
			for(int i=0;i<Array.getLength(obj1);i++)
			{
				Object obj11 = Array.get(obj1, i);

				for(int j=0;j<Array.getLength(obj11);j++)
				{
					System.out.println(Array.get(obj11,j));
				}
			}

			System.out.println("===========MARKOV========");
			Object[] obj2=getMarkovPR();
			for(int i=0;i<Array.getLength(obj2);i++)
			{
				Object obj11 = Array.get(obj2, i);

				for(int j=0;j<Array.getLength(obj11);j++)
				{
					System.out.println(Array.get(obj11,j));
				}
			}
			
		} catch (Exception e)
		{ 
			// TODO Auto-generated catch block
			DebugLog.log(e);
			
		}
		logEnd(); 
	
	}
	
	private static SELEXStats stats = null;
	
	/**
	 * Loads the XML configuration when the workspace is set.
	 */
	private static void loadStats()
	{
		//xml stats
		String stats_path = getStatsXMLOutputPath();
		if(new File(stats_path).exists())
		{
			try
			{
				JAXBContext context = JAXBContext.newInstance("config");
		        Unmarshaller marshaller = context.createUnmarshaller();
		        stats = (SELEXStats) marshaller.unmarshal(new FileInputStream(stats_path));
			} 
			catch (Exception e)
			{
				DebugLog.log(e);
				//throw new RuntimeException("Can't load config file:"+stats_path);
			}
			try
			{
				SELEX.getConfigReader().readAdditionalDS(stats.getAdditionalDataSet());
			} 
			catch (Exception e)
			{
				DebugLog.log(e);
				//throw new RuntimeException("Can't load config file:"+stats_path);
			}
		}
		if(stats==null)
		{
			stats=new SELEXStats();
		}
		
	}
	
	/** Returns all the calculated information gain values
	 * 
	 * Called in R.
	 * 
	 */
	public static Object[] getInformationGain()
	{
		//checkConfigFilePath();
	    Pattern r = Pattern.compile("(.+)\\.(\\d+)\\.gain");
	    
	    ArrayList<Double> informationGains=new ArrayList<Double>();
	    ArrayList<Integer> informationGainRounds=new ArrayList<Integer>();
	    ArrayList<String> samplesList=new ArrayList<String>();
	    ArrayList<String> markovList=new ArrayList<String>();
	    ArrayList<String> filterList=new ArrayList<String>();
	    ArrayList<String> methodList=new ArrayList<String>();

		for(InformationGainStats infoStats : SELEX.stats.getInformationGainStats())
		{
			String key = getSampleID(infoStats.getExperimentReference());
			
			informationGains.add(infoStats.getInformationGainValue());
			informationGainRounds.add(infoStats.getLength());
			samplesList.add(key);
			markovList.add(infoStats.getModelID());
			filterList.add(infoStats.getFilter());
			methodList.add(infoStats.getMarkovModelMethod());
		}
		
		double[] values=new double[informationGains.size()];
		for(int i=0;i<values.length;i++)values[i]=informationGains.get(i);
		int[] rounds=new int[informationGainRounds.size()];
		for(int i=0;i<values.length;i++)rounds[i]=informationGainRounds.get(i);
		String[] samples=new String[informationGainRounds.size()];
		for(int i=0;i<values.length;i++)samples[i]=samplesList.get(i);

		String[] markovListArray=new String[informationGainRounds.size()];
		for(int i=0;i<values.length;i++)markovListArray[i]=markovList.get(i);

		String[] filterListArray=new String[informationGainRounds.size()];
		filterList.toArray(filterListArray);
		
		String[] methodListArray=new String[informationGainRounds.size()];
		methodList.toArray(methodListArray);
		
		
		
		return new Object[]{samples, rounds, values, markovListArray, filterListArray, methodListArray};
	}
	
	/**
	 * Returns all the calculated Markov model statistics.
	 * 
	 * Called in R.
	 * 
	 * @return
	 */
	public static Object[] getMarkovPR()
	{
		//checkConfigFilePath();
	    
	    ArrayList<Double> prs=new ArrayList<Double>();
	    ArrayList<Integer> roundNums=new ArrayList<Integer>();
	    ArrayList<String> samplesList=new ArrayList<String>();
	    ArrayList<String> cvList=new ArrayList<String>();
	    ArrayList<Integer> cvLength=new ArrayList<Integer>();
	    ArrayList<String> filterList=new ArrayList<String>();
	    ArrayList<String> methodsList=new ArrayList<String>();
	    
		for(MarkovModelStats markovStats:SELEX.stats.getMarkovModelStats())
		{
			String key = getSampleID(markovStats.getExperimentReference());
			
			//DebugLog.log(m.group(1)+ " = "+ entry.getValue());
			if(markovStats.getPearsonCoefficient()==null)
				continue;
			prs.add(markovStats.getPearsonCoefficient());
			roundNums.add(markovStats.getModelLength());
			samplesList.add(key);
			CrossValidationSetting cv=markovStats.getCrossValidationSetting();
			if(cv!=null)
			{
				cvList.add(getSampleID(cv.getExperimentReference()));
				cvLength.add(cv.getLength());
			}
			else
			{
				cvList.add("NA");
				cvLength.add(0);
			}
			filterList.add(markovStats.getFilter());
			methodsList.add(markovStats.getMarkovModelMethod());
		}
		
		double[] values=new double[roundNums.size()];
		for(int i=0;i<values.length;i++)
			values[i]=prs.get(i);
		int[] rounds=new int[roundNums.size()];
		for(int i=0;i<values.length;i++)
			rounds[i]=roundNums.get(i);
		String[] samples=new String[roundNums.size()];
		for(int i=0;i<values.length;i++)
			samples[i]=samplesList.get(i);
		
		String[] cvListArray=new String[roundNums.size()];
		for(int i=0;i<values.length;i++)
			cvListArray[i]=cvList.get(i);
		int[] cvLengthArray=new int[roundNums.size()];
		for(int i=0;i<values.length;i++)
			cvLengthArray[i]=cvLength.get(i);

		String[] filterArray=new String[roundNums.size()];
		filterList.toArray(filterArray);
		
		String[] methodArray=new String[roundNums.size()];
		methodsList.toArray(methodArray);
		
		return new Object[]{samples, rounds, values, cvListArray, cvLengthArray,filterArray,methodArray};
	}
	
	/**
	 * Returns all the K-mer counting statistics
	 * 
	 * Called in R.
	 * 
	 * @return
	 */
	public static Object[] getCountStats()
	{
		ArrayList<String> sampleList=new ArrayList<String>();
		ArrayList<Integer> kList=new ArrayList<Integer>();
		ArrayList<Integer> offSetList=new ArrayList<Integer>();
		ArrayList<Long> minCountList=new ArrayList<Long>();
		ArrayList<Long> maxCountList=new ArrayList<Long>();
		ArrayList<Long> numOfReadsList=new ArrayList<Long>();
	    ArrayList<String> filterList=new ArrayList<String>();
	    
		for(KmerCountStats kmerStats : SELEX.stats.getKmerCountStats())
		{
			sampleList.add(getSampleID(kmerStats.getExperimentReference()));
			kList.add(kmerStats.getLength());
			offSetList.add(kmerStats.getOffset());
			minCountList.add(kmerStats.getLowestCount());
			maxCountList.add(kmerStats.getHighestCount());
			numOfReadsList.add(kmerStats.getValidReadCount());
			filterList.add(kmerStats.getFilter());
		}
		
		String[] sampleListArray=new String[sampleList.size()];
		sampleList.toArray(sampleListArray);
		int[] kListArray=new int[sampleList.size()];
		for(int i=0;i<sampleList.size();i++)kListArray[i]=kList.get(i);
		int[] offSetListArray=new int[sampleList.size()];
		for(int i=0;i<sampleList.size();i++)offSetListArray[i]=offSetList.get(i);
		long[] minCountListArray=new long[sampleList.size()];
		for(int i=0;i<sampleList.size();i++)minCountListArray[i]=minCountList.get(i);
		long[] maxCountListArray=new long[sampleList.size()];
		for(int i=0;i<sampleList.size();i++)maxCountListArray[i]=maxCountList.get(i);
		long[] numOfReadsListArray=new long[sampleList.size()];
		for(int i=0;i<sampleList.size();i++)numOfReadsListArray[i]=numOfReadsList.get(i);

		String[] filterArray=new String[sampleList.size()];
		filterList.toArray(filterArray);
		
		return new Object[]{sampleListArray, kListArray, offSetListArray, 
				minCountListArray, maxCountListArray,numOfReadsListArray, filterArray};
	}
	
	/**
	 * Splits the sample into multiple pieces
	 * 
	 * Called in R.
	 * 
	 * @param dataSet
	 * @param newDataSetPrerix
	 * @param ratios
	 * @return
	 */
	public static Object[] splitDataSet( ExperimentReference dataSet, String newDataSetPrerix, ArrayList<Double> ratios)
	{
		Double[] proportions=null;
		if(ratios==null || ratios.size() ==0)
		{
			proportions = new Double[]{0.5, 0.5};
		}
		else
		{
			proportions = new Double[ratios.size()];
			ratios.toArray(proportions);
		}
		SimpleKmerCount obj=new SimpleKmerCount();
		obj.initTraining(SELEX.getConfigReader(), dataSet);
		obj.setTempFolder(TEMP_FOLDER);
		
		Sample sample = reader.getSample(dataSet );
		Integer variableLength = sample.getVariableRegionLength();
		Object[] results = splitDataSet( dataSet, variableLength, 
				 obj,proportions , newDataSetPrerix);
		//DebugLog.log(Arrays.toString(newDS));
		return results;
	}
	

	public static void loadConfigFile(String configPath)
	{
		if(!new File(configPath).exists())
		{
			throw new RuntimeException("File not found:"+configPath);
		}
		loadProperties(configPath, null);
	}
	
	/**
	 * Reads in the XML configuration.
	 * 
	 * Call in R.
	 * 
	 * @param settingPath
	 * @param dataFolder If specified, it is used as the root folder of the relative paths in the XML config.
	 *
	 */
	public static void loadConfigFile(String configPath, String dataFolder)
	{
		if(!new File(configPath).exists())
		{
			throw new RuntimeException("File not found:"+configPath);
		}
		loadProperties(configPath, dataFolder);
	}

	/**
	 * Output the currently loaded samples to a xml file
	 * 
	 * Called in R.
	 * 
	 * @param configPath
	 */
	public static void exportConfigFile(String configPath)
	{
		SELEX.getConfigReader().exportSamples(configPath);
	}

	/**
	 * Returns the Round0 object of a sample if it exists.
	 * 
	 * Called in R.
	 * 
	 * @param s
	 * @return
	 */
	public static ExperimentReference getRound0(ExperimentReference s)
	{
		return SELEX.getConfigReader().getRound0(s);
	}
	
	/**
	 * Returns the attributes of the currently loaded samples
	 * 
	 * Called in R.
	 * 
	 * @return
	 */
	public static Object[] showSamples()
	{
		if(reader!=null)
		{
			return reader.showSamples();
		}
		else
		{
			return new Object[0];
		}
	}
	
	/**
	 * Clear the sample references in the memory.
	 * 
	 * Called in R.
	 */
	public static void unloadSamples()
	{
		if(reader!=null)
		{
			reader.clear();
		}
	}
	
	/**
	 * Calculates PSFM table for K-mers
	 * 
	 * Called in R.
	 * 
	 * @param ds
	 * @param k
	 * @param offset
	 * @return
	 */
	public static Object[] getPSFM(ExperimentReference ds, Integer k, Integer offset)
	{
		doMinimalCounting( ds, k,  null, null, false, offset, false, -1, null); 

		KmerCountStats countStats = getKMerCountStats( ds, k, offset, null);
		String filePath = countStats.getKmerCountStatsTablePath();
		ArraysMergerHeap.MinHeapNode node=new ArraysMergerHeap.MinHeapNode(filePath);
		
		CountObject obj=null;
		long[][] counts = new long[4][k];
		while((obj=node.peek())!=null)
		{
			Sequence seq = (Sequence) obj.getKey();
			//
			String str = seq.getString();
			for(int i=0;i<k;i++)
			{
				char c=str.charAt(i);
				int idx=(int)Sequence.getCharCode(c);
				counts[idx][i]+=obj.getCount();
			}
			//
			node.pop();
		}
		
		double[][] ratios = new double[4][k];
		for(int j=0;j<k;j++)
		{
			long total =0;
			for(int i=0;i<4;i++)
			{
				total += counts[i][j];
			}

			for(int i=0;i<4;i++)
			{
				ratios[i][j] = counts[i][j]*1.0/total;
			}
		}
		
		Object[] results=new Object[4];
		for(int i=0;i<4;i++)
		{
			results[i] = ratios[i];
		}
		return results;
	}
	
	/**
	 * Calculates PSFM for the whole FASTQ file.
	 * 
	 * Called in R.
	 * 
	 * @param ref
	 * @return
	 */
	public static Object[] getSamplePSFM(ExperimentReference ref)
	{
		SequencingRunInfo runInfo = reader.getSequencingRunInfo(ref);
		String path = runInfo.getDataFile();
		
		Sample sample = reader.getSample(ref);
		
		Properties ps=new Properties();
		FastQFileScanner scanner=new FastQFileScanner();
		Object[] results;
		try
		{
			results = scanner.calculatePSFM(path,sample.getVariableRegionLength(), 
					sample.getLeftBarcode(), sample.getRightBarcode(), ps);
		} 
		catch (Exception e)
		{
			DebugLog.log(e);
			throw new RuntimeException(e);
		}
		
		return results;
	}
	
	/**
	 * Calculates PSFM for the whole FASTQ file.
	 * 
	 * Called in R.
	 * 
	 * @param seqName
	 * @return
	 */
	public static Object[] getFastqPSFM(String seqName)
	{
		String path = reader.getSequencingRunInfoDataFile(seqName);
		DebugLog.log("Fastq file:"+path);
		if(path==null)
		{
			throw new RuntimeException("Fastq file with seqName["+seqName+"] not found.");
		}
		
		Properties ps=new Properties();
		FastQFileScanner scanner=new FastQFileScanner();
		Object[] results;
		try
		{
			results = scanner.calculatePSFM(path,ps);
		} catch (Exception e)
		{
			DebugLog.log(e);
			throw new RuntimeException(e);
		}
		
		return results;
	}
	
	/**
	 * Calculates PSFM for the whole FASTQ file by sample
	 * 
	 * Called in R.
	 * 
	 * @param seqName
	 * @return
	 */
	public static Object[] getFastqPSFM(ExperimentReference ref)
	{
		SequencingRunInfo runInfo = reader.getSequencingRunInfo(ref);
		String path = runInfo.getDataFile();
		
		DebugLog.log("Fastq file:"+path);
		
		Properties ps=new Properties();
		FastQFileScanner scanner=new FastQFileScanner();
		Object[] results;
		try
		{
			results = scanner.calculatePSFM(path,ps);
		} catch (Exception e)
		{
			DebugLog.log(e);
			throw new RuntimeException(e);
		}
		
		return results;
	}
	
	/**
	 * Retrieves the K-mer counting results
	 */
	public static Object[] getKmerCount(ExperimentReference ds, Integer k, Integer offset,
			Integer minCount, Integer top, MarkovModelInfo mm, RegexOption regexOption )
	{
		return  getKmerCount( ds,  k,  offset,
				 minCount,  top,  mm,  regexOption, "");
	}

	/**
	 * Retrieves the K-mer counting results
	 * 
	 * Called from R
	 */
	public static Object[] getKmerCount(ExperimentReference ds, Integer k, Integer offset,
			Integer minCount, Integer top, MarkovModelInfo mm, RegexOption regexOption , String outputPath)
	{
		try
		{
		
		//checkConfigFilePath();
		//String filter = regexOption;
		//String readRegex = ;

		KmerCountStats kcountStats = getKMerCountStats( ds,  k, offset, regexOption);
		
		Sample sample = reader.getSample(ds);
		String leftFlank = sample.getLeftFlank();
		
		String filePath = kcountStats.getKmerCountStatsTablePath();
		//ArraysMergerHeap.MinHeapNode node=new ArraysMergerHeap.MinHeapNode(filePath, 1024*1024*100, 10*1000*1000);
		ArraysMergerHeap.MinHeapNode node=new ArraysMergerHeap.MinHeapNode(filePath );

		class LocalCountObjectComparator implements Comparator<ArrayList>
		{
		    public int compare( ArrayList x, ArrayList y )
		    {
		    	CountObject x0=(CountObject)x.get(0);
		    	CountObject y0=(CountObject)y.get(0);
		        return y0.compareTo(x0);
		    }
		}
		
		Pattern selectPattern =null;
		Pattern excludePattern = null;
		Pattern selectKmerPattern = null;
		if(regexOption!=null){
			selectPattern = regexOption.getViewIncludePattern();
			excludePattern = regexOption.getViewExcludePattern();
			selectKmerPattern = regexOption.getViewIncludeOnlyPattern();
			//minCount = -1;
		}
		DebugLog.log("minCount is set to "+minCount);
		boolean useTop=true;
		if(top==null || top <=0)
		{
			useTop=false;
		}
		DebugLog.log("useTop is set to "+useTop);
		PriorityQueue<ArrayList> list=new PriorityQueue<ArrayList>(100,new LocalCountObjectComparator());
		
		float[] marKovProbabilities = null;
		String morkovObjPath = null;
		String markovCountPath = null;
		int markovModelLength = 0;
		int[] markovCounts = null;
		long markovTotalCount = 0;
		long kNoOfValidRead = 0;
		boolean useMarkov=false;
		
		boolean useLeftFlank=false;
		
		if(mm!=null && mm.getMarkovObjPath()!=null)
		{
			useLeftFlank =  mm.getMarkovModelMethod().toUpperCase().contains(MarkovModelOption.WITH_LEFT_FLANK);
			useMarkov=true;
			KmerCountStats countStats = getKMerCountStats(ds, k, -1, regexOption);
			kNoOfValidRead = countStats.getValidReadCount();
			try
			{
				morkovObjPath = mm.getMarkovObjPath();
				markovCountPath = mm.getMarkovCountsPath();
				markovModelLength = mm.getMarkovLength();
				markovTotalCount = mm.getMarkovLengthTotalCount();
				
				// load markovProbabilities
				DebugLog.log("Reading Markov Prob file:" + morkovObjPath);
				FileInputStream fos = new FileInputStream(morkovObjPath);
				ObjectInputStream oos = new ObjectInputStream(fos);
				marKovProbabilities = (float[]) oos.readObject();

				DebugLog.log("Reading Markov Count file:" + markovCountPath);
				ArraysMergerHeap.MinHeapNode node2 = new ArraysMergerHeap.MinHeapNode(
						markovCountPath);
				CountObject obj = null;
				markovCounts = new int[1 << (2 * (markovModelLength))];

				while ((obj = node2.peek()) != null)
				{
					Sequence seq = (Sequence) obj.getKey();
					markovCounts[(int) seq.getValue()] = obj.getCount();
					node2.pop();
				}
				
			}
			catch(Exception ex)
			{
				DebugLog.log(ex);
				throw new RuntimeException(ex);
			}
		}

		String[] header1 = new String[]{"Kmer","ObservedCount"} ; 
 		String[] header2 = new String[]{"Kmer","ObservedCount","Probability","ExpectedCount"};
 		
		PrintWriter writer = null;
		if(!Util.isEmpty(outputPath))
		{
			writer =new PrintWriter(new BufferedOutputStream(new FileOutputStream(outputPath),1024*1024*10));
			DebugLog.log("Output counts to file: "+ outputPath);
		}
		
		Boolean directOutputWithoutBuildingArray = (writer !=null && !useTop) ;
		if(directOutputWithoutBuildingArray)
		{
			if(!useMarkov)
			{
				Util.output1DArray(writer, header1);
			}
			else
			{
				Util.output1DArray(writer, header2);	
			}
		}
		
		CountObject obj=null;
		int lineNumber = 0;
		int STEP = 10000000;
		int NEXT_STEP = 10000000;
		while((obj=node.peek())!=null)
		{
			lineNumber++;
			
			if(lineNumber==NEXT_STEP)
			{
				DebugLog.log("Lines processed:"+ lineNumber);
				NEXT_STEP += STEP;
			}
			
			if(obj.getCount() >= minCount)
			{
				if(selectPattern!=null)
				{
					Matcher m=selectPattern.matcher(obj.getKey().getString());
					if(!m.find())
					{
						//DebugLog.log("Skipping..."+obj.getKey().getString());
						node.pop();
						continue;
					}
				}
				if(excludePattern!=null)
				{
					Matcher m=excludePattern.matcher(obj.getKey().getString());
					if(m.find())
					{
						//DebugLog.log("Skipping..."+obj.getKey().getString());
						node.pop();
						continue;
					}
				}
				
				if(selectKmerPattern!=null)
				{
					Matcher m=selectKmerPattern.matcher(obj.getKey().getString());
					if(!m.find())
					{
						//DebugLog.log("Skipping..."+obj.getKey().getString());
						node.pop();
						continue;
					}
				}
//				if(selectKmers!=null)
//				{
//					if(!selectKmers.contains(obj.getKey().getString()))
//						continue;
//				}
				
				ArrayList results=new ArrayList();
				results.add(obj);
				results.add(obj.getKey());
				results.add(obj.getCount());
				
				//
				if(useMarkov)
				{
					Sequence seq = (Sequence) obj.getKey();
					double modelProb = 1;
					
					if(useLeftFlank)
					{
						modelProb = SimpleKmerCount.getPredictedCountWithLeftFlank(leftFlank,
								markovModelLength,seq,  marKovProbabilities);
					}
					else
					{
						modelProb = SimpleKmerCount.getPredictedCount(
								markovModelLength, markovTotalCount, seq,
								markovCounts, marKovProbabilities);
					}
					
					double expectedCount = modelProb * kNoOfValidRead;
					results.add(modelProb);
					results.add(expectedCount);
				}
				
				if(directOutputWithoutBuildingArray)
				{
					results.set(1, ((Sequence)results.get(1)).getString());
					Util.output1DArray(writer, results, 1);
				}
				else
				{
					list.add(results);
				}
				
				if(useTop && list.size() > top)
				{
					list.poll();
				}
				
			}
			node.pop();
		}
		
		DebugLog.log("Lines processed:"+ lineNumber);
		
		String[] seqs=new String[list.size()];
		int[] counts=new int[list.size()];
		double[] exps=new double[list.size()];
		double[] probs=new double[list.size()];
		
		int i=0;
		while(!list.isEmpty())
		{
			ArrayList result = (ArrayList) list.poll();
			seqs[i]=   (String)  ((Sequence)result.get(1)).getString();
			counts[i]= (Integer) result.get(2);
			if(useMarkov)
			{
				probs[i] = (Double) result.get(3);
				exps[i] = (Double) result.get(4);
			}
			i++;
		}


		if(directOutputWithoutBuildingArray)
		{
			writer.close();
		}
		
		Object[] data=null;
		if(!useMarkov)
		{
			data=new Object[]{seqs, counts};
			
			if(writer!=null && useTop)
			{
				Util.outputColumnBasedArrays(writer, data, header1);
				writer.close();
				return null;
			}
			
			return data;
		}
		else
		{
			data=new Object[]{seqs, counts,probs, exps};

			if(writer!=null && useTop)
			{
				Util.outputColumnBasedArrays(writer, data,  header2);
				writer.close();
				return null;
			}
			
			return data;
		}
		
		}catch(Exception ex)
		{
			DebugLog.log(ex);
			return null;
		}
	}
	
	/**
	 * Calculates the probability for a particular sequence in alphabets
	 * 
	 * Called in R.
	 * 
	 * @param sequence
	 * @param mm
	 * @return
	 */
	public static double getProbability(String sequence, MarkovModelInfo mm)
	{	
		float[] marKovProbabilities = null;
		String morkovObjPath = null;
		String markovCountPath = null;
		int markovModelLength = 0;
		int[] markovCounts = null;
		long markovTotalCount = 0;
		try
		{
			morkovObjPath = mm.getMarkovObjPath();
			markovCountPath = mm.getMarkovCountsPath();
			markovModelLength = mm.getMarkovLength();
			markovTotalCount = mm.getMarkovLengthTotalCount();
			
			// load markovProbabilities
			DebugLog.log("Reading Markov Prob file:" + morkovObjPath);
			FileInputStream fos = new FileInputStream(morkovObjPath);
			ObjectInputStream oos = new ObjectInputStream(fos);
			marKovProbabilities = (float[]) oos.readObject();

			DebugLog.log("Reading Markov Count file:" + markovCountPath);
			ArraysMergerHeap.MinHeapNode node2 = new ArraysMergerHeap.MinHeapNode(
					markovCountPath);
			CountObject obj = null;
			markovCounts = new int[1 << (2 * (markovModelLength))];

			while ((obj = node2.peek()) != null)
			{
				Sequence seq = (Sequence) obj.getKey();
				markovCounts[(int) seq.getValue()] = obj.getCount();
				node2.pop();
			}
			
		}
		catch(Exception ex)
		{
			DebugLog.log(ex);
			throw new RuntimeException(ex);
		}
		

		int k = sequence.length();
		Sequence seq = new Sequence(sequence,0, k );
		
		double modelProb = SimpleKmerCount.getPredictedCount(
				markovModelLength, markovTotalCount, seq,
				markovCounts, marKovProbabilities);
		return modelProb;
	}
	
	/**
	 * Creates and verifies a ExperimentReference object 
	 * 
	 * Called in R.
	 * 
	 * @param sequencingRunName
	 * @param sampleName
	 * @param round
	 * @return
	 */
	public static ExperimentReference getExperimentReference(String sequencingRunName, 
										String sampleName, Integer round)
	{
		ExperimentReference ref=new ExperimentReference();
		ref.setSequencingName(sequencingRunName);
		ref.setSampleName(sampleName);
		ref.setSampleRound(round);
		
		Sample sample = reader.getSample(ref);
		if(sample==null)
		{
			throw new RuntimeException("Sample ["+getSampleID(ref)+"] not found.");
		}
		return ref;
		
	}
	
	/**
	 * Returns the attributes of a sample
	 * 
	 * Called in R.
	 * 
	 * @param ref
	 * @return
	 */
	public static Object[] getExperimentReferenceInfo(ExperimentReference ref)
	{
		ArrayList<String> attributes =new ArrayList<String>();
		ArrayList<String> values =new ArrayList<String>();
		
		SequencingRunInfo runInfo = reader.getSequencingRunInfo(ref);
		Sample sample = reader.getSample(ref);
		if(sample==null)
		{
			throw new RuntimeException("Sample ["+getSampleID(ref)+"] not found.");
		}
		attributes.add("FilePath");
		values.add(runInfo.getDataFile());
		attributes.add("LeftBarcode");
		values.add(sample.getLeftBarcode());
		attributes.add("RightBarcode");
		values.add(sample.getRightBarcode());
		attributes.add("VariableRegionLength");
		values.add(sample.getVariableRegionLength()+"");
		
		String as[]=new String[attributes.size()];
		String vs[]=new String[attributes.size()];
		attributes.toArray(as);
		values.toArray(vs);
		return new Object[]{as, vs };
	}
	
	public static void removeTempFiles(String tempFolder)
	{
		File f=new File(tempFolder);
		for(File file:f.listFiles())
		{
			if(!file.getName().endsWith(".txt"))
			{
				DebugLog.log("Removing temp file:"+file.getAbsolutePath());
				file.delete();
			}
		}
	}
}
