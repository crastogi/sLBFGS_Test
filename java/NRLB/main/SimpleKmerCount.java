package main;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Properties;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Semaphore;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import config.ExperimentReference;
import config.InputDataSetStats;
import config.SELEXSequencingConfig;
import config.Sample;
import config.SequencingRunInfo;

import base.ArraysMerger;
import base.ArraysMergerHeap;
import base.CountObject;
import base.CountObjectIterator;
import base.DebugLog;
import base.FatalException;
import base.MarkovModelInfo;
import base.MarkovModelOption;
import base.RegexOption;
import base.SELEXConfigReader;
import base.Sequence;
import base.Util;

/**
 * SimpleKmerCount is designed to do efficient k-mer counting, Markov model building, 
 * and information gain calculation.
 */
public class SimpleKmerCount
{
	private String path;
	private int kmerLength;
	private int variableLength;
	private long kmerMask = 0;
	private String outputPath;
	private Boolean saveAsText=false;

	private String leftBarcode;
	private String leftFlank;
	private Boolean leftStartsWith;
	private String rightBarcode;

	private static int NUCLEOTIDE_PER_LONG = 30;
	private boolean externalSorting = false;
	private int MAX_THREAD_NUM = 1;
	private int KMER_PER_THREAD =   2000 * 1000;
	private int EXTERNAL_SORT_SIZE = 2000 * 1000;
//	private int LINES_PER_THREAD = 200 * 1000;

	private int[] counts;
	private ArraysMerger merger;

	private String tempFolder ="./SELEX/";
	private long minimalCount;
	private long maximalCount;
	private long noOfValidRead;
	private Boolean isBinaryFile=false;
	private Boolean isFASTQFile=true; //else FASTA
	private int leftOffset=0;
	private Boolean useSlidingWindow = true;
	private RegexOption regexOption;
	public static Boolean includeReversed = false;
	
	public static final String DATA_FILE_TYPE_TXT="FASTQ.TXT";
	public static final String DATA_FILE_FASTA_TYPE_TXT="FASTA.TXT";
	public static final String DATA_FILE_TYPE_BINARY_CACHE="FASTQ.BINARY.CACHE";
	public static final String DATA_FILE_TYPE_BINARY_COUNT="FASTQ.BINARY.COUNT"; //not used yet

	private Long fastqTotalLineCount;
	private Long fastqTotalReadCount;
	private Long fastqBarcodeMatchedReadCount;
	private Long fastqValidReadCount;
	
	private String fastqCacheFileName;
	private Boolean fastqFileProcessed;
	private ExperimentReference dataSet;
	
	public SimpleKmerCount()
	{
	}
	
	/**
	 * Setup variables related to sample dataset
	 * @param reader
	 * @param dataSet
	 */
	public void initTraining(SELEXConfigReader reader, ExperimentReference dataSet)
	{
		try
		{
			leftOffset= 0;//reset the left offset
			useSlidingWindow=true;//
			this.dataSet = dataSet;
			
			//includeReversed=false;
			
			DebugLog.log("Reading configurtion of [" + SELEX.getSampleID(dataSet) + "]");
			
			SequencingRunInfo runInfo = reader.getSequencingRunInfo(dataSet);
			if(runInfo==null)
			{
				throw new RuntimeException("Can't find SequencingRunInfo["+dataSet.getSequencingName()+"]");
			}
			Sample sample = reader.getSample(dataSet );
			if(sample==null)
			{
				throw new RuntimeException("Can't find Sample["+dataSet.getSampleName()+"]");
			}

			isBinaryFile = false;
			if(DATA_FILE_TYPE_BINARY_CACHE.equals(runInfo.getDataFileType()))
			{
				isBinaryFile=true;
			}
			
			if(DATA_FILE_FASTA_TYPE_TXT.equals(runInfo.getDataFileType()))
			{
				isFASTQFile=false;
			}
			
			path = runInfo.getDataFile();

			//runInfo==null || sample==null
			variableLength = sample.getVariableRegionLength() ; 
			leftBarcode = sample.getLeftBarcode();
			leftFlank = sample.getLeftFlank();
			
			if(reader.getModelSettings()==null)
			{
				leftStartsWith = true;
			}
			else
			{
				leftStartsWith = reader.getModelSettings().getBarcodeOffset()==0; 
			}
			rightBarcode = sample.getRightBarcode() ; 

			trimBarcodes();

			initKMerLength(this.kmerLength);
		}
		catch(Exception ex)
		{
			DebugLog.log(ex); 
			throw new RuntimeException(ex);
		}
	}
	
	public void setTempFolder(String path)
	{
		this.tempFolder = path;
	}
	
	public Boolean getIsBinaryFile()
	{
		return isBinaryFile;
	}

	public Boolean getIsFASTQFile()
	{
		return this.isFASTQFile;
	}

	public RegexOption getRegexOption()
	{
		return regexOption;
	}

	public void setRegexOption(RegexOption regexOption)
	{
		this.regexOption = regexOption;
	}

	private void trimBarcodes()
	{
		if(leftBarcode==null)
		{
			leftBarcode = "";
		}
		else
		{
			leftBarcode =leftBarcode.trim();
		}
		if(rightBarcode==null)
		{
			rightBarcode = "";
		}
		else
		{
			rightBarcode=rightBarcode.trim();
		}

		DebugLog.log("left = " + leftBarcode);
		DebugLog.log("right = " + rightBarcode);
		DebugLog.log("Left Barcode Starts at the beginning = " + leftStartsWith);
	}

	public void initKMerLength(String k)
	{
		 initKMerLength(Integer.parseInt(k));
	}

	public void initKMerLength(Integer k)
	{
		this.minimalCount = Long.MAX_VALUE;
		this.maximalCount = Long.MIN_VALUE;
		this.kmerLength = k;

		DebugLog.log("KMerLength:" + this.kmerLength);

		kmerMask = 0;
		for (int i = 0; i < this.kmerLength; i++)
		{
			kmerMask <<= 2;
			kmerMask |= 3;
		}

		// Switch to external sorting when K >= 14
		externalSorting = kmerLength >= 14;
		//externalSorting = kmerLength >= 16; //for large RAM

		if (externalSorting)
		{
			String tempKey= getFileName()+":"+getFastqKey()+":"+kmerLength;
			merger = new ArraysMerger(tempFolder, "temp-"+Util.getMD5(tempKey),  EXTERNAL_SORT_SIZE);
		} else
		{
			long size = 1 << kmerLength * 2;
			long MB = ((size * 4) / 1000000);
			try
			{
				DebugLog.log("Requesting memory size of " + MB + "MB (" + size + ")elements");
				counts = null;
				counts = new int[(int)size];
				DebugLog.log("Memory allocated.");
			} catch (Exception e)
			{
				throw new FatalException(
						"JVM Memory is not large enough. At least " + MB + " MB needed.");
			}
		}
	}
	
	// convert to Sequence to Longs
	/** Example:
	 *  Sequence:
	 *        G  C  A  A  G  A  G  C  C  C  T  G  C  A  A  
	 *        C  C
	 *  Output Longs:
	 *  0b xx 10 01 00 00 10 00 10 01 01 01 11 10 01 00 00
	 *  0b xx 01 00 00 00 00 00 00 00 00 00 00 00 00 00 00
	 */
	private void cacheOutput(DataOutputStream os, SubString s)
			throws IOException
	{
		int seqIdx = 0;
		long tempResult = 0;
		for (int i = s.startIdx; i < s.endIdx; i++, seqIdx++)
		{
			char c = s.base.charAt(i);
//			System.out.print(c);
			long n = Sequence.getCharCode(c);
			tempResult = (tempResult << 2 | n);

			if (seqIdx == NUCLEOTIDE_PER_LONG - 1)
			{
//				System.out.println(" Output:"
//						+ new Sequence(tempResult, 15).getString());
				os.writeLong(tempResult);
				seqIdx = -1;
				tempResult = 0;
			}
		}
 
		if (seqIdx != 0)
		{
			//Push the data to the left most position
			//DebugLog.log("Org="+tempResult +"  --> pushed:" + (tempResult << (30-seqIdx*2) ));
			os.writeLong(tempResult << (NUCLEOTIDE_PER_LONG- seqIdx)*2 );
		}

		//System.out.println("");
	}

	public int getLeftOffset()
	{
		return leftOffset;
	}

	public void setLeftOffset(int leftOffset)
	{
		this.leftOffset = leftOffset;
	}

	public void setUseSlidingWindow(Boolean useSlidingWindow)
	{
		this.useSlidingWindow = useSlidingWindow;
	}

	private void cacheOutput(DataOutputStream os, BinarySequence seq)throws IOException
	{
		for(int i=0;i<seq.data.length;i++)
		{
			os.writeLong(seq.data[i]);
		}
	}
	
	public String getCacheFileName()
	{
		return this.tempFolder +getFileName (this.path)+"-"+Util.getMD5(getFastqKey());
	}
	

	private String getFileName()
	{
		return this.path;
	}
	
	public String getFastqKey()
	{
		String key=this.leftBarcode+"-"+ this.rightBarcode+"-"+this.variableLength + "-" + 
		Util.getVariableRegionFilterString(this.regexOption);
		DebugLog.log("Fastqkey: "+ key);
		return key;
	}
	
	private String getFileName(String path)
	{
//		int index = path.lastIndexOf('/');
//		if(index==-1)
//			return path;
//		else
//			return path.substring(index+1);
		return new File(path).getName();
	}
	
	private static Random randomObj = new Random();
	
	private int getRandomIdx(long counts[], long[] avgNum, int k)
	{
		int idx= randomObj.nextInt(k);
		while(counts[idx] >= avgNum[idx])
		{
			idx--;
			if(idx<0)
			{
				idx=k-1;
			}
		}
		return idx;
	}

	/**
	 * Returns a list of newly created data set names and their counts.
	 * The split is done 
	 */
	public Object[] splitCacheFile(SELEXConfigReader reader, ExperimentReference dataSet, Properties newProp, 
		    String prefix, Long totalReads, Double[] proportions, String outputFolder,
			List<InputDataSetStats> ds)
	{
		String[] newDataSets1 = new String[]{};
		String[] newDataSets2 = new String[]{};
		Integer[] newDataSets3 = new Integer[]{};
		int[] rounds = new int[]{};
		long counts[] = new long[]{};
		try
		{
			 String binaryFileCache = getCacheFileName();
			 if(!new File(binaryFileCache).exists())
			 {
				 throw new RuntimeException("Can't find binary cache file to split:"+binaryFileCache);
			 }
			
			DebugLog.log("Opening cache file:" +binaryFileCache);
			DataInputStream cacheInput = new DataInputStream(new
					 				BufferedInputStream( new FileInputStream(binaryFileCache)));

			int len=(int)Math.ceil(this.variableLength*1.0/NUCLEOTIDE_PER_LONG);
			DebugLog.log("Length of each sequence:"+len);
			int NUM_SPLITS = proportions.length;
			DebugLog.log("Spliting into ["+NUM_SPLITS+"] parts.");
			
			//
			SequencingRunInfo runInfo = reader.getSequencingRunInfo(dataSet);
			Sample sample = reader.getSample(dataSet );
			//output files
			DataOutputStream os[]=new DataOutputStream[NUM_SPLITS];
			counts =new long[NUM_SPLITS];
			long capCounts[] = new long[NUM_SPLITS];
			newDataSets1 = new String[NUM_SPLITS];
			newDataSets2 = new String[NUM_SPLITS];
			newDataSets3 = new Integer[NUM_SPLITS];
			rounds = new int[NUM_SPLITS];
			DebugLog.log("Total Size:"+ totalReads);
			for(int i=0;i<NUM_SPLITS;i++)
			{
				counts[i] =0L;
				capCounts[i] = (long) Math.ceil( (totalReads*proportions[i]));
				DebugLog.log("Bucket Cap["+(i)+"] Size:"+ capCounts[i]);
				String path = outputFolder+ "/"+(i)+".part" ;
				String xmlPath = outputFolder+ "/"+(i)+".xml" ;
				String newDataSetIdx = sample.getName()+"."+ prefix+"."+(i+1);
				
				//copyDataSetProperties(reader, dataSet, newProp , newDataSetIdx , path, true);
				os[i] = new DataOutputStream( new BufferedOutputStream(new FileOutputStream(path),1024*1024*5));
				DebugLog.log("Outputting to : "+ path);
				//
				SELEXSequencingConfig config=new SELEXSequencingConfig();
				SequencingRunInfo clonedRunInfo = Util.cloneWithoutSample(runInfo);
				clonedRunInfo.setDataFile(path);
				clonedRunInfo.setName(clonedRunInfo.getName()+"."+ prefix+"."+(i+1));
				Sample clonedSample = Util.cloneSample(sample);
				clonedSample.setName(newDataSetIdx);
				clonedRunInfo.getSample().add(clonedSample);
				clonedRunInfo.setDataFileType(DATA_FILE_TYPE_BINARY_CACHE);
				config.getSequencingRunInfo().add(clonedRunInfo);
				SELEXConfigReader.writeConfig(xmlPath,config);
				DebugLog.log("XML Config : "+ xmlPath);
				newProp.setProperty(newDataSetIdx, xmlPath);
				//
				newDataSets1[i]= (clonedRunInfo.getName());
				newDataSets2[i]= (clonedSample.getName());
				newDataSets3[i]= (clonedSample.getRound());

				ExperimentReference ref=new ExperimentReference();
				ref.setSequencingName(clonedRunInfo.getName());
				ref.setSampleName(clonedSample.getName());
				ref.setSampleRound(clonedSample.getRound());
				InputDataSetStats dataSetStats=new InputDataSetStats();
				dataSetStats.setId(SELEX.getSampleID(ref));
				dataSetStats.setOriginalDataFilePath(clonedRunInfo.getDataFile());
				dataSetStats.setCachedDataFilePath(clonedRunInfo.getDataFileType());
				dataSetStats.setExperimentReference(dataSet);
				dataSetStats.setLeftBarcode(clonedSample.getLeftBarcode());
				dataSetStats.setRightBarcode(clonedSample.getRightBarcode());
				
				ds.add(dataSetStats);
			}

			int q=0;
			while(true)
			{
				BinarySequence binSeq=new BinarySequence(len);
//				if(cacheInput.available()==0)
//				{
//					break;
//				}
				for(int k=0;k<len;k++)
				{
					long l=cacheInput.readLong();
					binSeq.write(l);
				}
				
				q++;
				int idx= getRandomIdx(counts , capCounts, NUM_SPLITS);
				counts[idx]++;
				
				cacheOutput(os[idx], binSeq);
				
				if(q==totalReads)
				{
					break;
				}
			}

			for(int i=0;i<NUM_SPLITS;i++)
			{
				os[i].flush();
				DebugLog.log("Bucket["+(i)+"] Size:"+ counts[i]);
				
				InputDataSetStats dataSetStats=ds.get(i);
				dataSetStats.setTotalLineCount(counts[i]);
				dataSetStats.setTotalReadCount(counts[i]);
				dataSetStats.setTotalBarcodeMatchedReadCount(counts[i]);
				dataSetStats.setValidReadCount(counts[i]);
			}
			
		}catch(Exception ex)
		{
			DebugLog.log(ex);
			throw new RuntimeException(ex);
		}
		
		for(int i=0;i<newDataSets3.length;i++)
		{
			rounds[i]=newDataSets3[i];
		}
		return new Object[]{newDataSets1,newDataSets2,rounds,counts };

	}
	
	public void processCacheFile(String binaryFileCache) throws Exception
	{
		DebugLog.log("Opening cache file:" +binaryFileCache);
		DataInputStream cacheInput = new DataInputStream(new
				 				BufferedInputStream( new FileInputStream(binaryFileCache)));
		 
		final ExecutorService service = Executors.newFixedThreadPool(MAX_THREAD_NUM);
		final Semaphore semaphore = new Semaphore(MAX_THREAD_NUM + 1);

		int totalReads = 0;
		int numberOfValidSeq = 0;
		int numberOfValidSeqTemp = 0;
		List<BinarySequence> lines = new LinkedList<BinarySequence>();

		DebugLog.log("Reading...");
		int len=(int)Math.ceil(this.variableLength*1.0/NUCLEOTIDE_PER_LONG);
		DebugLog.log("Length of each sequence:"+len);
		
		int kmerMultiplier = this.variableLength - this.kmerLength + 1;
		int LINES_PER_THREAD = KMER_PER_THREAD/kmerMultiplier;
		
		DebugLog.log("Reading...");
		boolean doneReading=false;
		while(!doneReading)
		{
			doneReading = (cacheInput.available()==0);
			if(!doneReading)
			{
				BinarySequence binSeq=new BinarySequence(len);
				for(int k=0;k<len;k++)
				{
					long l=cacheInput.readLong();
					binSeq.write(l);
				}
				
				lines.add(binSeq);
				
				totalReads++;
				
				numberOfValidSeq++;
				numberOfValidSeqTemp++;
			}
			
			if (numberOfValidSeqTemp == LINES_PER_THREAD || doneReading)
			{
				final List<BinarySequence> linesLocal = lines;
				lines = new LinkedList<BinarySequence>();
				semaphore.acquire();
				numberOfValidSeqTemp = 0;

//				DebugLog.log("Submitting a job: ");
				service.submit(new Runnable()
				{
					public void run()
					{
						DebugLog.log("Processsing....trunk of size:"
								+ linesLocal.size());
						try
						{
							doProcessBinarySeqence(linesLocal,leftOffset);
						} catch (Throwable ex)
						{
							DebugLog.log(ex);
							service.shutdown();
							throw new RuntimeException(ex);
						}
						finally
						{
							semaphore.release();
						}
					}
				});

			}

			if (numberOfValidSeq % 1000000 == 0)
			{
				DebugLog.log("Read Count:" + numberOfValidSeq);
				Util.printMemeoryUsage();
			}
		}

		service.shutdown();
		// now wait for the jobs to finish
		boolean result  = service.awaitTermination(1000, TimeUnit.SECONDS);
		DebugLog.log("service.awaitTermination result:" + result);

		DebugLog.log("Total sequence lines:" + totalReads);
		DebugLog.log("Valid k-mer sequences:" + numberOfValidSeq);

		lines = null; //release
		try
		{
			save(this.saveAsText);
		} catch (Exception e)
		{
			DebugLog.log(e);
		}

	}
	
	public void setSaveAsText(Boolean val)
	{
		this.saveAsText=val;
	}
	
	public final String INCOMPLETE_FLAG = "-incomplete";
	
//	private Pattern excludePattern = Pattern.compile("TGA[T|C].{5,6}TGA[T|C]|TGA[T|C].{5,6}[A|G]TCA|[A|G]TCA.{5,6}TGA[T|C]|[A|G]TCA.{5,6}[A|G]TCA|TAATTA.{1,4}TAATTA");

	public void process() throws Exception
	{
		DebugLog.log("Input file["+getFileName()+"] is binary["+isBinaryFile+"] isFASTQFile["
				+isFASTQFile+"] LeftOffset["+this.leftOffset+"]");

		if(!this.isFASTQFile)
		{
			processFASTAFile();
			return;
		}

		//for a fastq file
		if(this.kmerLength > this.variableLength)
		{
			throw new RuntimeException("Kmer count length["+this.kmerLength 
					+"] is greater than variable region length["+this.variableLength+"].");
		}

		DebugLog.log("UseSlidingWindow = " + this.useSlidingWindow);
		DebugLog.log("LeftOffSet       = " + this.leftOffset);
		DebugLog.log("IncludeReversed  = " + this.includeReversed);
		
		if(isBinaryFile)
		{
			 fastqFileProcessed = false;
			 processCacheFile(getFileName());
			 return;
		}
		
		 //create cache file here
		 String binaryFileCache = getCacheFileName();
		 fastqCacheFileName = binaryFileCache;
		 fastqFileProcessed = true;
		 
		 if(new File(binaryFileCache).exists())
		 {
			 fastqFileProcessed = false;
			 processCacheFile(binaryFileCache);
			 return;
		 }
		 else
		 {
			 binaryFileCache += INCOMPLETE_FLAG;
		 }
		 
		 DebugLog.log("Cache to file:" +binaryFileCache);
		 DataOutputStream cacheOutput = new DataOutputStream(new
				 				BufferedOutputStream( new FileOutputStream(binaryFileCache),1000*1000*10));
		 
		BufferedReader br = null;
		String path=getFileName();
		if (path.endsWith(".gz"))
		{
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(
					new FileInputStream(path))));
		} else
		{
			br = new BufferedReader(new InputStreamReader(new FileInputStream(
					path)));
		}


		int kmerMultiplier = this.variableLength - this.kmerLength + 1;
		int LINES_PER_THREAD = KMER_PER_THREAD/kmerMultiplier;
		
		final ExecutorService service = Executors.newFixedThreadPool(MAX_THREAD_NUM);
		final Semaphore semaphore = new Semaphore(MAX_THREAD_NUM * 2);

		long lineCount = 0;
		long totalReads = 0;
		String line = null;
		long numberOfValidSeq = 0;
		long numberOfValidSeqTemp = 0;
		long numberOfBarcodeMatched = 0;
		long numberOfFiltered = 0;
		List<SubString> lines = new LinkedList<SubString>();
		
		Pattern excludePattern =null;
		Pattern includePattern =null;
		Pattern groupPattern = null;
		if(regexOption!=null)
		{
			excludePattern=regexOption.getVariableRegionExcludePattern();
			includePattern=regexOption.getVariableRegionIncludePattern();
			groupPattern = regexOption.getVariableRegionGroupPattern();
		}

		DebugLog.log("Reading...");
		boolean doneReading=false;
		while(!doneReading)
		{
			doneReading = ((line = br.readLine()) == null);
			if(!doneReading)
			{
				//variable region processing
				lineCount++;

				if (lineCount % 4 != 2)
					continue;
				boolean valid = true;

				totalReads++;
				
				int idx1 = -1;
				if (leftStartsWith)
				{
					if (line.startsWith(this.leftBarcode))
						idx1 = 0;
				} else
				{
					idx1 = line.indexOf(this.leftBarcode);
				}

				// DebugLog.log(idx1);
				// DebugLog.log((idx1!=-1));
				int rightIdx = idx1 + this.leftBarcode.length()
						+ this.variableLength;
				if (idx1 != -1
						&& line.regionMatches(rightIdx, this.rightBarcode, 0,
								this.rightBarcode.length()))
				{
					// DebugLog.log("Matched:"+line);
					numberOfBarcodeMatched++;

					int leftIdx = idx1 + this.leftBarcode.length();
					// DebugLog.log(line);
					for (int i = leftIdx; i < rightIdx; i++)
					{
						if (!Sequence.isValidSymbol(line.charAt(i)))
						{
							valid = false;
							break;
						}
					}
					
					//apply filters here
					if(valid)
					{
						String str=null;
						
						if(groupPattern!=null)
						{
							str=line.substring(leftIdx,rightIdx);
							valid = false;
							Matcher m = groupPattern.matcher(str);
							if(m.find() && m.groupCount()>0)
							{
								rightIdx = leftIdx + m.end(1);
								leftIdx =  leftIdx + m.start(1);
								valid = true;
							}
						}
						else
						{
							if(includePattern!=null)
							{
								str=line.substring(leftIdx,rightIdx);
								if(!includePattern.matcher(str).find())
								{
									valid = false;
									//DebugLog.log("Filtered ( not matched by the includePattern ): "+ str);
								}
							}
							else if(excludePattern!=null)
							{
								str=line.substring(leftIdx,rightIdx);
								if(excludePattern.matcher(str).find())
								{
									valid = false;
									//DebugLog.log("Filtered ( matched by the excludePattern ): "+ str);
								}
							}
						}
						
						
						if(!valid)
						{
							numberOfFiltered ++;
						}
					}
					
					if (valid)
					{
						numberOfValidSeq++;
						numberOfValidSeqTemp++;
						
						// handle valid reads here
						// DebugLog.log("Matched:"+line);
//						DebugLog.log("Matched:"+line.substring(leftIdx,
//								 rightIdx));

						SubString subStr = new SubString(line, leftIdx, rightIdx);
						lines.add(subStr);
						cacheOutput(cacheOutput, subStr);
					}

				}
			}
			
			if (numberOfValidSeqTemp == LINES_PER_THREAD || doneReading)
			{
				DebugLog.log("Read Count:" + numberOfValidSeq);
				Util.printMemeoryUsage();
				
				final List<SubString> linesLocal = lines;
				lines = new LinkedList<SubString>();
				semaphore.acquire();
				numberOfValidSeqTemp = 0;

//				DebugLog.log("Submitting a job: ");
				service.submit(new Runnable()
				{
					public void run()
					{
						DebugLog.log("Processsing....trunk of size:"
								+ linesLocal.size());
						try
						{
							doProcess(linesLocal,leftOffset);
						} catch (Throwable ex)
						{
							DebugLog.log(ex);
							service.shutdown();
							throw new RuntimeException(ex);
						}
						finally
						{
							semaphore.release();
						}
					}
				});

			}

		}

		service.shutdown();
		// now wait for the jobs to finish
		boolean result = service.awaitTermination(1000, TimeUnit.SECONDS);
		DebugLog.log("service.awaitTermination result:" + result);
		
		//Thread.sleep(1000*60);
		
		while(!service.isTerminated())
		{
			DebugLog.log("waiting...");
			Thread.sleep(1000);
		}

		lines = null; //release memory
		
		DebugLog.log("Total lines:" + lineCount);
		DebugLog.log("Total sequence lines:" + totalReads);
		DebugLog.log("Total barcode-matched lines:" + numberOfBarcodeMatched);
		DebugLog.log("Total filtered lines:" + numberOfFiltered);
		DebugLog.log("Valid k-mer sequences:" + numberOfValidSeq);

		fastqTotalLineCount = lineCount;
		fastqTotalReadCount = totalReads;
		fastqBarcodeMatchedReadCount = numberOfBarcodeMatched ;
		fastqValidReadCount = numberOfValidSeq;
		
		// this.noOfValidRead = numberOfValidSeq;
		// save output
		try
		{
			cacheOutput.close();
			String cleanBinaryFileCache = binaryFileCache.substring(0,binaryFileCache.length() - INCOMPLETE_FLAG.length());
			DebugLog.log("Rename cache file ["+binaryFileCache+"] to ["+cleanBinaryFileCache+"].");
			new File(binaryFileCache).renameTo( new File(cleanBinaryFileCache) );
			save(this.saveAsText);
		} catch (Exception e)
		{
			// TODO Auto-generated catch block
			DebugLog.log(e);
		}

	}
	
	public void processFASTAFile() throws Exception
	{
		DebugLog.log("IncludeReversed  = " + this.includeReversed);
		
		fastqFileProcessed = true;
		 
		BufferedReader br = null;
		String path=getFileName();
		if (path.endsWith(".gz"))
		{
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(
					new FileInputStream(path))));
		} else
		{
			br = new BufferedReader(new InputStreamReader(new FileInputStream(
					path)));
		}

		//KMER_PER_THREAD

		long charCounts = 0;
		long totalReads = 0;
		int character = -1;
		long numberOfValidSeq = 0;
		long numberOfValidSeqTemp = 0;
		long numberOfBarcodeMatched = 0;
		long numberOfFiltered = 0;
		List<SubString> lines = new LinkedList<SubString>();
		
		DebugLog.log("Reading...");
		DebugLog.log("Header Info:\t"+ br.readLine());
		
		long numInvalidChar = 0;
		
		HashSet<Integer> invalidSymbols =new HashSet<Integer>();
		boolean doneReading=false;
		Sequence lastSequence = null;
		StringBuffer sb=new StringBuffer();
		
		while(!doneReading)
		{
			doneReading = ((character = br.read()) == -1);
			char c = (char)character;
			//DebugLog.log(c);
			if(!Sequence.isValidSymbolEnhanced(c))
			{
				numInvalidChar++;
				invalidSymbols.add((int)c);
			}
			else
			{
				if(lastSequence==null)
				{
					sb.append(c);
					if(sb.length() == this.kmerLength)
					{
						lastSequence = new Sequence(sb.toString(), 0, this.kmerLength) ;
					}
					else
					{
						continue;
					}
				}
				else
				{
					long v = lastSequence.getValue(); // result[i - 1].getValue();
					long currentChar = Sequence.getCharCode(c);
					long newValue = (v << 2 | currentChar) & kmerMask;
					lastSequence = new Sequence(newValue, this.kmerLength);
				}
				
				if(this.regexOption!=null)
				{
					boolean toInclude=true;
					if(this.regexOption.getKmerIncludePattern()!=null )
					{
						toInclude = this.regexOption.getKmerIncludePattern().matcher(lastSequence.getString()).find() ;
					}
					if(this.regexOption.getKmerExcludePattern()!=null )
					{
						toInclude = !this.regexOption.getKmerExcludePattern().matcher(lastSequence.getString()).find() ;
					}
					if(this.regexOption.getKmerIncludeOnlyPattern()!=null )
					{
						toInclude = this.regexOption.getKmerIncludeOnlyPattern().matcher(lastSequence.getString()).find() ;
					}
					if(!toInclude)
					{
						continue;
					}
				}
				
				numberOfValidSeq++;
				if (this.externalSorting)
				{
					merger.add(new CountObject(lastSequence, 1));
					if(includeReversed)
					{
						numberOfValidSeq++;
						merger.add(new CountObject(lastSequence.getReverseComplement(), 1));
					}
				} else
				{
					updateCounts(lastSequence, 1);
					if(includeReversed)
					{
						numberOfValidSeq++;
						updateCounts(lastSequence.getReverseComplement(), 1);
					}
				}
				
			}
			
			if(!doneReading)
			{
				//variable region processing
				charCounts++;
			}
			
		
		}

		DebugLog.log("Total charCounts: "+charCounts);
		DebugLog.log("numInvalidChar: "+numInvalidChar);
		DebugLog.log("invalidSymbols: "+invalidSymbols);

		DebugLog.log("Valid k-mer sequences:" + numberOfValidSeq);

		fastqTotalLineCount = 2L;
		fastqTotalReadCount = 2L;
		fastqBarcodeMatchedReadCount = 2L ;
		fastqValidReadCount = numberOfValidSeq;
		fastqCacheFileName = "NA";
		
		// this.noOfValidRead = numberOfValidSeq;
		// save output
		try
		{
			save(this.saveAsText);
		} catch (Exception e)
		{
			// TODO Auto-generated catch block
			DebugLog.log(e);
		}

	}
	
	

	/**
	 * A helper class to manage a String and its indices
	 */
	class SubString
	{
		String base;
		int startIdx;
		int endIdx;

		public SubString(String base, int startIdx, int endIdx)
		{
			super();
			this.base = base;
			this.startIdx = startIdx;
			this.endIdx = endIdx;
		}
	}
	
	/**
	 * A helper class to manage a String and its indices
	 */
	class BinarySequence
	{
		long data[];
		int writeIndex=0;
		int readIndex=0;
		
		long readMask = 0;
		int shiftOffset = 0;

		public BinarySequence(int len)
		{
			data=new long[len];
			writeIndex=0;
			resetMask();
		}
		
		public void write(long d)
		{
			data[writeIndex++]=d;
		}
		
		private void resetMask()
		{
			shiftOffset = (NUCLEOTIDE_PER_LONG - 1)*2 ; //28;
			readMask =  3L <<shiftOffset;
		}
		
		private void next()
		{
			shiftOffset-=2;
			readMask >>= 2;
			if(shiftOffset<0)
			{
				resetMask();
				readIndex++;
			}
		}
		
		public long readChar()
		{
			long val = (data[readIndex] & readMask) >> shiftOffset; 
			next();
			return val;
		}
		
	}
	
	private void doProcess(List<SubString> lines,int leftOffset)
	{
		//DebugLog.log("Thread #"+slotIdx+" started.");
		DebugLog.log("doProcessing:"+lines.size());
		int i = 0;
		for (SubString line : lines)
		{
			i++;
			try
			{
				if (this.externalSorting)
				{
					for (Sequence subseq : slidingWindow(line.base,
							line.startIdx, line.endIdx,leftOffset, useSlidingWindow))
					{
						merger.add(new CountObject(subseq, 1));
						if(includeReversed)
						{
							merger.add(new CountObject(subseq.getReverseComplement(), 1));
						}
					}
				} else
				{
					for (Sequence subseq : slidingWindow(line.base,
							line.startIdx, line.endIdx,leftOffset, useSlidingWindow))
					{
						updateCounts(subseq, 1);
						if(includeReversed)
						{
							updateCounts(subseq.getReverseComplement(), 1);
						}
					}
				}

//				if (i % 100000 == 0)
//				{
//					DebugLog.log((i * 100.0 / size) + " % ");
//				}
			} catch (FatalException e)
			{
				DebugLog.log(e);
				throw e;
			} catch (OutOfMemoryError e)
			{
				DebugLog.log(e);
				throw e;
			} catch (Throwable ex)
			{
				System.err.println("Error processing:" + line);
				DebugLog.log(ex);
			}
		}
	}
	
	private void doProcessBinarySeqence(List<BinarySequence> lines,int leftOffset)
	{
		// DebugLog.log("Thread #"+slotIdx+" started.");
//		DebugLog.log("doProcessing:");
		int size = lines.size();
		int i = 0;
		for (BinarySequence line : lines)
		{
			i++;
			try
			{
				// String[] tokens=line.split("\\s+");
				// Integer c=Integer.parseInt(tokens[1]);

				if (this.externalSorting)
				{
					for (Sequence subseq : slidingWindowBinary(line,leftOffset,useSlidingWindow))
					{
						merger.add(new CountObject(subseq, 1));
						if(includeReversed)
						{
							merger.add(new CountObject(subseq.getReverseComplement(), 1));
						}
					}
				} else
				{
					for (Sequence subseq : slidingWindowBinary(line,leftOffset,useSlidingWindow))
					{
						// System.out.println(subseq.getString());
						updateCounts(subseq, 1);
						if(includeReversed)
						{
							updateCounts(subseq.getReverseComplement(), 1);
						}
					}
				}

//				if (i % 100000 == 0)
//				{
//					DebugLog.log((i * 100.0 / size) + " % ");
//				}
			} catch (FatalException e)
			{
				DebugLog.log(e);
				throw e;
			} catch (OutOfMemoryError e)
			{
				DebugLog.log(e);
				throw e;
			} catch (Throwable ex)
			{
				System.err.println("Error processing:" + line);
				DebugLog.log(ex);
			}
		}

		// System.out.println("Thread #"+slotIdx+" ended.");
	}
	
	public List<Sequence> slidingWindowBinary(BinarySequence binSeq,int leftOffset, boolean useSlidingWindow)
	{
		int len = this.variableLength;
		int c = this.kmerLength;
		// size of output items
		int size = 1;
		if(useSlidingWindow)
		{
			size = len - c + 1 - leftOffset;
		}

		// DebugLog.log("slidingWindow:"+base.substring(start, end));
		if (size <= 0)
		{
			// DebugLog.log(subString);
			throw new RuntimeException("Raw input [" + binSeq + "](length:"
					+ len + ") not long enough for required " + c
					+ "-mer counting.");
		}
		LinkedList<Sequence> result = new LinkedList<Sequence>();
		StringBuffer sb=new StringBuffer(c);
		
		for (int i = 0; i < size; i++)
		{
			long currentChar = 0;
			if(i==0)
			{
				long val=0;
				for(int j=0;j<leftOffset;j++) //skip a few alphabets
				{
					binSeq.readChar();
				}
				for(int j=0;j<c;j++)
				{
					val =  (val <<2) | binSeq.readChar();
				}
				result.add( new Sequence(val, c) );
			}
			else
			{
				long v = result.peekLast().getValue(); //result[i - 1].getValue();
				currentChar = binSeq.readChar();
				long newValue = (v << 2 | currentChar) & kmerMask;
				result.add( new Sequence(newValue, c) );
			}
			
			if(this.regexOption!=null)
			{
				if(i==0)
				{
					sb.append( result.peekLast() );
				}
				else
				{
					for(int j=sb.length()-1;j>0;j--)
					{
						sb.setCharAt(j,sb.charAt(j-1));
					}
					sb.setCharAt(0, Sequence.getChar((int)currentChar));
				}

				boolean include=true;
				if(this.regexOption.getKmerIncludePattern()!=null )
				{
					include = this.regexOption.getKmerIncludePattern().matcher(sb).find() ;
				}
				if(this.regexOption.getKmerExcludePattern()!=null )
				{
					include = !this.regexOption.getKmerExcludePattern().matcher(sb).find() ;
				}
				if(this.regexOption.getKmerIncludeOnlyPattern()!=null )
				{
					include = this.regexOption.getKmerIncludeOnlyPattern().matcher(sb).find() ;
				}
				if(!include)
				{
					result.pollLast();
				}
			}
			//System.out.println(result[i].getString());
			
		}
		//System.out.println(Arrays.toString(result));
		return result;
	}

	public List<Sequence> slidingWindow(String base, int start, int end,int leftOffset,boolean useSlidingWindow)
	{
		int len = end - start;
		int c = this.kmerLength;
		// size of output items
		int size = len - c + 1  - leftOffset; //of sliding window
		
		// DebugLog.log("slidingWindow:"+base.substring(start, end));
		if (size <= 0)
		{
			String subString = base.substring(start, end);
			// DebugLog.log(subString);
			throw new FatalException("Raw input [" + subString + "](length:"
					+ len + ") not long enough for required " + c
					+ "-mer counting. Offset=["+leftOffset+"] UseSlidingWindow=["+ useSlidingWindow+ "]");
		}
		
		if(!useSlidingWindow)
		{
			size=1;
		}
		
		StringBuffer sb=new StringBuffer();
		LinkedList<Sequence> result = new LinkedList<Sequence>();
		long currentChar=0;
		
		for (int i = 0; i < size; i++)
		{
			// result[i]=new Sequence(seq, i, c);

			if (i == 0)
			{
				result.add( new Sequence(base, start + i + leftOffset, c) );
			} else
			{
				long v = result.peekLast().getValue(); // result[i - 1].getValue();
				currentChar = Sequence.getCharCode(base.charAt(start + i + c - 1 + leftOffset));
				long newValue = (v << 2 | currentChar) & kmerMask;
				result.add( new Sequence(newValue, c) );
			}
			if(this.regexOption!=null)
			{
				if(i==0)
				{
					sb.append( result.peekLast() );
				}
				else
				{
					for(int j=sb.length()-1;j>0;j--)
					{
						sb.setCharAt(j,sb.charAt(j-1));
					}
					sb.setCharAt(0, Sequence.getChar((int)currentChar));
				}
				
				boolean include=true;
				if(this.regexOption.getKmerIncludePattern()!=null )
				{
					include = this.regexOption.getKmerIncludePattern().matcher(sb).find() ;
				}
				if(this.regexOption.getKmerExcludePattern()!=null )
				{
					include = !this.regexOption.getKmerExcludePattern().matcher(sb).find() ;
				}
				if(this.regexOption.getKmerIncludeOnlyPattern()!=null )
				{
					include = this.regexOption.getKmerIncludeOnlyPattern().matcher(sb).find() ;
				}
				if(!include)
				{
					result.pollLast();
				}
			}
		}
		return result;
	}

	public String setOutputPath(String outputPath)
	{
		this.outputPath = outputPath;
	
		if(this.getRegexOption()!=null)
		{
			this.outputPath += "_"+Util.getMD5( this.getRegexOption().getVariableRegionRegexFormattedString() );
		}
		
		return this.outputPath;
	}

	public void save(boolean asText) throws Exception
	{
		this.noOfValidRead = 0;
		if (this.externalSorting)
		{
			merger.finish();
			// merger.output(this.outputText, this.outputPath);
			Properties props = merger.output(asText, this.outputPath);
			DebugLog.log(props);
			this.noOfValidRead = (Long) props.get("noOfValidRead");
			this.minimalCount = (Long) props.get("lowestCount");
			this.maximalCount = (Long) props.get("highestCount");
		} else
		{
			DebugLog.log("Outputting to file:" + this.outputPath
					+ " [textOutput:" + asText + "]");

			PrintWriter writer = null;
			DataOutputStream out = null;
			if (asText)
			{
				writer = new PrintWriter(this.outputPath);
			} else
			{
				out = new DataOutputStream(new BufferedOutputStream(
						new FileOutputStream(this.outputPath),1000*1000*10));
			}

			int j = 0;
			
			CountObjectIterator iterator =new CountObjectIterator(this.counts,this.kmerLength);
			while(iterator.hasNext())
			{
				CountObject obj=iterator.next();
				
				if(obj==null)
				{
					break;
				}
				if (asText)
				{
					writer.println(obj.getKey() + " " + obj.getCount());
				} else
				{
					CountObject.steamOutSequence(out, obj);
				}

				this.noOfValidRead += obj.getCount();
				if (this.minimalCount > obj.getCount())
				{
					this.minimalCount = obj.getCount();
				}
				if (this.maximalCount < obj.getCount())
				{
					this.maximalCount = obj.getCount();
				}
				j++;
			}

			if (asText)
			{
				writer.close();
			} else
			{
				CountObject.steamOutSequenceClose(out);
				out.flush();
			}
			DebugLog.log("Number of Unique read:" + j);
			DebugLog.log("Closed file:" + this.outputPath);

		}

		Util.printMemeoryUsage();
	}

	/**
	 * This method is not used. Pearson Coefficient is used instead.
	 * @param kmax
	 * @param testTotalCount
	 * @return
	 * @throws Exception
	 */
	/*
	public double crossValidateRSquared(int kmax, int testTotalCount)
			throws Exception
	{
		if (this.externalSorting)
		{
			throw new RuntimeException("Not supported yet.");
		} else
		{
			int k = this.kmerLength;
			ArraysMergerHeap.MinHeapNode node = new ArraysMergerHeap.MinHeapNode(
					this.tempFolder+"/kmax-" + kmax + ".dat");
			CountObject obj = null;

			String statOuptut = this.tempFolder+"/model-" + k + ".dat";
			PrintWriter writer = new PrintWriter(statOuptut);
			DebugLog.log("Outputinig stats to " + statOuptut);
			writer.println("sequence\tobserved\tpredicted\tprob");

			long N = 0;
			double YObservedSum = 0;

			while ((obj = node.peek()) != null)
			{
				N++;
				YObservedSum += obj.getCount();
				// next one
				node.pop();
			}

			double YObservedAvg = YObservedSum / N;
			double SSE = 0;
			double SST = 0;

			node.reset();
			while ((obj = node.peek()) != null)
			{
				double prob = getPredictedCountInternal(kmax,
						(Sequence) obj.getKey());
				long expectedCount = (long) (prob * testTotalCount);
				writer.println(((Sequence) obj.getKey()).getString() + "\t"
						+ obj.getCount() + "\t" + (prob * testTotalCount)
						+ "\t" + prob);
				// update stats
				SSE += Math.pow(obj.getCount() - expectedCount, 2);
				SST += Math.pow(obj.getCount() - YObservedAvg, 2);
				// next one
				node.pop();
			}
			writer.close();
			DebugLog.log("File closed: " + statOuptut);

			DebugLog.log("Model Length : " + this.kmerLength);
			DebugLog.log("N : " + N);
			DebugLog.log("SSE: " + SSE);
			DebugLog.log("SST: " + SST);
			DebugLog.log("YObservedAvg : " + YObservedAvg);
			double R2 = 1 - (SSE / SST) * ((N - 1) / (N - 2));
			DebugLog.log("R squared : " + R2);

			return R2;
		}
	}
	*/

	/**
	 * This method is not used. A division based method is used instead.
	 * @param kmax
	 * @param seq
	 * @return
	 */
	@SuppressWarnings("unused")
	private double getPredictedCountEfficient(int kmax, Sequence seq)
	{
		int size = kmax - this.kmerLength + 1;
		Long mask = (1 << this.kmerLength * 2) - 1L;
		int leftIdx = kmax - this.kmerLength;
		mask <<= leftIdx * 2;
		long val = seq.getValue();
		long seqStub = ((val & mask) >> ((leftIdx + 1) * 2));
		double prob = this.probabilities2[(int) seqStub];
		// System.out.print(new Sequence(seqStub,this.kmerLength-1)+ ":");
		for (int i = 0; i < size; i++)
		{
			long subSeq = (val & mask) >> (leftIdx - i) * 2;
			// System.out.print(new Sequence(subSeq,this.kmerLength)+ " ");
			prob *= this.probabilities[(int) subSeq];
			mask >>= 2;
		}
		return prob;
	}

	public static double getPredictedCountWithLeftFlank(String leftFlank, 
			int modelLength, Sequence seq,
			float[] probabilitiesTable)
	{
		//seq = new Sequence("TCGTACGT");
		
		int predictLength = seq.getLength();
		
		int size = predictLength - modelLength + 1;
		long mask = (1 << modelLength * 2) - 1L;
		int leftIdx = predictLength - modelLength;
		mask <<= leftIdx * 2;
		long val = seq.getValue();
		int order = modelLength-1 ;

		double prob = 1.0;
		if (modelLength != 1) // example:  P(TGG ABCDE) = P(GG) * P(A|GG) * P(B|GA) * P(C|AB) * P(D|BC) * P(E|CD) = P(ABC) * P(D|BC) * P(E|CD)
							  //           P(TGG ABCDE) = 1     * P(A|GG) * P(B|GA) * P(C|AB) * P(D|BC) * P(E|CD) = P(ABC) * P(D|BC) * P(E|CD)
		{
			if(leftFlank.length() < order)
			{
				throw new RuntimeException("Left flank["+leftFlank+"] is too short for calculating Markov model expected counts.");
			}
			
			long leftFlankVal = new Sequence(leftFlank, leftFlank.length() - order, order).getValue();
			long singleCharMask = (3L << (predictLength * 2- 2 )) ;
			long leftFlankMask  = (1  << (2*order+2)) -1 ;
			
			long subSeq = leftFlankVal;
		    //System.out.println(new Sequence(subSeq,order)+ ":");
		    
			for(int i=0;i<=order;i++)//for the first few al
			{ 
				subSeq<<=2;
				subSeq = subSeq & leftFlankMask; 
				subSeq |= ((singleCharMask & val ) >> ((predictLength -1 - i) * 2)) ;
				//System.out.println("xxx "+new Sequence(subSeq,modelLength)+ " ");
				prob *= probabilitiesTable[(int) subSeq]; // P(A|GG), P(B|GA) 
				singleCharMask>>=2;
			}
			
			for (int i = 1; i < size; i++)
			{
				mask >>= 2;
				subSeq = (val & mask) >> (leftIdx - i) * 2; // P(D|BC), and then P(E|CD)
				//System.out.println("--> "+new Sequence(subSeq,modelLength)+ " ");
				prob *= probabilitiesTable[(int) subSeq];
			}
		} 
		else
		{
			for (int i = 0; i < size; i++)
			{
				long subSeq = (val & mask) >> (leftIdx - i) * 2;
				// System.out.print(new Sequence(subSeq,this.kmerLength)+ " ");
				prob *= probabilitiesTable[(int) subSeq];
				mask >>= 2;
			}
		}

		// System.out.println();
		return prob;
	}
	
	public static double getPredictedCount( int modelLength,
			long totalCount, Sequence seq, int[] countsTable,
			float[] probabilitiesTable)
	{
		int predictLength = seq.getLength();
		if(predictLength< modelLength)
		{
			throw new RuntimeException("The sequence["+seq.getString()+"] is too short for the "+(modelLength-1)+"-order Markov model).");
		}
		
		int size = predictLength - modelLength + 1;
		Long mask = (1 << modelLength * 2) - 1L;
		int leftIdx = predictLength - modelLength;
		mask <<= leftIdx * 2;
		long val = seq.getValue();

		double prob = 1.0;
		if (modelLength != 1) // example:  P(ABCDE) = P(AB) * P(C|AB) * P(D|BC) * P(E|CD) = P(ABC) * P(D|BC) * P(E|CD)
		{
			long seqStub = ((val & mask) >> ((leftIdx) * 2)); 
		    //System.out.print(new Sequence(seqStub,modelLength)+ ":");

			prob = 1.0 * countsTable[(int) seqStub] / totalCount; // P(ABC) =  count(ABC) / # 3-mer
			for (int i = 1; i < size; i++)
			{
				mask >>= 2;
				long subSeq = (val & mask) >> (leftIdx - i) * 2; // P(D|BC), and then P(E|CD)
				// System.out.print(new Sequence(subSeq,modelLength)+ " ");
				prob *= probabilitiesTable[(int) subSeq];
			}
		} 
		else
		{
			for (int i = 0; i < size; i++)
			{
				long subSeq = (val & mask) >> (leftIdx - i) * 2;
				// System.out.print(new Sequence(subSeq,this.kmerLength)+ " ");
				prob *= probabilitiesTable[(int) subSeq];
				mask >>= 2;
			}
		}

		// System.out.println();
		return prob;
	}

	public double calculateInformationGain(MarkovModelInfo mm) throws Exception
	{
		int markovModelLength = mm.getMarkovLength();	
		long markovTotalCount= mm.getMarkovLengthTotalCount(); 
		String markovCountPath = mm.getMarkovCountsPath();
		String morkovObjPath = mm.getMarkovObjPath() ;
		
		boolean useLeftFlank = mm.getMarkovModelMethod().toUpperCase().contains(MarkovModelOption.WITH_LEFT_FLANK);
		
			// load markovProbabilities
			DebugLog.log("Reading Markov Prob file:" + morkovObjPath);
			FileInputStream fos = new FileInputStream(morkovObjPath);
			ObjectInputStream oos = new ObjectInputStream(fos);
			float[] marKovProbabilities = (float[]) oos.readObject();

			int[] counts = null;
			DebugLog.log("Reading Markov Count file:" + markovCountPath);
			ArraysMergerHeap.MinHeapNode node = new ArraysMergerHeap.MinHeapNode(
					markovCountPath);
			CountObject obj = null;
			counts = new int[1 << (2 * (markovModelLength))];

			while ((obj = node.peek()) != null)
			{
				Sequence seq = (Sequence) obj.getKey();
				counts[(int) seq.getValue()] = obj.getCount();
				node.pop();
			}
			
			//
			String infoGainCountPath=this.tempFolder +"/infoGain-model-M"+markovModelLength+"-L"+this.kmerLength+".txt";
			PrintWriter writer=new PrintWriter(infoGainCountPath);
			writer.println("Seqence\tExpectedProb\tExpectedCount\tObservedProb\tObservedCount");
			// go through counts
			int j = 0;
			double S100_SUM = 0;
			double S100_MODEL_PROB_SUM = 0;
			double S100_EXP_PROB_SUM = 0;
			boolean cover100percent = true;
			
			//
			CountObjectIterator iterator = null;
			if (this.externalSorting)
			{
				iterator=new CountObjectIterator(this.outputPath);
			}
			else
			{
				iterator=new CountObjectIterator(this.counts,this.kmerLength);
			}

			while(iterator.hasNext())
			{
				CountObject objIn=iterator.next();
				Sequence seq = (Sequence) objIn.getKey();
				//old
				double modelProb = 1 ;
				
				if(useLeftFlank)
				{
					modelProb = getPredictedCountWithLeftFlank(this.leftFlank,
							markovModelLength,seq,  marKovProbabilities);
				}
				else
				{
					modelProb = getPredictedCount( markovModelLength, markovTotalCount, seq,
							counts, marKovProbabilities);
				}

				
				double expectedCount = modelProb * this.noOfValidRead;
				 
				int observedCount = objIn.getCount();
				double observedProb = 1.0 * observedCount / this.noOfValidRead;

				if (observedCount >= 100)
				{
					S100_SUM += observedProb
							* Util.log2(observedProb / modelProb);
					S100_MODEL_PROB_SUM += modelProb;
					S100_EXP_PROB_SUM += observedProb;
					writer.println(seq.getString()+"\t"+modelProb+"\t"+expectedCount+"\t"+observedProb+"\t"+observedCount);
					
				} else
				{
					cover100percent = false;
				}
				//
				j++;
			}
			//
			
			writer.close();
			
			DebugLog.log("noOfValidRead = " + noOfValidRead);

			S100_MODEL_PROB_SUM = Math.min(S100_MODEL_PROB_SUM, 1.0);
			S100_EXP_PROB_SUM = Math.min(S100_EXP_PROB_SUM, 1.0);

			double INFO_GAIN = S100_SUM;

			if (!cover100percent)
			{
				INFO_GAIN += (1 - S100_EXP_PROB_SUM)
						* Util.log2((1 - S100_EXP_PROB_SUM) / (1 - S100_MODEL_PROB_SUM));

			}
			DebugLog.log("Predition Length = " + this.kmerLength);
			DebugLog.log("cover100percent = " + cover100percent);
			DebugLog.log("S100_SUM = " + S100_SUM);
			DebugLog.log("S100_MODEL_PROB_SUM = " + S100_MODEL_PROB_SUM);
			DebugLog.log("S100_EXP_PROB_SUM = " + S100_EXP_PROB_SUM);
			DebugLog.log("INFO_GAIN = " + INFO_GAIN);

			return INFO_GAIN;
	}

	/**
	 * Calculates the Pearson Coefficient
	 * 
	 * @param testCountTablePath 
	 * @param testTotalCount
	 * @param useLeftFlank
	 * @return
	 * @throws Exception
	 */
	public double crossValidatePR(String testCountTablePath, long testTotalCount, boolean useLeftFlank)
			throws Exception
	{
		if (this.externalSorting)
		{
			throw new RuntimeException("Not supported yet.");
		} else
		{
			int k = this.kmerLength;
			ArraysMergerHeap.MinHeapNode node = new ArraysMergerHeap.MinHeapNode(testCountTablePath);
			CountObject testObject = null;

			boolean debug=false;
			PrintWriter writer = null;
			if(debug)
			{
				String statOuptut = this.tempFolder+"/markovModel-" + k + ".txt";
				DebugLog.log("Outputinig stats to " + statOuptut);
				writer = new PrintWriter(statOuptut);
				writer.println("sequence\tobserved\tpredicted\tprob");
			}

			long N = 0;
			double YObservedSum = 0;
			double YPredictedSum = 0;

			while ((testObject = node.peek()) != null)
			{
				Sequence seq = (Sequence) testObject.getKey();
				N++;
				//
				YObservedSum += testObject.getCount();
				//
				double prob = 1;
				
				if(useLeftFlank)
				{
					prob = getPredictedCountWithLeftFlank(this.leftFlank, 
							this.kmerLength, seq,  this.probabilities);
				}
				else
				{
					prob = getPredictedCount(this.kmerLength, this.noOfValidRead,
							seq, this.counts, this.probabilities);
				}
				
				long expectedCount = (long) (prob * testTotalCount);
				YPredictedSum += expectedCount;
				//
				if(debug)
				{
					writer.println(seq.getString() + "\t" + testObject.getCount() + "\t"
							+ (prob * testTotalCount) + "\t" + prob);
				}
				
				// next one
				node.pop();
			}
			if(debug)
			{
				writer.close();
			}

			double YObservedAvg = YObservedSum / N;
			double YPredictedAvg = YPredictedSum / N;

			double PRODUCT = 0;
			double SQUARE_X = 0;
			double SQUARE_Y = 0;

			node.reset();
			while ((testObject = node.peek()) != null)
			{

				double prob = 1;
				
				if(useLeftFlank)
				{
					prob =  getPredictedCountWithLeftFlank(this.leftFlank, 
							this.kmerLength, (Sequence) testObject.getKey(),  this.probabilities);
				}
				else
				{
					prob = getPredictedCount(this.kmerLength, this.noOfValidRead,
							 (Sequence) testObject.getKey(), this.counts, this.probabilities);
				}
							
				long expectedCount = (long) (prob * testTotalCount);
				// update stats
				PRODUCT += (expectedCount - YPredictedAvg)
						* (testObject.getCount() - YObservedAvg);
				SQUARE_X += Math.pow((expectedCount - YPredictedAvg), 2);
				SQUARE_Y += Math.pow((testObject.getCount() - YObservedAvg), 2);
				// next one
				node.pop();
			}

			DebugLog.log("Model Length : " + this.kmerLength);
			DebugLog.log("N : " + N);
			DebugLog.log("YObservedAvg  : " + YObservedAvg);
			DebugLog.log("YPredictedAvg : " + YPredictedAvg);
			DebugLog.log("SQUARE_X: " + SQUARE_X);
			DebugLog.log("SQUARE_Y: " + SQUARE_Y);
			double R = PRODUCT / (Math.sqrt(SQUARE_X) * Math.sqrt(SQUARE_Y));
			DebugLog.log("R: " + R);

			return R;
		}
	}

	public Properties getStats()
	{
		Properties prop = new Properties();
		prop.put("lowestCount", minimalCount);
		prop.put("highestCount", maximalCount);
		prop.put("noOfValidRead", noOfValidRead);
		prop.put("fastqFileProcessed", fastqFileProcessed);
		
		if(fastqFileProcessed)
		{
			prop.put("fastqTotalLineCount", fastqTotalLineCount);
			prop.put("fastqTotalReadCount", fastqTotalReadCount);
			prop.put("fastqBarcodeMatchedReadCount", fastqBarcodeMatchedReadCount);
			prop.put("fastqValidReadCount", fastqValidReadCount);
			
			prop.put("fastqCacheFileName", fastqCacheFileName);
		}
		
		return prop;
	}

	private float[] probabilities;
	private float[] probabilities2;

	/**
	 * Calculates the transition table for K-mers. The calculation bases only on K-mer count table.
	 * Exmple :
	 * 	   P(A|TGG) = P( counts(TGGA) ) / P( counts(TGGx) )
	 * 
	 */
	public void buildMarkovModelEfficient(MarkovModelOption mmOption)
	{
		if (this.externalSorting)
		{
			// To use iterator/visitor pattern to consolidate both models.
			// DebugLog.log(e); throw new RuntimeException("Not supported yet.");
			System.err.println("buildMarkovModel not suppported for 16+ mer");
		} else
		{
			// Arrays.sort(temp);

			int mininalCount = -1;
			this.outputPath = this.outputPath + ".prob";
			mmOption.setModelLogFile(this.outputPath);
			mmOption.setModelOutputFile(this.outputPath+".obj");

			DebugLog.log("Outputting to file:" + this.outputPath);
			PrintWriter writer = null;
			try
			{
				writer = new PrintWriter(this.outputPath);
			} catch (FileNotFoundException e)
			{
				// TODO Auto-generated catch block
				DebugLog.log(e);
			}

			// transition probabilities
			probabilities = new float[this.counts.length];
			// k-1 frequency/probabilities
			probabilities2 = new float[this.counts.length / 4];

			Long prefixMask = ((1L << (2 * this.kmerLength - 2)) - 1) << 2;
			DebugLog.log("prefixMask:" + prefixMask);
			int tmpSum = 0;
			long allCounts = 0;

			for (int i = 0; i < counts.length; i++)
			{
				tmpSum += counts[i];

				if (i % 4 == 3)
				{
					for (int k = 0; k < 4; k++)
					{
						int l = i - 3 + k;
						probabilities[l] = counts[l] * 1F / tmpSum;
						writer.println(new Sequence(l, this.kmerLength)
								.getString()
								+ " c="
								+ counts[l]
								+ " s="
								+ tmpSum + " p=" + (probabilities[l]));
					}
					allCounts += tmpSum;
					probabilities2[i >> 2] = tmpSum;
					tmpSum = 0;
				}

			}

			DebugLog.log("Total " + kmerLength + "-mer length Unique read:"
					+ allCounts);
			for (int i = 0; i < probabilities2.length; i++)
			{
				probabilities2[i] /= allCounts;
			}

			DebugLog.log("Mininal Count:" + mininalCount);

			writer.close();
			DebugLog.log("Closed file:" + this.outputPath);

			// output the probability table
			try
			{
				FileOutputStream fos = new FileOutputStream(mmOption.getModelOutputFile());
				ObjectOutputStream oos = new ObjectOutputStream(fos);
				oos.writeObject(probabilities);
				oos.writeObject(probabilities2);
				oos.flush();
				oos.close();
			} catch (Exception ex)
			{
				DebugLog.log(ex);
			}
		}
	}


	/**
	 * Calculates the transition table for K-mers. The calculation bases only on K-mer count table and (K-1)-mer count table.
	 * Exmple :
	 * 	   P(A|TGG) = P( counts(TGGA) in K-mer table ) / P( counts(TGG) in (K-1)-mer table )
	 * 
	 */
	
	public void buildMarkovModelDivisionMethod(String subKmerCountPath,
			long subKmerTotalCount, MarkovModelOption mmOption)
	{
		if (this.externalSorting)
		{
			// To use iterator/visitor pattern to consolidate both models.
			// DebugLog.log(e); throw new RuntimeException("Not supported yet.");
			System.err.println("buildMarkovModel not suppported for 16+ mer");
		} else
		{
			//loads the k-1 mer count table
			int[] subCounts = null;
			if (subKmerCountPath != null)
			{
				DebugLog.log("Reading subKmerCount file:" + subKmerCountPath);
				ArraysMergerHeap.MinHeapNode node = new ArraysMergerHeap.MinHeapNode(
						subKmerCountPath);
				CountObject obj = null;
				subCounts = new int[1 << (2 * (this.kmerLength - 1))];

				while ((obj = node.peek()) != null)
				{
					Sequence seq = (Sequence) obj.getKey();
					subCounts[(int) seq.getValue()] = obj.getCount();
					node.pop();
				}
			}

			// transition probabilities
			probabilities = null;
			probabilities = new float[this.counts.length];
			String probOutputPath = this.outputPath + ".prob";

			mmOption.setModelLogFile(probOutputPath);
			mmOption.setModelOutputFile(probOutputPath + ".obj2");
			
			DebugLog.log("Outputting to file:" + probOutputPath);
			PrintWriter writer = null;
			try
			{
				writer = new PrintWriter(probOutputPath);
			} catch (FileNotFoundException e)
			{
				DebugLog.log(e);
				throw new RuntimeException(e);
			}
			// division method
			long shortMask = ((1 << (2 * (this.kmerLength - 1))) - 1) << 2;
			for (int i = 0; i < counts.length; i++)
			{
				long total = 0;

				if (subCounts != null) //P(C|AB) = P(ABC)/P(AB) = (count(ABC)/#3-mer) / (count(BC)/#2-mer)
				{
					long valShort = ((shortMask & i) >> 2);
					// DebugLog.log("Long:"+new Sequence(i,
					// this.kmerLength).getString()+" Sub:"+new
					// Sequence(valShort, this.kmerLength-1).getString() );
					this.probabilities[i] = (1F * this.counts[i] / this.noOfValidRead)
							/ (1F * subCounts[(int) valShort] / subKmerTotalCount);
					total = subCounts[(int) valShort];
				} else
				{
					this.probabilities[i] = 1F * this.counts[i]
							/ this.noOfValidRead;
					total = this.noOfValidRead;
				}

				writer.println(new Sequence(i, this.kmerLength).getString()
						+ " c=" + counts[i] + " s= " + total + " p="
						+ (probabilities[i]));
			}
			//

			writer.close();
			DebugLog.log("Closed file:" + this.outputPath);

			// output the probability table
			try
			{
				FileOutputStream fos = new FileOutputStream(mmOption.getModelOutputFile());
				ObjectOutputStream oos = new ObjectOutputStream(fos);
				oos.writeObject(probabilities);
				oos.flush();
				oos.close();
			} catch (Exception ex)
			{
				DebugLog.log(ex);
			}
		}
	}


	private void updateCounts(Sequence str, int num)
	{
		counts[(int) str.getValue()] += num;
	}

}
