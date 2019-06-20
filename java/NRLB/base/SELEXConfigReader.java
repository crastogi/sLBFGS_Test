package base;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;

import config.AdditionalDataSet;
import config.ExperimentReference;
import config.Round0;
import config.SELEXProcessConfig;
import config.SELEXSequencingConfig;
import config.Sample;
import config.SequencingRunInfo;
import config.Settings;

/**
 * SELEXConfigReader reads and holds the configuration of samples
 */
public class SELEXConfigReader
{
	private SELEXProcessConfig processConfig;
	private HashSet<SeqSampleRound> configSet;
	private HashMap<String,String> seqFileMapping;
	
	public SELEXConfigReader()
	{
		configSet=new HashSet<SeqSampleRound>();
		seqFileMapping = new HashMap<String,String>();
	}
	
	public void clear()
	{
		configSet.clear();
		seqFileMapping.clear();
		processConfig=null;
	}
	
	public static void main(String[] args) throws Exception
	{
		SELEXConfigReader reader= new SELEXConfigReader();
		String xmlPath ="config-Process.xml";
		//reader.readConfig(xmlPath,null);
		//reader.readAdditionalDS("./tmp/stats.txt");
		reader.readConfig(xmlPath, null);
	}
	
	/**
	 * Outputs an object in the 'config' package to a XML file
	 * @param xmlPath
	 * @param obj
	 * @throws Exception
	 */
	public static void writeConfig(String xmlPath,Object obj) throws Exception
	{
        JAXBContext context = JAXBContext.newInstance("config");
        Marshaller marshaller = context.createMarshaller();
        marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
        marshaller.marshal(obj, new FileOutputStream(xmlPath));
	}
	
	public void readAdditionalDS(List<AdditionalDataSet> additionalDS)
	{	
		reloadAdditionalSequenceConfig(additionalDS);
	}
	
	/**
	 * XML file and the data files are assumed to be in the same folder.
	 * If not, dataFolder can be specified.
	 */
	public void readConfig(String xmlPath, String dataFolder) throws Exception
	{
        JAXBContext context = JAXBContext.newInstance("config");
        Unmarshaller marshaller = context.createUnmarshaller();
        
        File f=new File(xmlPath);
    	String p=f.getParent();
    	
    	if(!Util.isEmpty(dataFolder))
    	{
    		p=dataFolder;
    	}
    	
    	Object configObj = marshaller.unmarshal(new FileInputStream(xmlPath));
    	
    	if(configObj instanceof SELEXProcessConfig)
    	{
			processConfig=(SELEXProcessConfig)configObj;
			
			for(String path:processConfig.getSequencingConfigPath())
			{
				if(!Util.isAbsolutePath(path))
				{
					path = f.getParent()+"/"+ path;
				}
				DebugLog.log("Reading config .. "+ path);
				SELEXSequencingConfig config=(SELEXSequencingConfig)marshaller.unmarshal(
			    		new FileInputStream(path));
	    		for(SequencingRunInfo info:config.getSequencingRunInfo())
	    		{
	        		//fix the file names
	        		String path2 = info.getDataFile();
	        		if(!Util.isAbsolutePath(path2))
	        		{
	        			String newPath=p+File.separator+path2;
	        			info.setDataFile(newPath);
	        			DebugLog.log("Replacing path from ["+path2+"]to ["+newPath+"].");
	        			if(!new File(newPath).exists())
	        			{
	        				throw new RuntimeException("File ["+newPath+"] doesn't exist.");
	        			}
	        		}
					addSequencingConfig(info);
	    		}
			//  DebugLog.log("Sequencing Info:"+config.getSequencingRunInfo().getDataFile());
			}
    	}
    	else if (configObj instanceof SELEXSequencingConfig)
    	{
    		SELEXSequencingConfig seqCon = (SELEXSequencingConfig)configObj;
    		
    		for(SequencingRunInfo info:seqCon.getSequencingRunInfo())
    		{
        		//fix the file names
        		String path = info.getDataFile();
        		if(!Util.isAbsolutePath(path))
        		{
        			String newPath=p+File.separator+path;
        			info.setDataFile(newPath);
        			DebugLog.log("Replacing path from ["+path+"]to ["+newPath+"].");
        			if(!new File(newPath).exists())
        			{
        				throw new RuntimeException("File ["+newPath+"] doesn't exist.");
        			}
        		}
        		addSequencingConfig(info);
    		}
    	}
    	else
    	{
    		throw new RuntimeException("Unknowd config object:"+configObj);
    	}
	}
	
	/**
	 * Loads other XML configurations in the [stats.xml]
	 * @param additionDS
	 */
	public void reloadAdditionalSequenceConfig(List<AdditionalDataSet> additionDS)
	{
		for(AdditionalDataSet entry:additionDS)
		{
			String key=entry.getConfigurationPath();
			try
			{
				readConfig(key, null);
			} catch (Exception e)
			{
				DebugLog.log(e);
			}
		}
	}

	/**
	 * Adds SequencingRunInfo objects to the memory.
	 * It replaces the relative paths with the absolute paths.
	 * @param runInfo
	 */
	public void addSequencingConfig(SequencingRunInfo runInfo)
	{
		String seqName = runInfo.getName();
		String dataFileOld= this.seqFileMapping.get(seqName);
		if(!Util.isEmpty(runInfo.getDataFile()))
		{
			try
			{
				String dataFile = new File(runInfo.getDataFile()).getCanonicalPath();
				if(dataFileOld!= null )
				{
					dataFileOld = new File(dataFileOld).getCanonicalPath();
					if(!dataFile.equals(dataFileOld))
					{
						throw new RuntimeException("SeqName["+seqName+"] is already mapped to the fastq file["+dataFileOld+"]");
					}
				}	
				seqFileMapping.put(seqName, dataFile);
			}catch(Exception ex)
			{
				throw new RuntimeException(ex);
			}
		}
		else
		{
			if(dataFileOld==null)
			{
				throw new RuntimeException("SeqName["+seqName+"] is not mapped to a fastq file.");
			}
			runInfo.setDataFile(dataFileOld);
		}
		
		fillInDefaultValues(runInfo);
		
		for(Sample sample : runInfo.getSample())
		{
			fillInDefaultValues(sample);
			SeqSampleRound item = new SeqSampleRound();
			item.setSample(sample);
			item.setSeqRunInfo(runInfo);
			
			Round0 r0=sample.getRound0();
			if(r0!=null)
			{
				item.round0SeqName = r0.getSequencingName();
				item.round0SampleName  = r0.getSampleName();
			}
			
			boolean isDuplicate = configSet.contains(item);
			if(isDuplicate)
			{
				configSet.remove(item);
				DebugLog.log("Duplicate sample found. Overwriting : "+item);
			}
			DebugLog.log("Adding:"+item);
			configSet.add(item);
		}
	}
	
	/**
	 * Returns the round 0 sample.
	 * @param s
	 * @return
	 */
	public ExperimentReference getRound0(ExperimentReference s)
	{
		SeqSampleRound item=new SeqSampleRound();
		item.setExperimentReference(s);

		SeqSampleRound target = null;
		for(SeqSampleRound sample:this.configSet)
		{
			if(sample.equals(item))
			{
				target=new SeqSampleRound();
				target.round = 0;
				target.sampleName = sample.round0SampleName;
				target.seqName = sample.round0SeqName;
			}
		}

		if(target!=null && target.seqName!=null &&  target.sampleName!=null)
		{
			ExperimentReference result=new ExperimentReference();
			result.setSampleRound(0);
			result.setSequencingName(target.seqName);
			result.setSampleName(target.sampleName);
			Sample targetSample = getSample(result);
			if(targetSample==null)
			{
				DebugLog.log("Can't find the Round-0 sample of ["+item+"]");
				return null;
			}
			else
			{
				return result;
			}
		}
		else
		{
			DebugLog.log("Round-0 sample of ["+item+"] is not specified.");
			return null;
		}
	}
	
	/**
	 * Puts empty strings for all non-essential fields.
	 * @param runInfo
	 */
	private void fillInDefaultValues(SequencingRunInfo runInfo)
	{
		if(runInfo.getDescription()==null)	runInfo.setDescription("");
		if(runInfo.getNotes()==null)		runInfo.setNotes("");
		if(runInfo.getResearcherEmail()==null)		runInfo.setResearcherEmail("");
		if(runInfo.getResearcherName()==null)		runInfo.setResearcherName("");
		if(runInfo.getSequencingFacilityEmail()==null)		runInfo.setSequencingFacilityEmail("");
		if(runInfo.getSequencingFacilityName()==null)		runInfo.setSequencingFacilityName("");
		if(runInfo.getSequencingPlatform()==null)			runInfo.setSequencingPlatform("");
	}

	/**
	 * Puts empty strings for all non-essential fields.
	 * @param sample
	 */
	private void fillInDefaultValues(Sample sample)
	{
		if(sample.getProtein()==null)	sample.setProtein("");
		if(sample.getConcentration()==null)	sample.setConcentration("");
		if(sample.getNotes()==null)	sample.setNotes("");
	}
	
	/**
	 * Constructs a SELEXSequencingConfig object and outputs it to a xml file
	 * @param outputPath
	 */
	public void exportSamples(String outputPath)
	{
		SELEXSequencingConfig seqCon = new SELEXSequencingConfig();
		HashSet<String> exported = new HashSet<String>();
		for(SeqSampleRound item:configSet)
		{
			if(exported.contains(item.seqName))
			{
				continue;
			}
			exported.add(item.seqName);
			seqCon.getSequencingRunInfo().add(item.seqRunInfo);
		}
		try
		{
			writeConfig(outputPath,seqCon);
			DebugLog.log("File output:"+outputPath);
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}
	}
	
	public SELEXProcessConfig getProcessConfig()
	{
		return processConfig;
	}

	public Settings getModelSettings()
	{
		if(processConfig==null)
			return null;
		if(processConfig.getProcess()==null)
			return null;
		if(processConfig.getProcess().getModel()==null)
			return null;
		return processConfig.getProcess().getModel().getSettings();
	}
	
	public SequencingRunInfo getSequencingRunInfo(ExperimentReference expRef)
	{
		SeqSampleRound item=new SeqSampleRound();
		item.setExperimentReference(expRef);
		for(SeqSampleRound sample:this.configSet)
		{
			if(sample.equals(item))
			{
				return sample.seqRunInfo;
			}
		}
		
		DebugLog.log("Can't find SequencingRunInfo["+item.seqName+"]");
		return null;
	}
	

	public String getSequencingRunInfoDataFile(String seqName)
	{
		String dataFile = this.seqFileMapping.get(seqName);
		if(dataFile==null)
			DebugLog.log("Can't find SequencingRunInfo["+seqName+"]");
		return dataFile;
	}
	
	
	public Sample getSample(ExperimentReference expRef)
	{
		SeqSampleRound item=new SeqSampleRound();
		item.setExperimentReference(expRef);
		for(SeqSampleRound sample:this.configSet)
		{
			if(sample.equals(item))
			{
				return sample.sample;
			}
		}
		
		DebugLog.log("Can't find Sample["+item+"]");
		return null;
	}

	/**
	 * 
	 * @return The summary information of the loaded samples.
	 */
	public Object[] showSamples()
	{
		LinkedList<String> list1=new LinkedList<String>();
		LinkedList<String> list2=new LinkedList<String>();
		LinkedList<Integer> list3=new LinkedList<Integer>();
		LinkedList<String> list4=new LinkedList<String>();
		LinkedList<String> list5=new LinkedList<String>();

		LinkedList<String> list6_1=new LinkedList<String>();
		LinkedList<String> list6_2=new LinkedList<String>();
		
		LinkedList<String> list6=new LinkedList<String>();
		for(SeqSampleRound info: configSet)
		{
			Sample sample = info.sample;
			list1.add(info.seqName);
			list2.add(sample.getName());
			list3.add(sample.getRound());
			list4.add(sample.getLeftBarcode());
			list5.add(sample.getRightBarcode());
			list6.add(info.seqRunInfo.getDataFile());

			list6_1.add(sample.getLeftFlank());
			list6_2.add(sample.getRightFlank());
		}
		
		String[] list1Array = new String[list3.size()];
		String[] list2Array = new String[list3.size()];
		list1.toArray(list1Array);
		list2.toArray(list2Array);
		int[] list3Array=new int[list3.size()];
		for(int i=0;i<list3.size();i++)
		{
			list3Array[i]=list3.get(i);
		}
		String[] list4Array = new String[list3.size()];
		String[] list5Array = new String[list3.size()];
		String[] list6Array = new String[list3.size()];

		String[] list6Array_1 = new String[list3.size()];
		String[] list6Array_2 = new String[list3.size()];
		list4.toArray(list4Array);
		list5.toArray(list5Array);
		list6.toArray(list6Array);
		list6_1.toArray(list6Array_1);
		list6_2.toArray(list6Array_2);
		return new Object[]{list1Array,list2Array,  list3Array,list4Array,
				list5Array,list6Array_1, list6Array_2, list6Array};
	}

	/**
	 * A helper class to organize SequencingRunInfo and Sample objects. 
	 *
	 */
	private class SeqSampleRound
	{
		String seqName;
		String sampleName;
		Integer round;
		SequencingRunInfo seqRunInfo;
		Sample sample;
		
		String round0SeqName;
		String round0SampleName;
		
		@Override
		public String toString()
		{
			return  seqName + "|" + sampleName + "|" + round;
		}
		
		public void setSeqRunInfo(SequencingRunInfo seqRunInfo)
		{
			this.seqRunInfo = seqRunInfo;
			this.seqName = seqRunInfo.getName();
		}

		public void setSample(Sample sample)
		{
			this.sample = sample;
			this.sampleName = sample.getName();
			this.round = sample.getRound();
		}
		
		public void setExperimentReference(ExperimentReference ref)
		{
			this.seqName = ref.getSequencingName();
			this.sampleName = ref.getSampleName();
			this.round = ref.getSampleRound();
		}
		
		@Override
		public int hashCode()
		{
			final int prime = 31;
			int result = 1;
			result = prime * result + getOuterType().hashCode();
			result = prime * result + ((round == null) ? 0 : round.hashCode());
			result = prime * result
					+ ((sampleName == null) ? 0 : sampleName.hashCode());
			result = prime * result
					+ ((seqName == null) ? 0 : seqName.hashCode());
			return result;
		}
		
		@Override
		public boolean equals(Object obj)
		{
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			SeqSampleRound other = (SeqSampleRound) obj;
			if (!getOuterType().equals(other.getOuterType()))
				return false;
			if (round == null)
			{
				if (other.round != null)
					return false;
			} else if (!round.equals(other.round))
				return false;
			if (sampleName == null)
			{
				if (other.sampleName != null)
					return false;
			} else if (!sampleName.equals(other.sampleName))
				return false;
			if (seqName == null)
			{
				if (other.seqName != null)
					return false;
			} else if (!seqName.equals(other.seqName))
				return false;
			return true;
		}
		
		private SELEXConfigReader getOuterType()
		{
			return SELEXConfigReader.this;
		}
	}
}
