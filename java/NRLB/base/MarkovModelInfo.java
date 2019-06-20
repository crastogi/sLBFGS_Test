package base;

import config.ExperimentReference;

public class MarkovModelInfo
{
	private Integer markovLength;
	private Long markovLengthTotalCount;
	private Double markovR2;
	private String markovCountsPath;
	private String markovObjPath;
	private ExperimentReference sample;
	private String markovModelMethod;
	private ExperimentReference crossValidationSample;
	private String filter;
	
	
	@Override
	public String toString()
	{
		return "MarkovModelInfo [markovLength=" + markovLength
				+ ", markovLengthTotalCount=" + markovLengthTotalCount
				+ ", markovR2=" + markovR2 + ", markovCountsPath="
				+ markovCountsPath + ", markovObjPath=" + markovObjPath
				+ ", sample=" + sample + ", markovModelMethod="
				+ markovModelMethod + ", crossValidationSample="
				+ crossValidationSample + ", filter=" + filter + "]";
	}

	public int getMarkovLength()
	{
		return markovLength;
	}

	public Long getMarkovLengthTotalCount()
	{
		return markovLengthTotalCount;
	}

	public double getMarkovR2()
	{
		if(markovR2==null)
			return 0;
		return markovR2;
	}

	public String getMarkovCountsPath()
	{
		return markovCountsPath;
	}

	public String getMarkovObjPath()
	{
		return markovObjPath;
	}

	public void setMarkovLength(Integer markovLength)
	{
		this.markovLength = markovLength;
	}

	public void setMarkovLengthTotalCount(Long markovLengthTotalCount)
	{
		this.markovLengthTotalCount = markovLengthTotalCount;
	}

	public void setMarkovR2(Double markovR2)
	{
		this.markovR2 = markovR2;
	}

	public void setMarkovCountsPath(String markovCountsPath)
	{
		this.markovCountsPath = markovCountsPath;
	}

	public void setMarkovObjPath(String markovObjPath)
	{
		this.markovObjPath = markovObjPath;
	}

	public ExperimentReference getSample()
	{
		return sample;
	}

	public void setSample(ExperimentReference sample)
	{
		this.sample = sample;
	}

	public ExperimentReference getCrossValidationSample()
	{
		return crossValidationSample;
	}

	public void setCrossValidationSample(ExperimentReference crossValidationSample)
	{
		this.crossValidationSample = crossValidationSample;
	}

	public String getMarkovModelMethod()
	{
		return markovModelMethod;
	}

	public void setMarkovModelMethod(String markovModelMethod)
	{
		if(markovModelMethod==null)
			return;
		this.markovModelMethod = markovModelMethod;
	}

	public String getFilter()
	{
		return filter;
	}

	public void setFilter(String filter)
	{
		this.filter = filter;
	}
	
	
}