package base;

public class MarkovModelOption
{
	public static final String DIVISION = "DIVISION";
	public static final String TRANSITION = "TRANSITION";

	public static final String WITH_LEFT_FLANK = "WITH_LEFT_FLANK";
	
	private Boolean buildMarkovModel = false;
	private String modelOutputFile;
	private String modelLogFile;
	private String method = DIVISION; //
	
	public Boolean getBuildMarkovModel()
	{
		return buildMarkovModel;
	}
	public void setBuildMarkovModel(Boolean buildMarkovModel)
	{
		this.buildMarkovModel = buildMarkovModel;
	}
	public String getModelOutputFile()
	{
		return modelOutputFile;
	}
	public void setModelOutputFile(String modelOutputFile)
	{
		this.modelOutputFile = modelOutputFile;
	}
	public String getModelLogFile()
	{
		return modelLogFile;
	}
	public void setModelLogFile(String modelLogFile)
	{
		this.modelLogFile = modelLogFile;
	}
	public String getMethod()
	{
		return method;
	}
	public void setMethod(String method)
	{
		this.method = method;
	}
}
