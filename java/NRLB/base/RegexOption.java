package base;

import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;

public class RegexOption
{
	//before and during k-mer counting
    private String variableRegionIncludeRegex=null;
    private String variableRegionExcludeRegex=null;
    private String variableRegionGroupRegex=null;
    private String kmerIncludeRegex=null;
    private String kmerExcludeRegex=null;
    private Set<String> kmerIncludeOnly=null;
    //after k-mer counting options
    private String viewIncludeRegex=null;
    private String viewExcludeRegex=null;
    private Set<String> viewIncludeOnly=null;
    //
    private String variableRegionRegexFormattedString;
    private String kmerRegexFomattedString;
    //
    private Pattern variableRegionIncludePattern = null;
    private Pattern variableRegionExcludePattern = null;
    private Pattern variableRegionGroupPattern = null;
    private Pattern kmerIncludePattern = null;
    private Pattern kmerExcludePattern = null;
    private Pattern kmerIncludeOnlyPattern = null;
    private Pattern viewIncludePattern =null;
    private Pattern viewExcludePattern =null;
    private Pattern viewIncludeOnlyPattern = null;
    
    public static boolean isInitalized(RegexOption obj)
    {
    	return obj !=null && 
    		  (obj.variableRegionIncludeRegex !=null ||
			   obj.variableRegionExcludeRegex !=null ||
			   obj.variableRegionGroupRegex !=null ||
			   obj.kmerIncludeRegex !=null ||
			   obj.kmerExcludeRegex !=null ||
			   obj.kmerIncludeOnly !=null);
    }
    
    public String getVariableRegionRegexFormattedString()
    {
    	if(variableRegionRegexFormattedString==null)
    	{
    		variableRegionRegexFormattedString =
    						"variableRegionIncludeRegex:" + variableRegionIncludeRegex +"," +
    						"variableRegionExcludeRegex:" + variableRegionExcludeRegex +"," +
    						"variableRegionGroupRegex:" + variableRegionGroupRegex ;
    	}
    	return variableRegionRegexFormattedString;
    }
    
	public String getKmerRegexFomattedString()
	{
		if(kmerRegexFomattedString==null)
		{
			kmerRegexFomattedString = getVariableRegionRegexFormattedString() + ",";
			kmerRegexFomattedString += 
							"kmerIncludeRegex:" + kmerIncludeRegex +",";
			kmerRegexFomattedString += 
							"kmerExcludeRegex:" + kmerExcludeRegex +",";
			kmerRegexFomattedString += 
							"kmerIncludeOnly:" + getKmerIncludeOnlyString();
		}
		return kmerRegexFomattedString;
	}

	public String getVariableRegionIncludeRegex()
	{
		return variableRegionIncludeRegex;
	}
	public void setVariableRegionIncludeRegex(String variableRegionIncludeRegex)
	{
		this.variableRegionIncludeRegex = variableRegionIncludeRegex;
	}
	public String getVariableRegionExcludeRegex()
	{
		return variableRegionExcludeRegex;
	}
	public void setVariableRegionExcludeRegex(String variableRegionExcludeRegex)
	{
		this.variableRegionExcludeRegex = variableRegionExcludeRegex;
	}
	public String getKmerIncludeRegex()
	{
		return kmerIncludeRegex;
	}
	public void setKmerIncludeRegex(String kmerIncludeRegex)
	{
		this.kmerIncludeRegex = kmerIncludeRegex;
	}
	public String getKmerExcludeRegex()
	{
		return kmerExcludeRegex;
	}
	public void setKmerExcludeRegex(String kmerExcludeRegex)
	{
		this.kmerExcludeRegex = kmerExcludeRegex;
	}
	public Set<String> getKmerIncludeOnly()
	{
		if(kmerIncludeOnly==null)
		{
			kmerIncludeOnly = new TreeSet<String>();
		}
		return kmerIncludeOnly;
	}
	public void setKmerIncludeOnly(Set<String> kmerIncludeOnly)
	{
		this.kmerIncludeOnly = kmerIncludeOnly;
	}
	public String getViewIncludeRegex()
	{
		return viewIncludeRegex;
	}
	public void setViewIncludeRegex(String viewIncludeRegex)
	{
		this.viewIncludeRegex = viewIncludeRegex;
	}
	public String getViewExcludeRegex()
	{
		return viewExcludeRegex;
	}
	public void setViewExcludeRegex(String viewExcludeRegex)
	{
		this.viewExcludeRegex = viewExcludeRegex;
	}
	public Set<String> getViewIncludeOnly()
	{
		if(viewIncludeOnly==null)
		{
			viewIncludeOnly = new TreeSet<String>();
		}
		return viewIncludeOnly;
	}
	public void setViewIncludeOnly(Set<String> viewIncludeOnly)
	{
		this.viewIncludeOnly = viewIncludeOnly;
	}
	
	public Pattern getViewIncludeOnlyPattern()
	{
		if(viewIncludeOnlyPattern==null && viewIncludeOnly!=null && viewIncludeOnly.size() > 0)
		{
			String pattern = Util.joinStrings(viewIncludeOnly, "|");
			
			DebugLog.log("Compiling filter : "+ pattern);
			viewIncludeOnlyPattern = Pattern.compile(pattern);
		}
		return viewIncludeOnlyPattern;
	}
	
	public Pattern getKmerIncludeOnlyPattern()
	{
		if(kmerIncludeOnlyPattern==null && kmerIncludeOnly!=null && kmerIncludeOnly.size() > 0)
		{
			String pattern = Util.joinStrings(kmerIncludeOnly, "|");
			
			DebugLog.log("Compiling filter : "+ pattern);
			kmerIncludeOnlyPattern = Pattern.compile(pattern);
		}
		return kmerIncludeOnlyPattern;
	}
	
	public String getKmerIncludeOnlyString()
	{
		if(kmerIncludeOnly!=null && kmerIncludeOnly.size() > 0)
		{
			return Util.joinStrings(kmerIncludeOnly, "|");
		}
		return null;
	}
	
	public String getVariableRegionGroupRegex()
	{
		return variableRegionGroupRegex;
	}

	public void setVariableRegionGroupRegex(String variableRegionGroupRegex)
	{
		this.variableRegionGroupRegex = variableRegionGroupRegex;
	}

	public Pattern getVariableRegionGroupPattern()
	{
		if(variableRegionGroupPattern==null && variableRegionGroupRegex!=null)
		{
			DebugLog.log("Compiling filter : "+ variableRegionGroupRegex);
			variableRegionGroupPattern = Pattern.compile(variableRegionGroupRegex);
		}
		return variableRegionGroupPattern;
	}
	
	public Pattern getVariableRegionIncludePattern()
	{
		if(variableRegionIncludePattern==null && variableRegionIncludeRegex!=null)
		{
			DebugLog.log("Compiling filter : "+ variableRegionIncludeRegex);
			variableRegionIncludePattern = Pattern.compile(variableRegionIncludeRegex);
		}
		return variableRegionIncludePattern;
	}

	public Pattern getVariableRegionExcludePattern()
	{
		if(variableRegionExcludePattern==null && variableRegionExcludeRegex!=null)
		{
			DebugLog.log("Compiling filter : "+ variableRegionExcludeRegex);
			variableRegionExcludePattern = Pattern.compile(variableRegionExcludeRegex);
		}
		return variableRegionExcludePattern;
	}

	public Pattern getKmerIncludePattern()
	{
		if(kmerIncludePattern==null && kmerIncludeRegex!=null)
		{
			DebugLog.log("Compiling filter : "+ kmerIncludeRegex);
			kmerIncludePattern = Pattern.compile(kmerIncludeRegex);
		}
		return kmerIncludePattern;
	}

	public Pattern getKmerExcludePattern()
	{
		if(kmerExcludePattern==null && kmerExcludeRegex!=null)
		{
			DebugLog.log("Compiling filter : "+ kmerExcludeRegex);
			kmerExcludePattern = Pattern.compile(kmerExcludeRegex);
		}
		return kmerExcludePattern;
	}

	public Pattern getViewIncludePattern()
	{
		if(viewIncludePattern==null && viewIncludeRegex!=null)
		{
			DebugLog.log("Compiling filter : "+ viewIncludeRegex);
			viewIncludePattern = Pattern.compile(viewIncludeRegex);
		}
		return viewIncludePattern;
	}

	public Pattern getViewExcludePattern()
	{
		if(viewExcludePattern==null && viewExcludeRegex!=null)
		{
			DebugLog.log("Compiling filter : "+ viewExcludeRegex);
			viewExcludePattern = Pattern.compile(viewExcludeRegex);
		}
		return viewExcludePattern;
	}

	@Override
	public String toString()
	{
		return "RegexOption [variableRegionIncludeRegex="
				+ variableRegionIncludeRegex + ", variableRegionExcludeRegex="
				+ variableRegionExcludeRegex + ", variableRegionGroupRegex="
				+ variableRegionGroupRegex + ", kmerIncludeRegex="
				+ kmerIncludeRegex + ", kmerExcludeRegex=" + kmerExcludeRegex
				+ ", kmerIncludeOnly=" + kmerIncludeOnly
				+ ", viewIncludeRegex=" + viewIncludeRegex
				+ ", viewExcludeRegex=" + viewExcludeRegex
				+ ", viewIncludeOnly=" + viewIncludeOnly + "]";
	}
	
	
	
}
