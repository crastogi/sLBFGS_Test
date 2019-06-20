package utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

import base.Sequence;

public class FA {
	public int nBases;
	public String name;
	private ArrayList<String> data;
	
	public FA(String file) {
		nBases = 0;
		String currLine;
		data = new ArrayList<String>();
		BufferedReader br;
		
		try {
			br		 = new BufferedReader(new FileReader(file));
			currLine = br.readLine();
			name	 = currLine.substring(1);
			while ((currLine = br.readLine())!=null) {
				data.add(currLine);
			}
			br.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		nBases = 50*(data.size()-1) + data.get(data.size()-1).length();
	}
	
	/**
	 * 
	 * @param start
	 * @param end
	 * @param status:	-1 indicates the region contains an N
	 * 					0  indicates all good
	 * 					+1 indicates repeats exist within the region
	 * @return
	 */
	
	public long getLong(int start, int end, int[] status) {
		int startBlock	= start/50;
		int endBlock	= end/50;
		String region, modded;
		
		if ((end-start)>30) {
			status[0] = -2;
			return -2;
		}
		
		if (startBlock!=endBlock) {
			region = data.get(startBlock).substring((start % 50));
			for (int i=startBlock+1; i<endBlock; i++) {
				region = region.concat(data.get(startBlock+i));
			}
			region =  region.concat(data.get(endBlock).substring(0, (end % 50)+1));
		} else {
			region =  data.get(startBlock).substring((start % 50), (end % 50)+1);
		}
		modded = region.toUpperCase();
		if (!region.equals(modded)) {
			status[0] = 1;
		} else {
			status[0] = 0;
		}
		if (modded.indexOf("N") >= 0) {
			status[0] = -1;
			return -1;
		}
		return (new Sequence(modded, 0, end-start+1)).getValue();
	}

	public String getRegion(int start, int end) {
		int startBlock	= start/50;
		int endBlock	= end/50;
		String output;
		
		if (startBlock!=endBlock) {
			output = data.get(startBlock).substring((start % 50));
			for (int i=startBlock+1; i<endBlock; i++) {
				output = output.concat(data.get(startBlock+i));
			}
			return output.concat(data.get(endBlock).substring(0, (end % 50)+1));
		} else {
			return data.get(startBlock).substring((start % 50), (end % 50)+1);
		}
	}
}
