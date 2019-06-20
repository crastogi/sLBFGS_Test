package base;

import java.util.Iterator;

public class CountObjectIterator implements Iterator<CountObject>
{
	private int[] counts;
	private int countIdx;
	private boolean externalSort = false;
	private int kmerLength;
	private ArraysMergerHeap.MinHeapNode node;
	
	public CountObjectIterator(int[] counts, int kmerLength)
	{
		this.counts = counts;
		this.countIdx = 0;
		this.kmerLength=kmerLength;
	}
	
	public CountObjectIterator(String dataFilePath)
	{
		node=new ArraysMergerHeap.MinHeapNode(dataFilePath);
		this.externalSort=true;
	}
	
	@Override
	public boolean hasNext()
	{
		if(!externalSort)
		{
			//goes to an index that is not 0
			while( countIdx<counts.length && counts[countIdx]==0)
			{
				countIdx++;
			}
			return countIdx<counts.length;
		}
		else
		{
			return node.peek()!=null;
		}
	}

	@Override
	public CountObject next()
	{
		if(!externalSort)
		{
			CountObject obj  = null;
			//goes to the next non-0 entry
			while(hasNext())
			{
				//TODO: no need to check this. 
				if(counts[countIdx]>0)
				{
					obj = new CountObject(new Sequence(countIdx,
							this.kmerLength), counts[countIdx]);
					countIdx++;
					break;
				}
				//countIdx++;
			}
			return obj;
		}
		else
		{
			CountObject obj = node.peek();
			node.pop();
			return obj;
		}
	}

	@Override
	public void remove()
	{
	}

}
