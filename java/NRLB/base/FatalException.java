package base;

public class FatalException extends RuntimeException
{
	private static final long serialVersionUID = 244931767831309359L;

	public FatalException()
	{
		super();
	}

	public FatalException(String msg)
	{
		super(msg);
	}
}
