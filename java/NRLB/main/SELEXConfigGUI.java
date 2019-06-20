package main;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.*;

class SELEXConfigGUI extends JFrame
{
	/**
	 * 
	 */
	private static final long serialVersionUID = -5937430892049935167L;

	public SELEXConfigGUI()
	{
		addWindowListener(new Terminator());
		setTitle("My Empty Frame");
		setSize(300, 200); // default size is 0,0
		setLocation(10, 200); // default is 0,0 (top left corner)
	}

	public static void main(String[] args)
	{
		JFrame f = new SELEXConfigGUI();
		f.setVisible(true);
	}

	class Terminator extends WindowAdapter
	{
		public void windowClosing(WindowEvent e)
		{
			System.exit(0);
		}
	}
}