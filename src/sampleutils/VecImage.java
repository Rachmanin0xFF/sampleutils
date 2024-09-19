package sampleutils;


import processing.core.*;

/**
 * PVectorImage is 
 */

public class VecImage implements PConstants {
	
	// myParent is a reference to the parent sketch
	Vecf[][] pixels;
	
	public final static String VERSION = "##library.prettyVersion##";
	

	/**
	 * a Constructor, usually called in the setup() method in your sketch to
	 * initialize and start the Library.
	 * 
	 * @example Hello
	 * @param theParent the parent PApplet
	 */
	public VecImage() {
		welcome();
	}
	
	
	private void welcome() {
		System.out.println("##library.name## ##library.prettyVersion## by ##author##");
	}
}

