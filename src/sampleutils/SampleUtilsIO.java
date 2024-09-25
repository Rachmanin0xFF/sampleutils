package sampleutils;

import processing.core.*;

/**
 * Instancing necessary for using VecImage with native Processing methods.
 * <p>
 * Creating and displaying PImages requires having the PApplet object. This class exists to decouple that
 * requirement from VecImage.
 */
public class SampleUtilsIO implements PConstants {
	PApplet myParent;
	
	/**
	 * Initialize this in setup() by passing <code>this</code> into this constructor.
	 * @param parentApplet the PApplet used by the parent Processing application.
	 */
	public SampleUtilsIO(PApplet parentApplet) {
		myParent = parentApplet;
	}
	/**
	 * A wrapper for Processing's <code>loadImage</code> method.
	 * @param path the path of the input image.
	 * @param includeAlpha whether or not to load the alpha channel
	 * @return the (4-component) VecImage loaded from the path.
	 */
	public VecImage loadImage(String path, boolean includeAlpha) {
		return toVecImage(myParent.loadImage(path), includeAlpha);
	}
	
	/**
	 * Converts a VecImage to a PImage
	 * @param v the VecImage
	 * @return a PImage with the same data
	 */
	public PImage toPImage(VecImage v) {
		PImage p = myParent.createImage(v.width, v.height, ARGB);
		p.loadPixels();
		p.pixels = v.toIntegerPixels();
		p.updatePixels();
		return p;
	}
	
	/**
	 * Converts a PImage to a VecImage
	 * @param p the PImage
	 * @param includeAlpha whether or not to pass along the alpha channel
	 * @return a (4-component) VecImage with the same data
	 */
	public VecImage toVecImage(PImage p, boolean includeAlpha) {
		VecImage v = new VecImage(p.width, p.height, includeAlpha ? 4 : 3);
		for (int x = 0; x < p.width; x++)
			for (int y = 0; y < p.height; y++) {
				int c = p.get(x, y);
				if(includeAlpha) v.pixels[x][y] = new Vecf(r(c), g(c), b(c), a(c));
				else v.pixels[x][y] = new Vecf(r(c), g(c), b(c));
			}
		return v;
	}
	
	protected int a(int c) {
		return (c >> 24) & 255;
	}

	protected int r(int c) {
		return (c >> 16) & 255;
	}

	protected int g(int c) {
		return (c >> 8) & 255;
	}

	protected int b(int c) {
		return c & 255;
	}
}
