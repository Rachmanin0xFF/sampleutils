package sampleutils;


import processing.core.*;

/**
 * PVectorImage is 
 */

public class VecImage implements PConstants {
	
	// myParent is a reference to the parent sketch
	Vecf[][] pixels;
	
	public final static String VERSION = "##library.prettyVersion##";
	
	public static enum FilterMode {
		NEAREST,
		BILINEAR,
		BICUBIC,
		BCSPLINE
	}
	public static enum BCSplineMode {
		HERMITE,            // B=0,   C=0
		BSPLINE,            // B=1,   C=0 (smoothest)
		SHARP_BICUBIC,      // B=0,   C=1
		MITCHELL_NETRAVALI, // B=1/3, C=1/3 (recommended)
		CATMULL_ROM,        // B=0,   C=1/2
		PHOTOSHOP_CUBIC     // B=0,   C=3/4
	}
	public static enum EdgeMode {
		CLAMP,
		WRAP,
		BLACK,
		WHITE
	}

	/**
	 * a Constructor, usually called in the setup() method in your sketch to
	 * initialize and start the Library.
	 * 
	 * @example Hello
	 * @param width the width of the image in pixels
	 * @param height the the height of the image in pixels
	 */
	public VecImage(int width, int height, int channels) {
		pixels = new Vecf[width][height];
		for(int x = 0; x < width; x++) for(int y = 0; y < height; y++) {
			pixels[x][y] = new Vecf(channels);
		}
	}
	
	public VecImage(PImage p) {
		pixels = new Vecf[p.width][p.height];
		for(int x = 0; x < p.width; x++) for(int y = 0; y < p.height; y++) {
			int c = p.get(x, y);
			pixels[x][y] = new Vecf(r(c), g(c), b(c), a(c));
		}
	}
	
	int a(int c) {return (c >> 24) & 255; }
	int r(int c) {return (c >> 16) & 255; }
	int g(int c) {return (c >> 8) & 255;}
	int b(int c) {return c & 255; }
	
	private void welcome() {
		System.out.println("##library.name## ##library.prettyVersion## by ##author##");
	}
}

