package sampleutils;

import processing.core.*;

public class SampleUtilsIO implements PConstants {
	PApplet myParent;
	
	public SampleUtilsIO(PApplet parentApplet) {
		myParent = parentApplet;
	}
	public VecImage loadImage(String path) {
		return toVecImage(myParent.loadImage(path));
	}
	
	public PImage toPImage(VecImage v) {
		PImage p = myParent.createImage(v.width, v.height, ARGB);
		p.loadPixels();
		p.pixels = v.toIntegerPixels();
		p.updatePixels();
		return p;
	}
	public VecImage toVecImage(PImage p) {
		VecImage v = new VecImage(p.width, p.height, 4);
		for (int x = 0; x < p.width; x++)
			for (int y = 0; y < p.height; y++) {
				int c = p.get(x, y);
				v.pixels[x][y] = new Vecf(r(c), g(c), b(c), a(c));
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
