package sampleutils;

import processing.core.*;

/**
 * PVectorImage is
 */

public class VecImage {
	// myParent is a reference to the parent sketch
	Vecf[][] pixels;

	public final static String VERSION = "##library.prettyVersion##";
	public int width;
	public int height;
	public int channels;

	EdgeModes edgeMode = EdgeModes.CLAMP;
	FilterModes filterMode = FilterModes.BILINEAR;
	BCSplineModes BCSplineMode = BCSplineModes.MITCHELL_NETRAVALI;

	public float BC_SPLINE_B = 0.f;
	public float BC_SPLINE_C = 0.f;

	public static enum FilterModes {
		NEAREST, BILINEAR, BICUBIC, BCSPLINE
	}

	public static enum BCSplineModes {
		HERMITE, // B=0, C=0
		BSPLINE, // B=1, C=0 (smoothest)
		SHARP_BICUBIC, // B=0, C=1
		MITCHELL_NETRAVALI, // B=1/3, C=1/3 (recommended)
		CATMULL_ROM, // B=0, C=1/2
		PHOTOSHOP_CUBIC, // B=0, C=3/4
		CUSTOM
	}

	public static enum EdgeModes {
		CLAMP, WRAP, BLACK, WHITE
	}

	/**
	 * a Constructor, usually called in the setup() method in your sketch to
	 * initialize and start the Library.
	 * 
	 * @example Hello
	 * @param width  the width of the image in pixels
	 * @param height the the height of the image in pixels
	 */
	public VecImage(int width, int height, int channels) {
		pixels = new Vecf[width][height];
		for (int x = 0; x < width; x++)
			for (int y = 0; y < height; y++) {
				pixels[x][y] = new Vecf(channels);
			}
		this.width = width;
		this.height = height;
		this.channels = channels;
	}

	public VecImage(PImage p) {
		pixels = new Vecf[p.width][p.height];
		for (int x = 0; x < p.width; x++)
			for (int y = 0; y < p.height; y++) {
				int c = p.get(x, y);
				pixels[x][y] = new Vecf(r(c) / 255.f, g(c) / 255.f, b(c) / 255.f, a(c) / 255.f);
			}
		this.width = p.width;
		this.height = p.height;
		this.channels = 4;
	}

	// Makes a shallow copy
	public VecImage(VecImage p) {
		pixels = new Vecf[p.width][p.height];
		this.width = p.width;
		this.height = p.height;
		this.channels = p.channels;
		for (int x = 0; x < p.width; x++)
			for (int y = 0; y < p.height; y++) {
				pixels[x][y] = new Vecf(p.pixels[x][y]);
			}
	}

	public final int color(int r, int g, int b, int alpha) {
		return ((alpha < 0 ? 0 : (alpha > 255 ? 255 : alpha)) << 24) | ((r < 0 ? 0 : (r > 255 ? 255 : r)) << 16)
				| ((g < 0 ? 0 : (g > 255 ? 255 : g)) << 8) | ((b < 0 ? 0 : (b > 255 ? 255 : b)));
	}

	// for use with Processing's PImage
	public int[] toIntegerPixels() {
		int[] output = new int[width * height];
		for (int x = 0; x < width; x++)
			for (int y = 0; y < height; y++) {
				if (channels > 3)
					output[x + y * width] = color((int) (pixels[x][y].components[0] * 255.f),
							(int) (pixels[x][y].components[1] * 255.f), (int) (pixels[x][y].components[2] * 255.f),
							(int) (pixels[x][y].components[3] * 255.f));
				else if (channels == 3)
					output[x + y * width] = color((int) (pixels[x][y].components[0] * 255.f),
							(int) (pixels[x][y].components[1] * 255.f), (int) (pixels[x][y].components[2] * 255.f),
							255);
			}
		return output;
	}

	public Vecf getPixel(int x, int y) {
		switch (edgeMode) {
		case CLAMP:
			x = PixelMath.clamp_inclusive(x, 0, width - 1);
			y = PixelMath.clamp_inclusive(y, 0, height - 1);
			break;
		case WRAP:
			x = PixelMath.mod(x, width);
			y = PixelMath.mod(x, height);
			break; // TODO
		case WHITE:
			if (x >= width || y >= height || x < 0 || y < 0)
				return Vecf.onesLike(pixels[x][y]);
			break;
		case BLACK:
			if (x >= width || y >= height || x < 0 || y < 0)
				return Vecf.zeroesLike(pixels[x][y]);
			break;
		}
		return pixels[x][y];
	}

	public Vecf sample(float x, float y) {
		switch (filterMode) {
		case NEAREST:
			return getPixel((int) x, (int) y);
		case BILINEAR:
			float dx = x % 1;
			float dy = y % 1;

			return Vecf.add(
					Vecf.add(Vecf.mult(getPixel((int) Math.floor(x), (int) Math.floor(y)), (1 - dx) * (1 - dy)),
							Vecf.mult(getPixel((int) Math.floor(x + 1), (int) Math.floor(y)), dx * (1 - dy))),
					Vecf.add(Vecf.mult(getPixel((int) Math.floor(x), (int) Math.floor(y + 1)), (1 - dx) * dy),
							Vecf.mult(getPixel((int) Math.floor(x + 1), (int) Math.floor(y + 1)), dx * dy)));
		case BICUBIC:
			return sampleCubic(x, y);
		case BCSPLINE:
			switch (BCSplineMode) {
			case HERMITE:
				return sampleBC(x, y, 0.f, 0.f);
			case BSPLINE:
				return sampleBC(x, y, 1.f, 0.f);
			case CATMULL_ROM:
				return sampleBC(x, y, 0.f, 1.f);
			case MITCHELL_NETRAVALI:
				return sampleBC(x, y, 1.f / 3.f, 1.f / 3.f);
			case PHOTOSHOP_CUBIC:
				return sampleBC(x, y, 0.f, 0.5f);
			case SHARP_BICUBIC:
				return sampleBC(x, y, 0.f, 0.75f);
			case CUSTOM:
				return sampleBC(x, y, BC_SPLINE_B, BC_SPLINE_C);

			}
		default:
			throw new IllegalStateException("Cannot sample; filterMode or BCSplineMode not initialized!");
		}
	}

	Vecf sampleBC(float x, float y, float B, float C) {
		Vecf p0, p1, p2, p3;
		float d = x - (int) Math.floor(x);

		if (d == 0.0) {
			p0 = sampleBCY((int) Math.floor(x) - 1, y, B, C);
			p1 = sampleBCY((int) Math.floor(x), y, B, C);
			p2 = sampleBCY((int) Math.floor(x) + 1, y, B, C);
			p3 = sampleBCY((int) Math.floor(x) + 2, y, B, C);
		} else {
			p0 = sampleBCY((int) Math.floor(x) - 1, y, B, C);
			p1 = sampleBCY((int) Math.floor(x), y, B, C);
			p2 = sampleBCY((int) Math.ceil(x), y, B, C);
			p3 = sampleBCY((int) Math.ceil(x) + 1, y, B, C);
		}
		return Vecf.add(
				Vecf.add(
						Vecf.mult(
								Vecf.add(
										Vecf.add(Vecf.mult(p0, -B / 6.f - C), Vecf.mult(p1,
												-1.5f * B - C + 2.f)),
										Vecf.add(Vecf.mult(p2, 1.5f * B + C - 2.f), Vecf.mult(p3, B / 6.f + C))),
								d * d * d),
						Vecf.mult(
								Vecf.add(Vecf.add(Vecf.mult(p0, 0.5f * B + 2.f * C), Vecf.mult(p1, 2.f * B + C - 3.f)),
										Vecf.add(Vecf.mult(p2, -2.5f * B - 2.f * C + 3.f), Vecf.mult(p3, -C))),
								d * d)),
				Vecf.add(Vecf.mult(Vecf.add(Vecf.mult(p0, -0.5f * B - C), Vecf.mult(p2, 0.5f * B + C)), d), Vecf
						.add(Vecf.mult(p0, B / 6.f), Vecf.add(Vecf.mult(p1, -B / 3.f + 1.f), Vecf.mult(p2, B / 6.f)))));
	}

	Vecf sampleBCY(int x, float y, float B, float C) {
		Vecf p0, p1, p2, p3;
		float d = y - (int) Math.floor(y);

		if (d == 0.0) {
			p0 = getPixel(x, (int) Math.floor(y) - 1);
			p1 = getPixel(x, (int) Math.floor(y));
			p2 = getPixel(x, (int) Math.floor(y) + 1);
			p3 = getPixel(x, (int) Math.floor(y) + 2);
		} else {
			p0 = getPixel(x, (int) Math.floor(y) - 1);
			p1 = getPixel(x, (int) Math.floor(y));
			p2 = getPixel(x, (int) Math.ceil(y));
			p3 = getPixel(x, (int) Math.ceil(y) + 1);
		}
		return Vecf.add(
				Vecf.add(
						Vecf.mult(
								Vecf.add(
										Vecf.add(Vecf.mult(p0, -B / 6.f - C), Vecf.mult(p1,
												-1.5f * B - C + 2.f)),
										Vecf.add(Vecf.mult(p2, 1.5f * B + C - 2.f), Vecf.mult(p3, B / 6.f + C))),
								d * d * d),
						Vecf.mult(
								Vecf.add(Vecf.add(Vecf.mult(p0, 0.5f * B + 2.f * C), Vecf.mult(p1, 2.f * B + C - 3.f)),
										Vecf.add(Vecf.mult(p2, -2.5f * B - 2.f * C + 3.f), Vecf.mult(p3, -C))),
								d * d)),
				Vecf.add(Vecf.mult(Vecf.add(Vecf.mult(p0, -0.5f * B - C), Vecf.mult(p2, 0.5f * B + C)), d), Vecf
						.add(Vecf.mult(p0, B / 6.f), Vecf.add(Vecf.mult(p1, -B / 3.f + 1.f), Vecf.mult(p2, B / 6.f)))));
	}

	// these next two functions are helpers for cubic interpolation during image
	// upscaling
	Vecf sampleCubic(float x, float y) {
		Vecf _p0 = sampleCubicY((int) Math.floor(x) - 1, y);
		Vecf _p1 = sampleCubicY((int) Math.floor(x), y);
		Vecf _p2 = sampleCubicY((int) Math.ceil(x), y);
		Vecf _p3 = sampleCubicY((int) Math.ceil(x) + 1, y);

		float _x = x - (int) Math.floor(x);

		// credit to https://www.paulinternet.nl/?page=bicubic for crunching the algebra
		// here

		return Vecf.add(_p1, Vecf.mult(
				Vecf.add(Vecf.sub(_p2, _p0),
						Vecf.mult((Vecf.add(
								Vecf.add(Vecf.sub(Vecf.mult(_p0, 2.f), Vecf.mult(_p1, 5.f)),
										Vecf.sub(Vecf.mult(_p2, 4.f), _p3)),
								Vecf.mult(Vecf.add(Vecf.mult(Vecf.sub(_p1, _p2), 3.f), Vecf.sub(_p3, _p0)), _x))), _x)),
				0.5f * _x));
	}

	Vecf sampleCubicY(int x, float y) {
		Vecf _p0 = getPixel(x, (int) Math.floor(y) - 1);
		Vecf _p1 = getPixel(x, (int) Math.floor(y));
		Vecf _p2 = getPixel(x, (int) Math.ceil(y));
		Vecf _p3 = getPixel(x, (int) Math.ceil(y) + 1);

		float _x = y - (int) Math.floor(y);

		// credit to https://www.paulinternet.nl/?page=bicubic for crunching the algebra
		// here

		return Vecf.add(_p1, Vecf.mult(
				Vecf.add(Vecf.sub(_p2, _p0),
						Vecf.mult((Vecf.add(
								Vecf.add(Vecf.sub(Vecf.mult(_p0, 2.f), Vecf.mult(_p1, 5.f)),
										Vecf.sub(Vecf.mult(_p2, 4.f), _p3)),
								Vecf.mult(Vecf.add(Vecf.mult(Vecf.sub(_p1, _p2), 3.f), Vecf.sub(_p3, _p0)), _x))), _x)),
				0.5f * _x));
	}

	int a(int c) {
		return (c >> 24) & 255;
	}

	int r(int c) {
		return (c >> 16) & 255;
	}

	int g(int c) {
		return (c >> 8) & 255;
	}

	int b(int c) {
		return c & 255;
	}

	private void welcome() {
		System.out.println("##library.name## ##library.prettyVersion## by ##author##");
	}
}
