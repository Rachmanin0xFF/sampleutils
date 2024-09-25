package sampleutils;

import processing.core.*;

/**
 * An image class supporting an arbitrary number of floating-point precision channels.
 * <p>
 * Images also carry information about their sampling/edge filtering preferences in enums and numbers.
 */
public class VecImage {
	// myParent is a reference to the parent sketch
	public Vecf[][] pixels;

	public final static String VERSION = "##library.prettyVersion##";
	public int width;
	public int height;
	public int channels;
	
	public EdgeModes edgeMode = EdgeModes.CLAMP;
	public FilterModes filterMode = FilterModes.BILINEAR;
	public BCSplineModes BCSplineMode = BCSplineModes.MITCHELL_NETRAVALI;

	public float BC_SPLINE_B = 0.f;
	public float BC_SPLINE_C = 0.f;
	public int LANCZOS_RADIUS = 3;
	
	/**
	 * Image filter modes (for both up and down-sampling).
	 */
	public static enum FilterModes {
		NEAREST, BILINEAR, BICUBIC, LANCZOS, BCSPLINE
	}

	/**
	 * Mitchell–Netravali filter (BC-spline) presets.
	 */
	public static enum BCSplineModes {
		HERMITE, // B=0, C=0
		BSPLINE, // B=1, C=0 (smoothest)
		SHARP_BICUBIC, // B=0, C=1
		MITCHELL_NETRAVALI, // B=1/3, C=1/3 (recommended)
		CATMULL_ROM, // B=0, C=1/2
		PHOTOSHOP_BICUBIC, // B=0, C=3/4
		CUSTOM
	}
	
	/**
	 * Texture wrap modes for sampling outside of an image's bounds.
	 */
	public static enum EdgeModes {
		CLAMP, WRAP, BLACK, WHITE
	}

	/**
	 * @param width  the width of the image in pixels
	 * @param height the the height of the image in pixels
	 * @param channels the number of color channel the image has (usually 3 or 4)
	 */
	public VecImage(int width, int height, int channels) {
		if(width < 1 || height < 1) {
			throw new IllegalArgumentException("Image dimensions must be at least (1px x 1px)!");
		}
		pixels = new Vecf[width][height];
		for (int x = 0; x < width; x++)
			for (int y = 0; y < height; y++) {
				pixels[x][y] = new Vecf(channels);
			}
		this.width = width;
		this.height = height;
		this.channels = channels;
	}

	/**
	 * Initializes a shallow copy of another image.
	 * @param p the image to copy data from
	 */
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
	
	public VecImage(float[][] data) {
		this.width = data.length;
		this.height = data[0].length;
		pixels = new Vecf[width][height];
		this.channels = 1;
		for(int x = 0; x < width; x++) {
			for(int y = 0; y < height; y++) {
				pixels[x][y] = new Vecf(data[x][y]);
			}
		}
	}
	
	private final int color(int r, int g, int b, int alpha) {
		return ((alpha < 0 ? 0 : (alpha > 255 ? 255 : alpha)) << 24) | ((r < 0 ? 0 : (r > 255 ? 255 : r)) << 16)
				| ((g < 0 ? 0 : (g > 255 ? 255 : g)) << 8) | ((b < 0 ? 0 : (b > 255 ? 255 : b)));
	}

	/**
	 * Generates an array compatible with Processing's PImage.pixels[].
	 * @return A flattened, width*height-length array of integers representing 24-bit ARGB colors.
	 */
	public int[] toIntegerPixels() {
		int[] output = new int[width * height];
		for (int x = 0; x < width; x++)
			for (int y = 0; y < height; y++) {
				// check if we have enough information for alpha. if not, make fully opaque
				if (channels > 3)
					output[x + y * width] = color((int) (pixels[x][y].components[0]),
							(int) (pixels[x][y].components[1]), (int) (pixels[x][y].components[2]),
							(int) (pixels[x][y].components[3]));
				else if (channels == 3)
					output[x + y * width] = color((int) (pixels[x][y].components[0]),
							(int) (pixels[x][y].components[1]), (int) (pixels[x][y].components[2]),
							255);
				else if(channels == 1)
					output[x + y * width] = color((int) (pixels[x][y].components[0]),
							(int) (pixels[x][y].components[0]), (int) (pixels[x][y].components[0]),
							255);
			}
		return output;
	}
	
	/**
	 * Returns the pixel at the integer coordinates (x, y).
	 * <p>
	 * If x or y are outside the image bounds, return value is based on the image's <code>edgeMode</code>.
	 * @param x
	 * @param y
	 * @return the Vecf associated with the input coordinates
	 */
	public Vecf getPixel(int x, int y) {
		switch (edgeMode) {
		case CLAMP:
			x = PixelMath.clamp_inclusive(x, 0, width - 1);
			y = PixelMath.clamp_inclusive(y, 0, height - 1);
			break;
		case WRAP:
			x = PixelMath.mod(x, width);
			y = PixelMath.mod(y, height);
			break; // TODO
		case WHITE:
			if (x >= width || y >= height || x < 0 || y < 0)
				return Vecf.onesLike(pixels[0][0]);
			break;
		case BLACK:
			if (x >= width || y >= height || x < 0 || y < 0)
				return Vecf.zeroesLike(pixels[0][0]);
			break;
		}
		return pixels[x][y];
	}
	
	/**
	 * Returns an interpolated sample of the image using the input floating-point coordinates.
	 * <p>
	 * The interpolation mode used is determined by the image's <code>filterMode</code>.
	 * Edge modes will also affect the output if sampling near the edges of the image.
	 * @param x
	 * @param y
	 * @return A Vecf obtained by reconstructing the imaged object at (x, y)
	 */
	public Vecf sample(float x, float y) {
		x -= 0.5f;
		y -= 0.5f;
		switch (filterMode) {
		case NEAREST:
			x += 0.5f;
			y += 0.5f;
			return getPixel((int) x, (int) y);
		case BILINEAR:
			float dx = (x+1) % 1;
			float dy = (y+1) % 1;

			return Vecf.add(
					Vecf.add(Vecf.mult(getPixel((int) Math.floor(x), (int) Math.floor(y)), (1 - dx) * (1 - dy)),
							Vecf.mult(getPixel((int) Math.floor(x + 1), (int) Math.floor(y)), dx * (1 - dy))),
					Vecf.add(Vecf.mult(getPixel((int) Math.floor(x), (int) Math.floor(y + 1)), (1 - dx) * dy),
							Vecf.mult(getPixel((int) Math.floor(x + 1), (int) Math.floor(y + 1)), dx * dy)));
		case BICUBIC:
			return sampleCubic(x, y);
		case LANCZOS:
			return sampleLanczos(x, y, LANCZOS_RADIUS);
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
			case PHOTOSHOP_BICUBIC:
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
		// get relevant samples
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
		// combine
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
	
	Vecf sampleLanczos(float x, float y, int a) {
		// TODO: add case for sampling directly on the point
		// maybe do that everywhere, lol
		int x_floor = (int)Math.floor(x);
		int y_floor = (int)Math.floor(y);
		
		// The 2D Lanczos kernel is, by def., separable
		// this speeds things up by a factor of (a)
		float[] LUTX = new float[2*a];
		float[] LUTY = new float[2*a];
		for(int t = 1-a; t <= a; t++) {
			LUTX[t+a-1] = PixelMath.lanczos(x_floor - x + t, a);
			LUTY[t+a-1] = PixelMath.lanczos(y_floor - y + t, a);
		}
		
		float weight_sum = 0.f;
		Vecf sample_sum = new Vecf(this.channels);
		for(int ix = x_floor + 1-a; ix <= x_floor + a; ix++) {
			for(int iy = y_floor + 1-a; iy <= y_floor + a; iy++) {
				float s_i = LUTX[ix - (x_floor + 1-a)]*LUTY[iy - (y_floor + 1-a)];
				sample_sum.add(Vecf.mult(getPixel(ix, iy), s_i));
				weight_sum += s_i;
			}
		}
		return Vecf.mult(sample_sum, 1.f/weight_sum);
	}

	private void welcome() {
		System.out.println("##library.name## ##library.prettyVersion## by ##author##");
	}
}
