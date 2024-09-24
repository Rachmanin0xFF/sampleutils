package sampleutils;

/**
 * static VecImage filters; includes up/downscaling
 */
public final class Filters {
	private Filters() {
	}
	
	/**
	 * Convolves a VecImage with a kernel. The kernel may be unnormalized, but it must have odd dimensions.
	 * @param input Any VecImage.
	 * @param kernel A (2*n+1) by (2*m+1) array of floats.
	 * @return The VecImage <a href="https://en.wikipedia.org/wiki/Convolution#Discrete_convolution">convolved</a> with the kernel.
	 */
	public static VecImage convolve(VecImage input, float[][] kernel) {
		// conveniently, this also checks that the lengths aren't zero
		if(kernel.length % 2 == 0 || kernel[0].length%2 == 0) {
			throw new IllegalArgumentException("Kernel dimensions must be odd numbers!");
		}
		int radiusX = (kernel.length-1)/2;
		int radiusY = (kernel.length-1)/2;
		VecImage output = new VecImage(input.width, input.height, input.channels);
		
		for(int x = 0; x < output.width; x++) {
			for(int y = 0; y < output.height; y++) {
				for(int i = -radiusX; i < radiusX; i++) {
					for(int j = -radiusY; j < radiusY; j++) {
						output.pixels[x][y].add(Vecf.mult(input.getPixel(x+i, y+j), kernel[i+radiusX][j+radiusY]));
					}
				}
			}
		}
		return output;
	}
	
	/**
	 * Upscales an image using its edge and sample mode properties.
	 * @param input the VecImage to upscale
	 * @param new_width
	 * @param new_height
	 * @return the upscaled VecImage
	 */
	public static VecImage upscale(VecImage input, int new_width, int new_height) {
		VecImage upscaled = new VecImage(new_width, new_height, input.channels);
		for (int x = 0; x < new_width; x++)
			for (int y = 0; y < new_height; y++) {
				upscaled.pixels[x][y] = input.sample((x / (float) new_width) * (float) input.width,
						(y / (float) new_height) * (float) input.height);
			}
		return upscaled;
	}
	
	/**
	 * Inverts an image given a white level.
	 * <p>
	 * Each component of each pixel is replaced with <code>(white_level - component)</code>.
	 * @param input
	 * @param white_level
	 * @return the inverted VecImage
	 */
	public static VecImage invert(VecImage input, float white_level) {
		VecImage output = new VecImage(input);
		for (int x = 0; x < output.width; x++)
			for (int y = 0; y < output.height; y++) {
				output.pixels[x][y] = Vecf.sub(Vecf.mult(Vecf.onesLike(output.pixels[x][y]), white_level),
						input.pixels[x][y]);
			}
		return output;
	}
	
	/**
	 * Calculates the medioid (higher-dimensional median analog) of a VecImage.
	 * <p>
	 * The medioid of an N-dimensional list of vectors is the vector with the
	 * smallest summed distance to every other vector in the list. This method uses the
	 * <a href="http://proceedings.mlr.press/v54/newling17a/newling17a.pdf">trimed</a> algorithm
	 * to compute the mediod, which runs in O(N^(3/2)).
	 * @param input
	 * @param rad the (circular) radius to take the medioid over
	 * @return the medioid-filtered image
	 */
	public static VecImage mediod(VecImage input, int rad) {
		VecImage output = new VecImage(input);
		boolean[][] kernel = new boolean[rad * 2 + 1][rad * 2 + 1];
		int samp_count = 0;
		for (int x = -rad; x <= rad; x++) {
			for (int y = -rad; y <= rad; y++) {
				kernel[x + rad][y + rad] = x * x + y * y <= rad * rad;
				if (kernel[x + rad][y + rad])
					samp_count++;
			}
		}
		Vecf[] samples = new Vecf[samp_count];
		int[][] coords = new int[samp_count][2];
		int i = 0;
		for (int x = -rad; x <= rad; x++) {
			for (int y = -rad; y <= rad; y++) {
				if (kernel[x + rad][y + rad]) {
					coords[i][0] = x;
					coords[i][1] = y;
					i++;
				}
			}
		}
		for (int x = 0; x < output.width; x++) {
			for (int y = 0; y < output.height; y++) {
				for (i = 0; i < coords.length; i++) {
					samples[i] = input.getPixel(coords[i][0] + x, coords[i][1] + y);
				}
				output.pixels[x][y] = medoid(samples);
			}
		}
		return output;
	}

	// This function selects the Vecf in the input array that best
	// represents the dataset as a whole (minimizing absolute error)
	// This function uses the "trimed" algorithm:
	// http://proceedings.mlr.press/v54/newling17a/newling17a.pdf
	// TODO: implement a more efficient mediod algorithm
	private static Vecf medoid(Vecf[] data) {
		// lower bound on the energy of each entry
		float[] lowers = new float[data.length];
		// used in loop
		float[] distances = new float[data.length];

		// minimum energy value and index
		float min_E = Float.MAX_VALUE;
		int min_E_i = 0;
		for (int i = 0; i < data.length; i++) {
			if (lowers[i] < min_E) {
				// compute distances
				float E_i = 0.f;
				for (int j = 0; j < data.length; j++) {
					if (i != j) {
						distances[j] = Vecf.magnitudeSquared(Vecf.sub(data[i], data[j]));
						E_i += distances[j];
					}
				}

				if (E_i < min_E) {
					// update medoid
					min_E = E_i;
					min_E_i = i;
				}

				// update lower bounds
				for (int j = i + 1; j < data.length; j++) {
					// multiply distances[j] by length since we don't normalize error
					lowers[j] = Math.max(lowers[j], (float) Math.abs(data.length * distances[j] - E_i));
				}
			}
		}

		return data[min_E_i];
	}

	// AdvMAME2x pixel art scaling algorithm
	// uses
	/**
	 * Scales an input image by 2x. Best for upscaling pixel art.
	 * <p>
	 * The <a href="https://en.wikipedia.org/wiki/Pixel-art_scaling_algorithms#EPX/Scale2%C3%97/AdvMAME2%C3%97">EPX2</a>
	 * algorithm was developed by Eric Johnston in 1992 as a way to upscale
	 * classic games to higher-resolution monitors. The algorithm works well on small photos
	 * with blobs of same-color pixels.
	 * @param input
	 * @return
	 */
	public static VecImage EPX2(VecImage input) {
		VecImage o = new VecImage(input.width * 2, input.height * 2, input.channels);
		for (int x = 0; x < input.width; x++) {
			for (int y = 0; y < input.height; y++) {
				Vecf P = input.getPixel(x, y);
				Vecf A = input.getPixel(x, y - 1);
				Vecf B = input.getPixel(x + 1, y);
				Vecf C = input.getPixel(x - 1, y);
				Vecf D = input.getPixel(x, y + 1);

				Vecf V1 = new Vecf(P);
				Vecf V2 = new Vecf(P);
				Vecf V3 = new Vecf(P);
				Vecf V4 = new Vecf(P);

				if (Vecf.equals(C, A) && !Vecf.equals(C, D) && !Vecf.equals(A, B))
					V1 = A;
				if (Vecf.equals(A, B) && !Vecf.equals(A, C) && !Vecf.equals(B, D))
					V2 = B;
				if (Vecf.equals(D, C) && !Vecf.equals(D, B) && !Vecf.equals(C, A))
					V3 = C;
				if (Vecf.equals(B, D) && !Vecf.equals(B, A) && !Vecf.equals(D, C))
					V4 = D;

				o.pixels[x * 2][y * 2] = V1;
				o.pixels[x * 2 + 1][y * 2] = V2;
				o.pixels[x * 2][y * 2 + 1] = V3;
				o.pixels[x * 2 + 1][y * 2 + 1] = V4;
			}
		}
		return o;
	}

}
