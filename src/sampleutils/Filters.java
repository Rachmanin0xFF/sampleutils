package sampleutils;

public final class Filters {
	private Filters() {
	}

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

	// This function selects the PVector in the input array that best
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

}
