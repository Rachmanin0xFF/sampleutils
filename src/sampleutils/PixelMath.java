package sampleutils;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

public final class PixelMath {

	static final float PI = 3.14159265359f;
	static final float TWO_PI = 6.28318530718f;
	static final float PHI = 1.61803398875f;

	public static float lanczos(float x, float a) {
		if (x < -a || x > a)
			return 0;
		if (x == 0.0)
			return 1;
		return (float) (a * Math.sin(PI * x) * Math.sin(PI * x / a) / (PI * PI * x * x));
	}

	public static float smoothstep(float x) {
		if (x < 0)
			return 0;
		if (x > 1)
			return 1;
		return 3 * x * x - 2 * x * x * x;
	}

	public static float smootherstep(float x) {
		if (x < 0)
			return 0;
		if (x > 1)
			return 1;
		return 6 * x * x * x * x * x - 15 * x * x * x * x + 10 * x * x * x;
	}

	public static int mod(int x, int n) {
		return ((x >> 31) & (n - 1)) + (x % n);
	}

	public static int clamp_inclusive(int x, int low, int high) {
		if (x <= low)
			return low;
		if (x >= high)
			return high;
		return x;
	}

	// for internal use -- I don't want this class to rely on Processing's
	// libraries, so we're rewriting this
	private static float random(float max) {
		return random(0, max);
	}

	private static float random(float min, float max) {
		return (float) (Math.random() * (max - min) + min);
	}

	private static int round(float x) {
		return (int) (x + 0.5);
	}

	// modified sunflower function
	// https://stackoverflow.com/a/28572551
	public static Vecf[] sunflower(int n) {
		return sunflower(n, 0.75f);
	}

	// higher alpha = more constant edge, default is 0.75
	public static Vecf[] sunflower(int n, float alpha) {
		Vecf[] o = new Vecf[n];
		int b = round(alpha * (float) Math.sqrt(n));
		for (int i = 0; i < n; i++) {
			float r = 1.0f;
			if (i < n - b)
				r = (float) (Math.sqrt((float) i - 0.5) / Math.sqrt(n - (b + 1.0) * 0.5));
			float theta = 2.0f * PI * (float) i / (PHI * PHI);
			o[i] = new Vecf(r * (float) Math.cos(theta), r * (float) Math.sin(theta));
		}
		return o;
	}

	// r - grid radius, k - sample attempts (default 20)
	// returns a disc
	Vecf[] poisson(int r, int k) {
		ArrayList<Vecf> dead = new ArrayList<Vecf>();
		HashSet<Vecf> active = new HashSet<Vecf>(); // O(1) removal time
		boolean[][] occupied = new boolean[r * 2 + 1][r * 2 + 1];
		Vecf[][] coord_grid = new Vecf[occupied.length][occupied[0].length];

		float x0 = random(-r / 2, r / 2);
		float y0 = random(-r / 2, r / 2);

		active.add(new Vecf(x0, y0));
		occupied[round(x0 + r)][round(y0 + r)] = true;
		coord_grid[round(x0 + r)][round(y0 + r)] = new Vecf(x0, y0);

		ArrayList<Vecf> new_active = new ArrayList<Vecf>();
		Iterator<Vecf> it = active.iterator();
		while (!(active.isEmpty() && new_active.isEmpty())) {
			while (it.hasNext()) {
				Vecf p = it.next();

				int i;
				float sd = random(TWO_PI);
				for (i = 0; i < k; i++) {
					float theta = sd + TWO_PI * ((float) i / (float) k);

					// the "maximal Poisson sampling" algorithm described here:
					// http://extremelearning.com.au/an-improved-version-of-bridsons-algorithm-n-for-poisson-disc-sampling/
					// produces veeeery minor rivers in the points if theta is non-uniformly
					// distributed (and maybe if it is, still?)
					// I counteract this by adding a small random value to the radius. There's still
					// some anisotropy, but it looks better, I think.
					float rd = 1.41421356237f + random(0.15f);
					float x = p.components[0] + (float) Math.cos(theta) * rd;
					float y = p.components[1] + (float) Math.sin(theta) * rd;
					int xc = round(x + r);
					int yc = round(y + r);

					if (xc > 0 && yc > 0 && xc < occupied.length - 1 && yc < occupied.length - 1) {
						boolean tooClose = false;

						if (x * x + y * y > r * r)
							tooClose = true;

						if (!tooClose)
							for (int xf = -1; xf <= 1; xf++) {
								for (int yf = -1; yf <= 1; yf++) {
									if (occupied[xc + xf][yc + yf]) {
										Vecf s = coord_grid[xc + xf][yc + yf];
										if (Vecf.magnitudeSquared(Vecf.sub(s, new Vecf(x, y))) <= 2.0) {
											tooClose = true;
											xf = 2;
											break;
										}
									}
								}
							}
						if (!tooClose) {
							new_active.add(new Vecf(x, y));
							dead.add(new Vecf(x, y));

							occupied[xc][yc] = true;
							coord_grid[xc][yc] = new Vecf(x, y);

							break;
						}
					}
				}
				if (i == k) {
					it.remove();
					dead.add(p);
				}
			}
			active.addAll(new_active);
			it = active.iterator();
			new_active.clear();
		}

		Vecf[] o = new Vecf[dead.size()];
		for (int i = 0; i < o.length; i++) {
			o[i] = Vecf.mult(dead.get(i), 1.f / (float) r);
		}

		return o;
	}
}
