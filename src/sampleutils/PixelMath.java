package sampleutils;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

/**
 * Algorithms and math for sampling, color, and other image-related utilities.
 */
public final class PixelMath {

	static final float PI = 3.14159265359f;
	static final float TWO_PI = 6.28318530718f;
	static final float PHI = 1.61803398875f;

	private PixelMath() {
	} // static use only
	
	/**
	 * Normalizes a 2D array of floats so that its entries sum to one.
	 * @param input an input array of floats
	 * @return The normalized array.
	 */
	public static float[][] normalize(float[][] input) {
		float[][] output = new float[input.length][input[0].length];
		float sum = 0.f;
		for(int x = 0; x < input.length; x++) for(int y = 0; y < input[0].length; y++) {
			sum += input[x][y];
		}
		if(sum == 0) throw new IllegalArgumentException("Input array sums to zero! Cannot normalize!");
		for(int x = 0; x < input.length; x++) for(int y = 0; y < input[0].length; y++) {
			output[x][y] = input[x][y] / sum;
		}
		return output;
	}
	
	/**
	 * Generates a normalized 2D Gaussian kernel.
	 * 
	 * @param stdev The standard deviation of the Gaussian.
	 * @param stdevs_to_include The extent of the Gaussian contained in the kernel, in units of standard deviation.
	 * @return A two-dimensional array of floating-point values that sum to one.
	 */
	public static float[][] gaussianKernel(float stdev, float stdevs_to_include) {
		if (stdev <= 0 || stdevs_to_include <= 0) {
			throw new IllegalArgumentException("Gaussian kernel input parameters cannot be <= 0!");
		}

		// generate separable part
		int sample_radius = (int) Math.ceil(stdevs_to_include * stdev);
		float[] kernelX = new float[2 * sample_radius + 1];

		for (int x = -sample_radius; x <= sample_radius; x++) {
			kernelX[x + sample_radius] = unnormalizedGaussian(x, stdev);
		}
		// combine
		float[][] outputKernel = new float[kernelX.length][kernelX.length];
		float sum = 0.f;
		for (int x = 0; x <= kernelX.length; x++)
			for (int y = 0; y <= kernelX.length; y++) {
				outputKernel[x][y] = kernelX[x] * kernelX[y];
				sum += outputKernel[x][y];
			}
		// normalize
		for (int x = 0; x <= kernelX.length; x++)
			for (int y = 0; y <= kernelX.length; y++) {
				outputKernel[x][y] /= sum;
			}
		return outputKernel;

	}
	
	/**
	 * Generates a normalized bivariate Gaussian kernel.
	 * <p>
	 * In addition to standard deviations in the x and y direction, this function also accepts a Pearson correlation parameter
	 * (higher values 'squish' the Gaussian along the diagonal). Together, these three values are sufficient to parameterize
	 * the space of symmetric covariance matrices (stdevx^2 and stdevy^2 on the diagonal, rho*stdevx*stdevy on the off-diagonal).
	 * <p>
	 * For axis-aligned distributions (i.e. vertical and horizontal blurs), rho should be zero.
	 * 
	 * @param stdevx The standard deviation of the distribution in the x direction.
	 * @param stdevy The standard deviation of the distribution in the y direction.
	 * @param rho The Pearson correlation of the Gaussian distribution.
	 * @param stdevs_to_include The extent of the Gaussian contained in the kernel, in units of standard deviation.
	 * @return A two-dimensional array of floating-point values that sum to one.
	 */
	public static float[][] bivariateGaussianKernel(float stdevx, float stdevy, float rho, float stdevs_to_include) {
		if (stdevx <= 0 || stdevy <= 0 || stdevs_to_include <= 0) {
			throw new IllegalArgumentException("Gaussian kernel input parameters cannot be <= 0!");
		}
		
		// while this kernel is still separable with a change of basis, it is not separable
		// on our discrete grid, so we need to do everything at once
		int sample_radius = (int) Math.ceil(stdevs_to_include * Math.max(stdevx, stdevy));
		float[][] kernel = new float[2 * sample_radius + 1][2 * sample_radius + 1];
		
		float sum = 0.f;
		for (int x = -sample_radius; x <= sample_radius; x++)
			for (int y = -sample_radius; y <= sample_radius; y++) {
				kernel[x][y] = unnormalizedBivariateGaussian(x, y, stdevx, stdevy, rho);
				sum += kernel[x][y];
			}
		// normalize
		for (int x = 0; x <= kernel.length; x++)
			for (int y = 0; y <= kernel.length; y++) {
				kernel[x][y] /= sum;
			}
		return kernel;

	}

	private static float unnormalizedBivariateGaussian(float x, float y, float stdevx, float stdevy, float rho) {
		float factor = -1.f/(2.f*(1.f - rho*rho));
		factor *= x*x/(stdevx*stdevx) - 2*rho*x*y/(stdevx*stdevy) + y*y/(stdevy*stdevy);
		return (float) Math.exp(factor);
	}

	private static float unnormalizedGaussian(float x, float stdev) {
		return (float) Math.exp(-x * x / (2.f * stdev * stdev));
	}

	/**
	 * Similar to Processing's map() function.
	 * <p>
	 * Maps x from the range (in_min, in_max) to (out_min, out_max) without
	 * clamping. The output is precisely
	 * <code>out_min + (out_max-out_min)*(x - in_min)/(in_max - in_min)</code>.
	 * 
	 * @param x       The number to map
	 * @param in_min  The minimum value of x
	 * @param in_max  The maximum value of x
	 * @param out_min The minimum output value
	 * @param out_max The maximum output value
	 * @return x mapped to the new range
	 */
	public static float map(float x, float in_min, float in_max, float out_min, float out_max) {
		return out_min + (out_max - out_min) * (x - in_min) / (in_max - in_min);
	}

	/**
	 * Computes a 1-dimensional Lanczos kernel.
	 * 
	 * @param x The sample location.
	 * @param a The width of the kernel (typically an integer).
	 * @return The value of L(x)
	 */
	public static float lanczos(float x, float a) {
		if (x < -a || x > a)
			return 0;
		if (x == 0.0)
			return 1;
		return (float) (a * Math.sin(PI * x) * Math.sin(PI * x / a) / (PI * PI * x * x));
	}
	
	/**
	 * The unnormalized Laplacian of the 2D Gaussian function.
	 * @param x
	 * @param y
	 * @param stdev The standard deviation of the Gaussian.
	 * @return
	 */
	float LoG(float x, float y, float stdev) {
	  float xp = -(x*x+y*y)/(2*stdev*stdev);
	  return -(1.f+xp)*(float)Math.exp(xp)/(PI*stdev*stdev); // NOT NORMALIZED PROPERLY FOR CONSTANT EDGE DETECTION VALUES WITH VARYING DEVIATIONS!!!
	}

	/**
	 * A map from [0,1] onto [0,1] with first derivatives equal to zero at 0 and 1.
	 * 
	 * @param x
	 * @return 3x^2-2x^3
	 */
	public static float smoothstep(float x) {
		if (x < 0)
			return 0;
		if (x > 1)
			return 1;
		return 3 * x * x - 2 * x * x * x;
	}

	/**
	 * A map from [0,1] onto [0,1] with first and second derivatives equal to zero
	 * at 0 and 1.
	 * 
	 * @param x
	 * @return 6x^5-15x^4+10x^3
	 */
	public static float smootherstep(float x) {
		if (x < 0)
			return 0;
		if (x > 1)
			return 1;
		return 6 * x * x * x * x * x - 15 * x * x * x * x + 10 * x * x * x;
	}

	// this can handle negatives (x>>31 checks if the int is negative)
	/**
	 * Similar to java's modulo, but wraps neatly into the negatives.
	 * <p>
	 * For example, mod(-1, 3) = 2. While this is not technically how modulo is
	 * expected to operate, it is useful for texture wrapping.
	 * 
	 * @param x
	 * @param n
	 * @return modulo with support for negative numbers.
	 */
	public static int mod(int x, int n) {
		return ((x >> 31) & n) + (x % n);
	}

	/**
	 * Clamps the input between low and high.
	 * 
	 * @param x
	 * @param low
	 * @param high
	 * @return (low if x &lt;= low), (high if x &lt;= high), otherwise x
	 */
	public static int clamp_inclusive(int x, int low, int high) {
		if (x <= low)
			return low;
		if (x >= high)
			return high;
		return x;
	}

	// for internal use -- I'd like this class to be decoupled from Processing's
	// libraries, so we're rewriting this
	private static float random(float max) {
		return random(0, max);
	}

	private static float random(float min, float max) {
		return (float) (Math.random() * (max - min) + min);
	}

	private static int round(float x) {
		return (int) (x + 0.5f);
	}

	/**
	 * Returns n 2D Vecfs evenly distributed in a disk with radius 1.
	 * <p>
	 * Uses a modified
	 * <a href="https://stackoverflow.com/a/28572551">'sunflower'</a> function.
	 * Calls <code>sunflower(n, 0.75)</code>.
	 * 
	 * @param n the number of points to return
	 * @return an array of n 2D Vecfs with magnitude &lt;= 1
	 */
	public static Vecf[] sunflower(int n) {
		return sunflower(n, 0.75f);
	}

	// higher alpha = more constant edge, default is 0.75
	/**
	 * Returns n 2D Vecfs evenly distributed in a disk with radius 1.
	 * <p>
	 * Uses a modified
	 * <a href="https://stackoverflow.com/a/28572551">'sunflower'</a> function.
	 * 
	 * @param n     the number of points to return
	 * @param alpha the edge uniformity. A higher alpha means a smoother edge.
	 * @return an array of n 2D Vecfs with magnitude &lt;= 1
	 */
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
	// returns an evenly-spaced anisotropic set of points in the disk r < 1 centered
	// at the origin
	/**
	 * Generates an evenly-spaced anisotropic array of points in a disk with radius
	 * 1.
	 * <p>
	 * This is based on the <a href=
	 * "http://extremelearning.com.au/an-improved-version-of-bridsons-algorithm-n-for-poisson-disc-sampling/">
	 * "maximal poisson disk sampling"</a> algorithm. Typically, this produces very
	 * small rivers in the points if theta is non-uniformly distributed (and maybe
	 * in general?). This is counteracted here by adding a small random value to the
	 * radius.
	 * 
	 * @param r the radius of the hashing grid
	 * @param k the number of sample attempts (recommended value is 20)
	 * @return A set of
	 */
	public static Vecf[] poisson(int r, int k) {
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
	
	private static final float RAD_TO_DEG = 57.2957795131f;
	private static final float DEG_TO_RAD = 0.01745329251f;

	//============================================= COLOR SPACES =====================================================//

	// Descriptions of the color spaces listed here:
	// NAME    sRGB GAMUT RANGE
	// sRGB:   0 to 255 
	// RGB:    0 to 1 
	// XYZ:    0 to (95, 100, 108) 
	// L*a*b*: (0, -86, -108) to 100          
	// LCH_ab: (0, 0, -180) to (100, 134, 180)
	// LCH_uv: (0, 0, -180) to (100, 179, 180) 
	// ΔE2000: 0 to 162                         

	// All conversions use the CIE D65 Standard Illuminant unless specified otherwise
	public static Vecf ref_white = new Vecf(95.047f, 100.0f, 108.883f);
	
	/**
	 * Fairchild-Pirrotta lightness (helmholtz-kohlrausch compensation) of an sRGB color.
	 * @param c The color in sRGB space (components 0-255)
	 * @return The lightness of the color
	 */
	public static float lightnessFPsRGB(Vecf c) {
	  return lightnessFP(LAB_to_LCH(XYZ_to_LAB(RGB_to_XYZ(sRGB_to_RGB(c)))));
	}
	
	/**
	 * Converts an sRGB color to LCH coordinates.
	 * @param c The color in sRGB.
	 * @return The same color in LCH.
	 */
	public static Vecf sRGB_to_LCH(Vecf c) {
	  return LAB_to_LCH(XYZ_to_LAB(RGB_to_XYZ(sRGB_to_RGB(c))));
	}

	// fairchild pirrotta lightness (helmholtz-kohlrausch compensation)
	// thanks to https://stackoverflow.com/questions/596216/formula-to-determine-brightness-of-rgb-color/59602392#59602392
	/**
	 * Fairchild-Pirrotta lightness (helmholtz-kohlrausch compensation) of an LCH color.
	 * @param c The color in LCH space
	 * @return The lightness of the color
	 */
	public static float lightnessFP(Vecf lch) {
	  return lch.components[0] +(2.5f - 0.025f*lch.components[0])
	               *(0.116f*(float)Math.abs((float)Math.sin(DEG_TO_RAD*((lch.components[2]-90.0f)/2.0f))) + 0.085f)
	               *lch.components[1];
	}

	public static Vecf LAB_to_LCH(Vecf c) {
	  return new Vecf(c.components[0], (float)Math.sqrt(c.components[1]*c.components[1] + c.components[2]*c.components[2]), (float)Math.atan2(c.components[2], c.components[1])*RAD_TO_DEG);
	}
	public static Vecf LUV_to_LCH(Vecf c) {
	  return LAB_to_LCH(c);
	}

	// LAB is L*a*b*
	public static Vecf XYZ_to_LAB(Vecf c) {
	  Vecf v = new Vecf(c);
	  
	  // D65 white point
	  v.components[0] /= ref_white.components[0];
	  v.components[1] /= ref_white.components[1];
	  v.components[2] /= ref_white.components[2];
	  
	  v.components[0] = v.components[0] > 0.008856f ? (float)Math.pow(v.components[0], 1.0f/3.0f) : v.components[0]*7.787f + 16.0f/116.0f;
	  v.components[1] = v.components[1] > 0.008856f ? (float)Math.pow(v.components[1], 1.0f/3.0f) : v.components[1]*7.787f + 16.0f/116.0f;
	  v.components[2] = v.components[2] > 0.008856f ? (float)Math.pow(v.components[2], 1.0f/3.0f) : v.components[2]*7.787f + 16.0f/116.0f;
	  
	  return new Vecf(116.0f*v.components[1] - 16.0f, 500.0f*(v.components[0] - v.components[1]), 200.0f*(v.components[1] - v.components[2]));
	}

	// LUV is L*u*v*
	public static Vecf XYZ_to_LUV(Vecf c) {
	  float v_rp = 9.0f*ref_white.components[1]/(ref_white.components[0] + 15.0f*ref_white.components[1] + 3.0f*ref_white.components[2]);
	  float u_rp = 4.0f*ref_white.components[0]/(ref_white.components[0] + 15.0f*ref_white.components[1] + 3.0f*ref_white.components[2]);
	  
	  float vp =  9.0f*c.components[1]/(c.components[0] + 15.0f*c.components[1] + 3.0f*c.components[2]);
	  float up =  4.0f*c.components[0]/(c.components[0] + 15.0f*c.components[1] + 3.0f*c.components[2]);
	  
	  float y_r =  c.components[1]/ref_white.components[1];
	  float L = y_r > 0.008856f ? 116.0f*(float)Math.pow(y_r, 1.0f/3.0f) - 16.0f : 903.3f*y_r;
	  
	  float u = 13.0f*L*(up - u_rp);
	  float v = 13.0f*L*(vp - v_rp);
	  
	  return new Vecf(L, u, v);
	}

	public static Vecf RGB_to_XYZ(Vecf c) {
	  return new Vecf(41.24f * c.components[0] + 35.76f * c.components[1] + 18.05f * c.components[2],
	                     21.26f * c.components[0] + 71.52f * c.components[1] + 07.22f * c.components[2],
	                     01.93f * c.components[0] + 11.92f * c.components[1] + 95.05f * c.components[2]);
	}

	public static Vecf sRGB_to_RGB(Vecf c) {
	  return new Vecf(sX_to_X(c.components[0]/255.0f),
	                     sX_to_X(c.components[1]/255.0f),
	                     sX_to_X(c.components[2]/255.0f));
	}

	// sRGB to RGB gamma transformation
	private static float sX_to_X(float x) {
	  if(x <= 0.04045f) return x/12.92f;
	  else return (float)Math.pow((x + 0.055f)/1.055f, 2.4f);
	}
	
	/**
	 * The [delta-E 2000](http://www2.ece.rochester.edu/~gsharma/ciede2000/ciede2000noteCRNA.pdf) color distance metric. 
	 * @param lab1 The first color in sRGB space (0-255)
	 * @param lab2 The second color in sRGB space (0-255)
	 * @return The perceptual distance between the two colors (as approximated by the delta-E 2000 standard)
	 */
	public static float DELTA_E_2000_sRGB(Vecf a, Vecf b) {
	  return DELTA_E_2000(XYZ_to_LAB(RGB_to_XYZ(sRGB_to_RGB(a))), XYZ_to_LAB(RGB_to_XYZ(sRGB_to_RGB(b))));
	}

	// this was a pain
	// math from:
	// http://www2.ece.rochester.edu/~gsharma/ciede2000/ciede2000noteCRNA.pdf
	// http://zschuessler.github.io/DeltaE/learn/
	
	/**
	 * The [delta-E 2000](http://www2.ece.rochester.edu/~gsharma/ciede2000/ciede2000noteCRNA.pdf) color distance metric. 
	 * @param lab1 The first color in L*a*b* space
	 * @param lab2 The second color in L*a*b* space
	 * @return The perceptual distance between the two colors (as approximated by the delta-E 2000 standard)
	 */
	public static float DELTA_E_2000(Vecf lab1, Vecf lab2) {
	  float C_ab = ((float)Math.sqrt(lab1.components[1]*lab1.components[1] + lab1.components[2]*lab1.components[2]) + (float)Math.sqrt(lab2.components[1]*lab2.components[1] + lab2.components[2]*lab2.components[2]))/2.0f;
	  float G = 0.5f*(1-(float)Math.sqrt((float)Math.pow(C_ab, 7.0f)/((float)Math.pow(C_ab, 7f) + 6103515625.0f)));
	  float ap1 = (1+G)*lab1.components[1];
	  float ap2 = (1+G)*lab2.components[1];
	  float Cp1 = (float)Math.sqrt(ap1*ap1 + lab1.components[2]*lab1.components[2]);
	  float Cp2 = (float)Math.sqrt(ap2*ap2 + lab2.components[2]*lab2.components[2]);
	  
	  float h1 = (float)Math.atan2(lab1.components[2], ap1)*RAD_TO_DEG;
	  float h2 = (float)Math.atan2(lab2.components[2], ap2)*RAD_TO_DEG;
	  
	  float del_Lp = lab2.components[0] - lab1.components[0];
	  float del_Cp = Cp2 - Cp1;
	  
	  float del_hp = 0f;
	  if(Cp1*Cp2 != 0.0f) {
	    if((float)Math.abs(h2 - h1) <= 180f) del_hp = h2-h1;
	    else if(h2 - h1 > 180f)  del_hp = h2-h1 - 360f;
	    else if(h2 - h1 < -180f) del_hp = h2-h1 + 360f;
	  }
	  
	  float del_H = 2*(float)Math.sqrt(Cp1*Cp2)*(float)Math.sin(DEG_TO_RAD*(del_hp/2.0f));
	  
	  float Hp = (float)Math.abs(h1-h2) > 180f ? (h1 + h2 + 360.0f)/2.0f :
	                                (h1 + h2)/2.0f;
	  float T = 1f - 0.17f*(float)Math.cos(DEG_TO_RAD*(Hp - 30.0f))
	              + 0.24f*(float)Math.cos(DEG_TO_RAD*(2.0f*Hp))
	              + 0.32f*(float)Math.cos(DEG_TO_RAD*(3.0f*Hp + 6.0f))
	              - 0.2f*(float)Math.cos(DEG_TO_RAD*(4f*Hp - 63.0f));
	  
	  float Lb = (lab1.components[1]+lab2.components[1])/2.0f;
	  float Cp = (Cp1 + Cp2)/2.0f;
	  
	  float S_L = 1f + 0.015f*(Lb-50.0f)*(Lb-50.0f)/(float)Math.sqrt(20.0f + (Lb-50.0f)*(Lb-50.0f));
	  float S_C = 1f + 0.045f*Cp;
	  float S_H = 1f + 0.015f*Cp*T;
	  
	  float RT = -2f*(float)Math.sqrt((float)Math.pow(Cp, 7f)/((float)Math.pow(Cp, 7f) + 6103515625.0f)*(float)Math.sin(DEG_TO_RAD*60.0f*(float)Math.exp(-((Hp - 275.0f)*(Hp - 275.0f)/625.0f))));
	  
	  float E00 = (float)Math.sqrt((del_Lp/S_L)*(del_Lp/S_L) + (del_Cp/S_C)*(del_Cp/S_C) + (del_H)/S_H*(del_H)/S_H + RT*del_Cp*del_H/(S_C*S_H));
	  
	  return E00;
	}
}
