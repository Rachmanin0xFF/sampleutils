package sampleutils;

public class Vecf {
	float[] components;
	
	public Vecf(int dimension) {
		components = new float[dimension];
	}
	public Vecf(float _x, float _y, float _z) {
		components = new float[] {_x, _y, _z};
	}
	public Vecf(float _x, float _y, float _z, float _w) {
		components = new float[] {_x, _y, _z, _w};
	}
	public static Vecf add(Vecf a, Vecf b) {
		if(a.components.length != b.components.length) {
			System.err.println("Error! Trying to add vectors of length " + a.components.length + " " + b.components.length);
		}
		Vecf sum = new Vecf(a.components.length);
		for(int i = 0; i < a.components.length; i++) {
			sum.components[i] = a.components[i] + b.components[i];
		}
		return sum;
	}
	public static Vecf sub(Vecf a, Vecf b) {
		if(a.components.length != b.components.length) {
			System.err.println("Error! Trying to add vectors of length " + a.components.length + " " + b.components.length);
		}
		Vecf sum = new Vecf(a.components.length);
		for(int i = 0; i < a.components.length; i++) {
			sum.components[i] = a.components[i] - b.components[i];
		}
		return sum;
	}
	public static Vecf mult(Vecf a, float c) {
		Vecf sum = new Vecf(a.components.length);
		for(int i = 0; i < a.components.length; i++) {
			sum.components[i] = a.components[i]*c;
		}
		return sum;
	}
}
