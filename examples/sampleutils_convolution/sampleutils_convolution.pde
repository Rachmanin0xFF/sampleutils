import sampleutils.*;

SampleUtilsIO sio;
VecImage v;

void setup() {
  size(1200, 440, P2D);
  background(255, 0, 255);
  noSmooth();
  
  // Give sampleutils access to the sketch's PApplet
  // (it needs this to call createImage / loadImage, which both reserve GPU memory)
  sio = new SampleUtilsIO(this);
  v = sio.loadImage("pyramid.png", false);
}

void draw() {
  background(0);
  fill(255);
  
  v.edgeMode = VecImage.EdgeModes.CLAMP;
  v.filterMode = VecImage.FilterModes.NEAREST;
	
  textAlign(CENTER, CENTER);
  textSize(15);

  float[][] k = PixelMath.differenceOfGaussiansKernel(8, 3, 2.5);
  k = PixelMath.multKernel(k, 10.f);
  
  image(sio.toPImage(v), 0, 0);
  text("Original", 200, 415);
  translate(400, 0);

  image(sio.toPImage(Filters.upscale(new VecImage(PixelMath.multKernel(k, 10.0*k.length*k[0].length)), 400, 400)), 0, 0);
  text("Kernel", 200, 415);
  translate(400, 0);

  image(sio.toPImage(Filters.convolve(v, k)), 0, 0);
  text("Convolved", 200, 415);
  translate(400, 0);

  saveFrame("convolution_example_output.png");
  noLoop();
}
