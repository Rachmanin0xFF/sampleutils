import sampleutils.*;

SampleUtilsIO sio;
VecImage v;

void setup() {
  size(1000, 240, P2D);
  background(255, 0, 255);
  noSmooth();
  
  // Give sampleutils access to the sketch's PApplet
  // (it needs this to call createImage / loadImage, which both reserve GPU memory)
  sio = new SampleUtilsIO(this);
  v = sio.loadImage("4x4.png");
}

void draw() {
  background(0);
  fill(255);
  
  // also try WRAP, BLACK, and WHITE
  v.edgeMode = VecImage.EdgeModes.CLAMP;
	
  textAlign(CENTER, CENTER);
  textSize(15);
  v.filterMode = VecImage.FilterModes.NEAREST;
  image(sio.toPImage(Filters.upscale(v, 200, 200)), 0, 0);
  text("Nearest", 100, 215);
  translate(200, 0);
	
  // Commonly used as the default way to zoom in on an image
  v.filterMode = VecImage.FilterModes.BILINEAR;
  image(sio.toPImage(Filters.upscale(v, 200, 200)), 0, 0);
  text("Bilinear", 100, 215);
  translate(200, 0);
  
  // A smoother option that pulls values from adjacent pixels to interpolate
  v.filterMode = VecImage.FilterModes.BICUBIC;
  image(sio.toPImage(Filters.upscale(v, 200, 200)), 0, 0);
  text("Bicubic", 100, 215);
  translate(200, 0);
  
  // A clamped version of the mathematically 'ideal' function (sinc)
  // In practice, Lanczos tends to create undesirable 'ringing' artifacts
  v.filterMode = VecImage.FilterModes.LANCZOS;
  v.LANCZOS_RADIUS = 3; // How many neighbors do we consider?
  image(sio.toPImage(Filters.upscale(v, 200, 200)), 0, 0);
  text("Lanczos (a=3)", 100, 215);
  translate(200, 0);

  // BC splines are a family of bidimensional cubisc splines often used for
  // reconstruction and interpolation in computer grpahics.
  v.filterMode = VecImage.FilterModes.BCSPLINE;
  
  // See the docs for additional modes!
  v.BCSplineMode = VecImage.BCSplineModes.PHOTOSHOP_BICUBIC;
  
  image(sio.toPImage(Filters.upscale(v, 200, 200)), 0, 0);
  text("BC Spline (Photoshop Cubic)", 100, 215);
  translate(200, 0);

  saveFrame("interpolation_example_output.png");
}
