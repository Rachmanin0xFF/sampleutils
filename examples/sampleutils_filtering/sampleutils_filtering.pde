import sampleutils.*;

PImage q;
SampleUtilsIO sio;
VecImage v;

void setup() {
  size(1600, 900, P2D);
  background(255, 0, 255);
  sio = new SampleUtilsIO(this);
  v = sio.loadImage("4x4.png");
}

void draw() {
  background(0);
  fill(255);
  v.edgeMode = VecImage.EdgeModes.CLAMP;

  textAlign(CENTER, CENTER);
  textSize(16);
  v.filterMode = VecImage.FilterModes.NEAREST;
  image(sio.toPImage(Filters.upscale(v, 200, 200)), 0, 0);
  text("Nearest", 100, 220);
  translate(200, 0);

  v.filterMode = VecImage.FilterModes.BILINEAR;
  image(sio.toPImage(Filters.upscale(v, 200, 200)), 0, 0);
  text("Bilinear", 100, 220);
  translate(200, 0);

  v.filterMode = VecImage.FilterModes.BICUBIC;
  image(sio.toPImage(Filters.upscale(v, 200, 200)), 0, 0);
  text("Bicubic", 100, 220);
  translate(200, 0);

  v.filterMode = VecImage.FilterModes.LANCZOS;
  v.LANCZOS_RADIUS = 3;
  image(sio.toPImage(Filters.upscale(v, 200, 200)), 0, 0);
  text("Lanczos (a=3)", 100, 220);
  translate(200, 0);

  v.filterMode = VecImage.FilterModes.BCSPLINE;
  v.BCSplineMode = VecImage.BCSplineModes.PHOTOSHOP_CUBIC;
  image(sio.toPImage(Filters.upscale(v, 200, 200)), 0, 0);
  text("BC Spline (Photoshop Cubic)", 100, 220);
  translate(200, 0);
}