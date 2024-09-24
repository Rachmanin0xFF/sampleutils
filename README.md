**sampleutils** is an image resampling, filtering, and utility library for [Processing](https://processing.org/).

Note: This is a WIP; I'm cobbling together this library from various functions and utility scripts I've written over the years. First release should be out soon.

## Interpolation/Resampling
![interpolation example](examples/sampleutils_interpolation/interpolation_example_output.png)
### Supported:
* Nearest-neighbor interpolation (no interpolation)
* Bilinear interpolation
* Bicubic interpolation
* Mitcheljl-netravali filters (BC-splines)
* Lanczos resampling
* EXX2 upscaling
* Poisson disk generation
* Sunflower disk generation
### Planned:
* Box filtering
* Gaussian fitlering
* Arbitrary kernel filtering
* FFT methods
### *Probably* Not planned
* CNN-based/DLSS-style upscaling

## Filters
### Supported:
* Circular medioid blurring
* Inversion
### Planned:
* Kernel blurs (convolutions)
* Unsharp masking / local contrast enhancement
* Tone mapping
* Gamut Masking / clipping
* Dithered curve adjustment
* Component extraction
* Saturation/luma/chromaticity, etc.
* Color space transformations
* LUT extraction from image pairs

## Quantization/Dithering
### Planned:
* Error diffusion matrix-based methods
* Blue noise
* Bayer matrices
* Feature-aware methods


