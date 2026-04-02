![Build Status](https://github.com/OpenPIV/openpiv-c--qt/actions/workflows/build_test.yml/badge.svg)

# OpenPIV (c++), a fast, open-source particle image velocimetry (PIV) library

## An implementation of a PIV analysis engine in C++ using as few dependencies as possible; the implementation requires a c++17 compliant compiler.

This project is the result of the collaborative effort of numerous researchers in order to provide one of the fastest PIV software on the market while remaining cross-platform and open-source. The software can do the following:

 * Load images
 * Save images
 * Pre-process and modify images
 * Perform digital PIV analysis including subpixel estimation
 * Augment into Python workflows
 * and more!
 
## Image Loaders

Loading and storing images are crucial for any PIV software. Due to this requirement, openpiv-c--qt implements image loaders that can load, convert, and store images.
Currently, there are a few extensions that are supported, but more are under development.

| Supported Extensions     | Decode | Encode |
|-------------------------|--------|--------|
| .b16 (PCO CamWare :tm: )| Planned| -      | 
| .bmp                    | Planned| Planned| 
| .jpeg                   | Planned| Planned| 
| .png                    | Planned| Planned| 
| .pnm (.pbm, .pgm, .ppm) | x      | x      |
| .tif                    | x      | x      |
| .webp                   | Planned| Planned|                    

## Python Interface
The core functions of the OpenPIV c++ library are wrapped in helper functions to allow for efficient use in a Python environment. This enables for the efficiency and performance of c++ to be used in the ease of use and flexibility of Python. Users can expect significant performance increased compared to Python/NumPy based implementations. 

> [!NOTE]
> Image loaders should be disabled when the Python wrapper is built to enable better portability to cloud-based resources.

## Examples

To demonstrate how this library can perform, two examples are provided in the /examples folder. The first example,
[average_subtract](examples/average_subtract/README.md), is a utility that reads n images, calculates the average, and
writes out n new images with the average subtracted. Additionally, a second example, [process](examples/process/README.md),
is a straight-forward PIV cross correlator that reads two images and performs cross-correlation on them.
    
## Build

There are some external dependencies under external/, so when cloning use:

```git clone --recursive <path to git repo>```

Building uses meson, and is simplified by using meson wrap files to specify the dependent packages. Building has some pre-requisites:

* a compiler (e.g. `apt install build-essentials`)
* cmake (optional)
* python
* (linux) pkg-config (e.g. `apt install pkg-config`)
* curl, zip, unzip, tar (e.g. `apt install curl zip unzip tar`)
* ninja (e.g. `apt install ninja-build`)
* meson (e.g., `pip install --user meson`)

Unix users can also use a similar method to the Windows build environment as detailed below.

On Windows, the following can be used:
* install Visual Studio 2019 or 2022. Alternatively, MinGW-64 and Intel OneAPI c/c++ compilers are also known to work
* install miniconda or python along with venv and setup virtual environment
* pip install cmake
* pip install ninja
* pip install meson

> [!NOTE]
> Conda environments allows for GNU GCC toolchains to be installed in a pain-free fashion.

To build:
* `meson setup builddir`
* `meson compile -C builddir`
> [!NOTE]
> It is good practice to setup `--prefix` flags so files are not installed on the current directory.
>
> SIMD can be enabled using GCC and passing `-march=native` like this: `meson setup builddir -Dcpp_args="-march=native"`
>
> The Python wrapper can be built by enabling the build_wrapper option like this: `meson setup builddir -Dbuild_wrapper=true -Ddisable_image_loaders=true`

Meson provides multiple build types such as debug, debugoptimized, and release. To change the build type, use the `--buildtype` flag. For example, `meson setup builddir --buildtype debugoptimized`.

To run tests:

* `meson test -C builddir'

To get binaries:
* `meson install -C builddir` if the prefix was set or
* `meson install -C builddir --destdir <some directory>` to install in a specific directory.

Sometimes you only want the runtime dynamic libraries and executables. Meson comes with a handy targeted installation using the following command:
 * `meson install -C builddir --tags runtime`

Make sure the prefix, or destdir, is set so binaries are not accidentally installed somewhere unexpected.

The binaries are located in the build or installation directory:

Build directory:
* builddir
  * examples
    * process
    * average_subtract
  * openpiv -> libopenpivcore.so

Install directory:
* prefix/destdir
  * bindir
    * libopenpivcore.so
    * all other dependent shared libraries
    * process (executable)
    * average_subtract (executable)
    
> [!WARNING]
> When using GCC on Windows, libstdc++-6.dll and libgcc_s_seh-1.dll must be in PATH or next to the dll libraries in order to avoid a missing dll error. When using a different environment from where the build took place, these dll libraries would usually need to be copy and pasted next to libopenpivcore.

## Dependencies
* c++17 compiler e.g. clang++-5.0, gcc8
* python3
* [meson](https://mesonbuild.com/index.html)
  * benchmark: used to run performance benchmarks
  * catch2: unit test framework
  * cxxopts: nice command line parsing
  * libtiff: TIFF IO support
    * libjpeg-turbo
    * zlib

## Performance

### Raspberry Pi 3B
```
$ lscpu
Architecture:                    aarch64
CPU op-mode(s):                  32-bit, 64-bit
Byte Order:                      Little Endian
CPU(s):                          4
On-line CPU(s) list:             0-3
Thread(s) per core:              1
Core(s) per socket:              4
Socket(s):                       1
Vendor ID:                       ARM
Model:                           4
Model name:                      Cortex-A53
Stepping:                        r0p4
CPU max MHz:                     1200.0000
CPU min MHz:                     1200.0000
BogoMIPS:                        38.40
Vulnerability Itlb multihit:     Not affected
Vulnerability L1tf:              Not affected
Vulnerability Mds:               Not affected
Vulnerability Meltdown:          Not affected
Vulnerability Spec store bypass: Not affected
Vulnerability Spectre v1:        Mitigation; __user pointer sanitization
Vulnerability Spectre v2:        Not affected
Vulnerability Srbds:             Not affected
Vulnerability Tsx async abort:   Not affected
Flags:                           fp asimd evtstrm crc32 cpuid

$ time ./process -s 32 -i ../../../examples/images/F_00001.tif ../../../examples/images/F_00002.tif > out.piv
[15:59:23 +00:00] [info] [thread 40843] size: 32 x 32
[15:59:23 +00:00] [info] [thread 40843] overlap: 0.5
[15:59:23 +00:00] [info] [thread 40843] input files: [../../../examples/images/F_00001.tif, ../../../examples/images/F_00002.tif]
[15:59:23 +00:00] [info] [thread 40843] execution: pool
[15:59:23 +00:00] [info] [thread 40843] loaded images: [1024,1024]
[15:59:23 +00:00] [info] [thread 40843] generated grid for image size: [1024,1024], ia: [32,32] (50% overlap)
[15:59:23 +00:00] [info] [thread 40843] grid count: 3969
[15:59:23 +00:00] [info] [thread 40843] processing using thread pool

real    0m1.248s
user    0m3.249s
sys     0m0.109s
```
This is about 1ms per interrogation area (3 cores, 3969 interrogation areas, 1.248s) without using SIMD optimizations.

### AMD Ryzen 7 Laptop (Windows 11)
```
> wmic cpu get name,NumberOfCores,MaxClockSpeed,NumberOfLogicalProcessors,L2CacheSize,L3CacheSize
L2CacheSize  L3CacheSize  MaxClockSpeed  Name                                    NumberOfCores  NumberOfLogicalProcessors
4096         16384        2000           AMD Ryzen 7 7730U with Radeon Graphics  8              16
```
MSVC 2022 compiler (no SIMD optimizations)
```
> Measure-Command { ./process -s 32 -i exp1_001_a.tiff exp1_001_b.tiff -f pocket > out.piv }
[404804396117100] (24040) INFO: size: 32 x 32
[404804396118200] (24040) INFO: overlap: 0.5
[404804396265000] (24040) INFO: input files: [exp1_001_a.tiff, exp1_001_b.tiff]
[404804396266800] (24040) INFO: execution: pool
[404804409030500] (24040) INFO: loaded images have size: [511,369]
[404804409090600] (24040) INFO: generated grid for image size: [511,369], ia: [32,32] (50% overlap)
[404804409091800] (24040) INFO: grid count: 660
[404804409123400] (24040) INFO: processing using thread pool
**[404804422928400] (24040) INFO: processing time: 13804.900000us, 20.916515us per interrogation area


Days              : 0
Hours             : 0
Minutes           : 0
Seconds           : 0
Milliseconds      : 133
Ticks             : 1338931
TotalDays         : 1.54968865740741E-06
TotalHours        : 3.71925277777778E-05
TotalMinutes      : 0.00223155166666667
TotalSeconds      : 0.1338931
TotalMilliseconds : 133.8931

```
MingW-64 compiler (SIMD optimizations)
```
> Measure-Command { ./process -s 32 -i exp1_001_a.tiff exp1_001_b.tiff -f pocket > out.piv }
[1774808548511394000] (0xb8) INFO: size: 32 x 32
[1774808548511398000] (0xb8) INFO: overlap: 0.5
[1774808548511425000] (0xb8) INFO: input files: [exp1_001_a.tiff, exp1_001_b.tiff]
[1774808548511431000] (0xb8) INFO: execution: pool
[1774808548523492000] (0xb8) INFO: loaded images have size: [511,369]
[1774808548523524000] (0xb8) INFO: generated grid for image size: [511,369], ia: [32,32] (50% overlap)
[1774808548523525000] (0xb8) INFO: grid count: 660
[1774808548523545000] (0xb8) INFO: processing using thread pool
**[1774808548531184000] (0xb8) INFO: processing time: 7643.000000us, 11.580303us per interrogation area


Days              : 0
Hours             : 0
Minutes           : 0
Seconds           : 0
Milliseconds      : 88
Ticks             : 889349
TotalDays         : 1.02933912037037E-06
TotalHours        : 2.47041388888889E-05
TotalMinutes      : 0.00148224833333333
TotalSeconds      : 0.0889349
TotalMilliseconds : 88.9349

```
As you can see, GCC with SIMD enabled can have very significant performance increases compared to using MSVC compilers.

# TODO

* build
  * [x] travis/github actions/CI
  * [x] add clang/windows/OS X build
* core
  * [x] logging
  * [ ] iostream ops for ImageLoaders
  * [ ] image allocator support
    * [ ] shared memory
    * [ ] pooled memory
  * [ ] image loading
    * [x] load/save PNM files
      * [x] save double image
    * [x] support loading multiple images from files: modify loader interface to read header, report number of images, extract indexed images
    * [ ] memory map files - check performance for large files
    * [ ] PNG - lodepng
    * [ ] RAW - libraw looks less than ideal but no alternative?
    * [ ] BMP - easyBMP?
    * [ ] JPEG - already incl in libjpeg-turbo
    * [ ] b16/PCO
  * utilities
    * [x] split RGB to 4xgreyscale
    * [x] join 4xgreyscale to RGB
    * [x] split complex to planes
    * [x] create complex from planes
  * [x] add ImageInterface data iterators
  * [x] remove data() method from ImageInterface, ImageView
  * [x] allow ImageView to be non-const
* algorithms
  * [x] transpose
  * [x] normalize
  * FFT
    * [x] spectral correlation
    * [x] optimize
    * [ ] openCL
    * [ ] apply kernel in Fourier space
    * [x] use SIMD?
    * [x] real -> complex FFT/correlation of real data
    * [ ] normalized minimum quadratic differences (very robust)
    * [ ] zero pad images
  * direct correlation
    * [ ] full window correlation
    * [ ] partial window correlation (for enhancing FFT correlations)
    * [ ] use SIMD?
  * image deformation
    * [ ] shifted linear image deformation interpolation (same as first degree polynomial for my impl)
    * [ ] polynomial interpolation (lookup tables-based, so super fast!!)
    * [ ] sinc (lookup table-based, so also quite fast; only 7x7 and 11x11 kernels supported)
    * [ ] lanczos (lookup table-based, so also quite fast; generally more stable than sinc)
  * [x] peak detection
  * peak fitting
    * [x] 3 point Gaussian peak fit
    * [ ] 3 point parabolic peak fit
    * [ ] 3 point centroid peak fit
    * [ ] 3x3 least squares Gaussian peak fit (optimized via pseudo-inverse)
    * [ ] nxn least squares Gaussian peak fit (optimized via pseudo-inverse)
    * [ ] nxn non-linear Gaussian peak fit (optimized via Levenberg-Marquardt or something similar)
    * [ ] nxn centroid peak fit (can be used for calibration marker detection)
  * multi-pass PIV
    * [x] First pass w/ multi-threading
    * [ ] Multi-pass image deform w/ multi-threading
    * [ ] Failed correlations can use larger correlation window or imputation
* image processing
  * filters
    * [ ] change image_view to use array of pointers for each row?
    * [ ] Gaussian low-pass filter
    * [ ] Gaussian high-pass filter
    * [ ] normalized variance filter
    * [ ] contrast stretch filter (based on local percentile kernels)
    * [ ] median filter
    * [ ] min/max normalization filter
    * [ ] smoothn from OpenPIV-Python (requires L-BFGS-B impl)
  * adjustment
    * [ ] crop (region of interest/ROI)
    * [ ] skew
    * [ ] stretch
    * [ ] translate
    * [ ] rotate
    * [ ] flip (x or y axis)
    * [ ] normalize (0..1 * scaling value)
    * [ ] deform (Using same deformation algos in PIV image deformation)
  * masking
    * [ ] read mask image
    * [ ] write mask image
    * [ ] automatically create mask image
* examples
  * [x] parallel cross-correlate
  * [x] image processing
* processing framework
  * [x] cartesian grid generator
  * [ ] further grid generators
  * [ ] median validation with secondary peak check and interpolation
  * [ ] store signal/noise value
  * [ ] processing
  * [ ] marking
  * [ ] iterative analysis
  * [ ] PIV guided PTV?
* data output
  * [ ] output registry
  * [ ] ASCII/CSV
  * [ ] gnuplot/pyplot?
  * [ ] tecplot
* interfacing
  * [ ] nanobind (faster pybind11) for Python interfacing
