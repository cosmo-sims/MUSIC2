MUSIC - multi-scale cosmological initial conditions
===================================================

MUSIC is a computer program to generate nested grid initial conditions for
high-resolution "zoom" cosmological simulations. A detailed description
of the algorithms can be found in [Hahn & Abel (2011)][1]. You can
download the user's guide [here][3], or [read the Wiki](https://bitbucket.org/ohahn/music/wiki/Home) instead. Please consider joining the
[user mailing list][2].

Current MUSIC key features are:

- Supports output for RAMSES, ENZO, Arepo, Gadget-2/3, ART, Pkdgrav/Gasoline 
and NyX via plugins. New codes can be added.

- Support for first (1LPT) and second order (2LPT) Lagrangian perturbation 
theory, local Lagrangian approximation (LLA) for baryons with grid codes.

- Pluggable transfer functions, currently CAMB, Eisenstein&Hu, BBKS, Warm 
Dark Matter variants. Distinct baryon+CDM fields.

- Minimum bounding ellipsoid and convex hull shaped high-res regions supported 
with most codes, supports refinement mask generation for RAMSES.

- Parallelized with OpenMP
    
- Requires FFTW (v2 or v3), GSL (and HDF5 for output for some codes)

## Building MUSIC
While we still supply the old Makefile, using CMake is now the preferred way of building. 
CMake use out-of-source build, i.e. you create a build directory, and then configure the code using CMake. Inside the `music` directory, do
```
  mkdir build
  cd build
  ccmake ..
  make -j
```
to configure the code (you will se a menu), and then start a parallel compilation. If CMake has trouble finding your FFTW or HDF5 installation,
you can add hints as follows
```
  FFTW3_ROOT=<path> HDF5_ROOT=<path> ccmake ..
```
If you want to build on macOS, then it is strongly recommended to use GNU (or Intel) compilers instead of Apple's Clang. Install them e.g. via homebrew and then configure cmake to use them instead of the macOS default compiler via
```
  CC=gcc-9 CXX=g++-9 ccmake ..
```
This is necessary since Apple's compilers haven't supported OpenMP for years.


## Running

There is an example parameter file 'example.conf' in the main directory. Possible options are explained in it, it can be run
as a simple argument, e.g. from within the build directory:

```
  ./MUSIC ../ics_example.conf
```


## Disclaimer

This program is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
or FITNESS FOR A PARTICULAR PURPOSE. By downloading and using MUSIC, you 
agree to the LICENSE, distributed with the source code in a text 
file of the same name.

## References
[1]: http://arxiv.org/abs/1103.6031
[2]: https://groups.google.com/forum/#!forum/cosmo_music
[3]: https://bitbucket.org/ohahn/music/downloads/MUSIC_Users_Guide.pdf
