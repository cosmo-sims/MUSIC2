MUSIC - multi-scale cosmological initial conditions
===================================================

MUSIC is a computer program to generate nested grid initial conditions for
high-resolution "zoom" cosmological simulations. A detailed description
of the algorithms can be found in [Hahn & Abel (2011)][1]. You can
download the user's guide [here][3]. Please consider joining the
[user mailing list][2].

Current MUSIC key features are:

- Supports output for RAMSES, ENZO, Gadget-2/3, ART, Pkdgrav/Gasoline 
and NyX via pugins. New codes can be added.

- Support for first (1LPT) and second order (2LPT) Lagrangian perturbation 
theory, local Lagrangian approximation (LLA) for baryons with grid codes.

- Pluggable transfer functions, currently CAMB, Eisenstein&Hu, BBKS, Warm 
Dark Matter variants. Distinct baryon+CDM fields.

- Minimum bounding ellipsoid and convex hull shaped high-res regions supported 
with most codes, supports refinement mask generation for RAMSES.

- Parallelized with OpenMP
    
- Requires FFTW (v2 or v3), GSL (and HDF5 for output for some codes)


This program is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
or FITNESS FOR A PARTICULAR PURPOSE. By downloading and using MUSIC, you 
agree to the LICENSE, distributed with the source code in a text 
file of the same name.


[1]: http://arxiv.org/abs/1103.6031
[2]: https://groups.google.com/forum/#!forum/cosmo_music
[3]: https://bitbucket.org/ohahn/music/downloads/MUSIC_Users_Guide.pdf
