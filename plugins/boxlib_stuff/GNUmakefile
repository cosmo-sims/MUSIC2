


# Do not change the following variables here, but in Make.local


# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# it is set via environment Variable otherwise it's set to this:
BOXLIB_HOME = ${HOME}/nyx_tot_sterben/BoxLib

#TOP defines the directory in which we will find Source, Exec, etc
TOP = ${HOME}/music-stuff/music/plugins/boxlib_stuff

# PRECISION  HAS TO BE DOUBLE!!! SINGLE won't work out of the box because of fftw3 issues
PRECISION = DOUBLE
DEBUG	  = FALSE
DIM       = 3
COMP	  = Intel
FCOMP     = Intel
USE_MPI   = FALSE
USE_OMP	  = FALSE
EBASE     = Nyx-ic
USE_NU    = FALSE
USE_NONGAUSS = FALSE

# on-site changes to these variables go to Make.local:
-include Make.local

IN_MUSIC = NO

#include ./Make.package

DEFINES += -DFFTW3


include ./Make.ic
