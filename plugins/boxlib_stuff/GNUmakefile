# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# it is set via environment Variable otherwise it's set to this:
#BOXLIB_HOME = ${HOME}/nyx_tot_sterben/BoxLib

#TOP defines the directory in which we will find Source, Exec, etc
TOP = ${PWD}

# PRECISION  HAS TO BE DOUBLE!!! SINGLE won't work out of the box because of fftw3 issues
PRECISION = DOUBLE
DEBUG	  = FALSE
DIM       = 3
COMP	  = Intel
FCOMP     = Intel
USE_MPI   = FALSE
USE_OMP	  = FALSE

IN_MUSIC = NO

#hmm....awful (needed for general.hh)
ifeq ($(FFTW3), yes)
  DEFINES += -DFFTW3
endif

ifeq ($(SINGLE), yes)
  PRECISION = SINGLE
endif


include ./Make.ic
