##############################################################################
### compile time configuration options
FFTW3		= yes
MULTITHREADFFTW	= yes
MULTITHREADFFTW	= no
MULTITHREADFFTW	= no
SINGLEPRECISION	= no
HAVEHDF5        = no
HAVEBOXLIB	= yes

##############################################################################
### compiler and path settings
CC      = mpiicpc
OPT     = -Wall -O3 -g -msse2
CFLAGS  =  
LFLAGS  = -lgsl -lgslcblas 
CPATHS  = -I. -I$(HOME)/local/include -I/opt/local/include -I/usr/local/include
LPATHS  = -L$(HOME)/local/lib -L/opt/local/lib -L/usr/local/lib

CPATHS  = -I. -I$(HOME)/music-stuff/inst/include 
LPATHS  = -L$(HOME)/music-stuff/inst/lib 

##############################################################################
# if you have FFTW 2.1.5 or 3.x with multi-thread support, you can enable the 
# option MULTITHREADFFTW
ifeq ($(MULTITHREADFFTW), yes)
  ifeq ($(CC), mpiicpc)
    CFLAGS += -openmp
    LFLAGS += -openmp
  else
    CFLAGS += -fopenmp
    LFLAGS += -fopenmp
  endif
  ifeq ($(FFTW3),yes)
	ifeq ($(SINGLEPRECISION), yes)
		LFLAGS  +=  -lfftw3f_threads
	else
		LFLAGS  +=  -lfftw3_threads
	endif
  else
    ifeq ($(SINGLEPRECISION), yes)
      LFLAGS  += -lsrfftw_threads -lsfftw_threads
    else
      LFLAGS  += -ldrfftw_threads -ldfftw_threads
    endif
  endif
else
  CFLAGS  += -DSINGLETHREAD_FFTW
endif

ifeq ($(FFTW3),yes)
  CFLAGS += -DFFTW3
endif

##############################################################################
# this section makes sure that the correct FFTW libraries are linked
ifeq ($(SINGLEPRECISION), yes)
  CFLAGS  += -DSINGLE_PRECISION
  ifeq ($(FFTW3),yes)
    LFLAGS += -lfftw3f
  else
    LFLAGS  += -lsrfftw -lsfftw
  endif
else
  ifeq ($(FFTW3),yes)
    LFLAGS += -lfftw3
  else
    LFLAGS  += -ldrfftw -ldfftw
  endif
endif

##############################################################################
#if you have HDF5 installed, you can also enable the following options
ifeq ($(HAVEHDF5), yes)
  OPT += -DH5_USE_16_API -DHAVE_HDF5
  LFLAGS += -lhdf5
endif

##############################################################################
CFLAGS += $(OPT)
TARGET  = MUSIC
OBJS    = output.o transfer_function.o Numerics.o defaults.o constraints.o random.o\
		convolution_kernel.o densities.o cosmology.o poisson.o log.o main.o \
		$(patsubst plugins/%.cc,plugins/%.o,$(wildcard plugins/*.cc))

##############################################################################
# stuff for BoxLib
BLOBJS = ""
ifeq ($(HAVEBOXLIB), yes)
  IN_MUSIC = YES
  BOXLIB_HOME ?= ${HOME}/nyx_tot_sterben/BoxLib
  TOP = ${HOME}/music-stuff/music/plugins/boxlib_stuff
  CCbla := $(CC)
  include plugins/boxlib_stuff/Make.ic
  CC  := $(CCbla)
  CPATHS += $(INCLUDE_LOCATIONS)
  LPATHS += -L$(objEXETempDir)
  BLOBJS = $(foreach obj,$(objForExecs),plugins/boxlib_stuff/$(obj))
#
endif

##############################################################################
all: $(OBJS) $(TARGET)
	cd plugins/boxlib_stuff; make

bla:
	echo $(BLOBJS)

#FIXME!!!
$(TARGET): $(OBJS) 
	cd plugins/boxlib_stuff; make
	$(CC) $(LPATHS) -o $@ $^ $(LFLAGS) $(BLOBJS) -lifcore

#%.o: %.cc *.hh Makefile 
%.o: %.cc *.hh
	$(CC) $(CFLAGS) $(CPATHS) -c $< -o $@

clean:
	rm -rf $(OBJS)

