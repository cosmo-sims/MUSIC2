##############################################################################
### compile time configuration options
FFTW3		= yes
MULTITHREADFFTW	= yes
SINGLEPRECISION	= no
<<<<<<< local
HAVEHDF5        = yes
=======
HAVEHDF5        = no
>>>>>>> other
HAVEBOXLIB	= no
BOXLIB_HOME     = ${HOME}/nyx_tot_sterben/BoxLib

##############################################################################
### compiler and path settings
CC      = g++
<<<<<<< local
OPT     = -Wall -Wno-unknown-pragmas -O3 -g -mtune=native
=======
OPT     = -Wall -Wno-unknown-pragmas -O3 -g -msse2
>>>>>>> other
CFLAGS  =  
LFLAGS  = -lgsl -lgslcblas 
<<<<<<< local
CPATHS  = -I. -I/opt/local/include -I/usr/local/include
LPATHS  = -L/opt/local/lib -L/usr/local/lib
=======
CPATHS  = -I. -I$(HOME)/local/include -I/opt/local/include -I/usr/local/include
LPATHS  = -L$(HOME)/local/lib -L/opt/local/lib -L/usr/local/lib
>>>>>>> other

##############################################################################
# if you have FFTW 2.1.5 or 3.x with multi-thread support, you can enable the 
# option MULTITHREADFFTW
ifeq ($(strip $(MULTITHREADFFTW)), yes)
  ifeq ($(CC), mpiicpc)
    CFLAGS += -openmp
    LFLAGS += -openmp
  else
    CFLAGS += -fopenmp
    LFLAGS += -fopenmp
  endif
  ifeq ($(strip $(FFTW3)),yes)
	ifeq ($(strip $(SINGLEPRECISION)), yes)
		LFLAGS  +=  -lfftw3f_threads
	else
		LFLAGS  +=  -lfftw3_threads
	endif
  else
    ifeq ($(strip $(SINGLEPRECISION)), yes)
      LFLAGS  += -lsrfftw_threads -lsfftw_threads
    else
      LFLAGS  += -ldrfftw_threads -ldfftw_threads
    endif
  endif
else
  CFLAGS  += -DSINGLETHREAD_FFTW
endif

ifeq ($(strip $(FFTW3)),yes)
  CFLAGS += -DFFTW3
endif

##############################################################################
# this section makes sure that the correct FFTW libraries are linked
ifeq ($(strip $(SINGLEPRECISION)), yes)
  CFLAGS  += -DSINGLE_PRECISION
  ifeq ($(FFTW3),yes)
    LFLAGS += -lfftw3f
  else
    LFLAGS  += -lsrfftw -lsfftw
  endif
else
  ifeq ($(strip $(FFTW3)),yes)
    LFLAGS += -lfftw3
  else
    LFLAGS  += -ldrfftw -ldfftw
  endif
endif

##############################################################################
#if you have HDF5 installed, you can also enable the following options
ifeq ($(strip $(HAVEHDF5)), yes)
  OPT += -DH5_USE_16_API -DHAVE_HDF5
  LFLAGS += -lhdf5
endif

##############################################################################
CFLAGS += $(OPT)
TARGET  = MUSIC
OBJS    = output.o transfer_function.o Numerics.o defaults.o constraints.o random.o\
		convolution_kernel.o region_generator.o densities.o cosmology.o poisson.o\
		densities.o cosmology.o poisson.o log.o main.o \
		$(patsubst plugins/%.cc,plugins/%.o,$(wildcard plugins/*.cc))

##############################################################################
# stuff for BoxLib
BLOBJS = ""
ifeq ($(strip $(HAVEBOXLIB)), yes)
  IN_MUSIC = YES
  TOP = ${PWD}/plugins/nyx_plugin
  CCbla := $(CC)
  include plugins/nyx_plugin/Make.ic
  CC  := $(CCbla)
  CPATHS += $(INCLUDE_LOCATIONS)
  LPATHS += -L$(objEXETempDir)
  BLOBJS = $(foreach obj,$(objForExecs),plugins/boxlib_stuff/$(obj))
#
endif

##############################################################################
all: $(OBJS) $(TARGET) Makefile
#	cd plugins/boxlib_stuff; make

bla:
	echo $(BLOBJS)

ifeq ($(strip $(HAVEBOXLIB)), yes)
$(TARGET): $(OBJS) plugins/nyx_plugin/*.cpp
	cd plugins/nyx_plugin; make BOXLIB_HOME=$(BOXLIB_HOME) FFTW3=$(FFTW3) SINGLE=$(SINGLEPRECISION)
	$(CC) $(LPATHS) -o $@ $^ $(LFLAGS) $(BLOBJS) -lifcore
else
$(TARGET): $(OBJS)
	$(CC) $(LPATHS) -o $@ $^ $(LFLAGS)
endif

#%.o: %.cc *.hh Makefile 
%.o: %.cc *.hh
	$(CC) $(CFLAGS) $(CPATHS) -c $< -o $@

clean:
	rm -rf $(OBJS)
ifeq ($(strip $(HAVEBOXLIB)), yes)
	oldpath=`pwd`
	cd plugins/nyx_plugin; make realclean BOXLIB_HOME=$(BOXLIB_HOME)
endif
	cd $(oldpath)
	
