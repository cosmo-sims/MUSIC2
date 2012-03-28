##############################################################################
### compile time configuration options
FFTW3		= no
MULTITHREADFFTW	= yes
SINGLEPRECISION	= no
HAVEHDF5        = yes

##############################################################################
### compiler and path settings
CC      = g++ 
OPT     = -Wall -O3 -g -msse2
CFLAGS  =  
LFLAGS  = -lgsl -lgslcblas 
CPATHS  = -I. -I$(HOME)/local/fftw-3.2.2_double -I$(HOME)/local/include -I/opt/local/include -I/usr/local/include
LPATHS  = -L$(HOME)/local/fftw-3.2.2_double -L$(HOME)/local/lib -L/opt/local/lib -L/usr/local/lib

##############################################################################
# if you have FFTW 2.1.5 or 3.x with multi-thread support, you can enable the 
# option MULTITHREADFFTW
ifeq ($(MULTITHREADFFTW), yes)
  CFLAGS += -fopenmp
  LFLAGS += -fopenmp
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
all: $(OBJS) $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(LPATHS) -o $@ $^ $(LFLAGS)

%.o: %.cc *.hh Makefile 
	$(CC) $(CFLAGS) $(CPATHS) -c $< -o $@

clean:
	rm -rf $(OBJS)

