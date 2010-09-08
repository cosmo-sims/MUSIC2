##############################################################################
### compile time configuration options
MULTITHREADFFTW	= yes
SINGLEPRECISION	= no
HAVEHDF5        = no

##############################################################################
### compiler and path settings
CC      = g++
OPT     = -O3 -msse2
CFLAGS  = -fopenmp 
LFLAGS  = -fopenmp -lgsl -lgslcblas 
CPATHS  = -I. -I$(HOME)/local/include -I/opt/local/include -I/usr/local/include
LPATHS  = -L$(HOME)/local/lib -L/opt/local/lib -L/usr/local/lib

##############################################################################
# if you have FFTW 2.1.5 with multi-thread support, you can enable the option
ifeq ($(MULTITHREADFFTW), yes)
  ifeq ($(SINGLEPRECISION), yes)
    LFLAGS  += -lsrfftw_threads -lsfftw_threads
  else
    LFLAGS  += -ldrfftw_threads -ldfftw_threads
  endif
else
  CFLAGS  += -DSINGLETHREAD_FFTW
endif

##############################################################################
# this section makes sure that the correct FFTW libraries are linked
ifeq ($(SINGLEPRECISION), yes)
  CFLAGS  += -DSINGLE_PRECISION
  LFLAGS  += -lsrfftw -lsfftw
else
  LFLAGS  += -ldrfftw -ldfftw
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
OBJS    = output.o transfer_function.o Numerics.o \
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

