PROGS1 = filter_profile
CC = gcc

# Uncomment these to use psrchive for reading input
PSRCHIVE_CFLAGS = $(shell psrchive --cflags)
PSRCHIVE_LIBS = $(shell psrchive --libs)
FILEIO_OBJ = cyclic_fileio_psrchive.o

# This didn't used to be necessary...
PSRCHIVE_LIBS += -lstdc++

# Uncomment these to use old fits file format
#PSRCHIVE_CFLAGS = 
#PSRCHIVE_LIBS = 
#FILEIO_OBJ = cyclic_fileio.o

CFLAGS = -g -Wall -O3
CXXFLAGS = $(PSRCHIVE_CFLAGS)
LDLIBS = -lcfitsio -lfftw3f_threads -lfftw3f -lnlopt -lm $(PSRCHIVE_LIBS)
all: $(PROGS1)
clean:
	rm -rf *.o
install: all
	cp -f $(PROGS1) $(LOCAL)/bin;

filter_profile: filter_profile.o merit_functions.o model_cyclic.o \
	filter_fileio.o cyclic_utils.o $(FILEIO_OBJ)
