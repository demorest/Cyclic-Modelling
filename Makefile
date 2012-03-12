PROGS1 = filter_profile
CC = gcc

PSRCHIVE_CFLAGS = $(shell psrchive --cflags)
PSRCHIVE_LIBS = $(shell psrchive --libs)

CFLAGS = -g -Wall -O3
CXXFLAGS = $(PSRCHIVE_CFLAGS)
LDLIBS = -lcfitsio -lfftw3f_threads -lfftw3f -lnlopt -lm $(PSRCHIVE_LIBS)
all: $(PROGS1)
clean:
	rm -rf *.o
install: all
	cp -f $(PROGS1) $(LOCAL)/bin;

#filter_profile: filter_profile.o merit_functions.o model_cyclic.o \
#	cyclic_utils.o cyclic_fileio.o
filter_profile: filter_profile.o merit_functions.o model_cyclic.o \
	cyclic_utils.o cyclic_fileio_psrchive.o
