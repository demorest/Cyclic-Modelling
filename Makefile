PROGS1 = filter_profile
CC = gcc
CFLAGS = -g -Wall -O3
LDLIBS = -lcfitsio -lfftw3f_threads -lfftw3f -lnlopt -lm
all: $(PROGS1)
clean:
	rm -rf *.o
install: all
	cp -f $(PROGS1) $(LOCAL)/bin;

filter_profile: filter_profile.o merit_functions.o model_cyclic.o \
	cyclic_utils.o cyclic_fileio.o
