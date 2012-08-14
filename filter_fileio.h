/* filter_fileio.h
 * simple funcs for reading/writing profile and filter data.
 */

#ifndef _FILTER_FILEIO_H
#define _FILTER_FILEIO_H

#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include "cyclic_utils.h"

/* Output simple text-based versions of various quantities			*/
void write_profile(const char *fname, struct profile_phase *p);
void write_fprofile(const char *fname, struct profile_harm *p);
void write_filter(const char *fname, struct filter_time *h);
void write_filter_freq(const char *fname, struct filter_freq *h);

/* Read in the pulse profile										*/
void read_profile(const char *fname, struct profile_phase *pp);

/* Read in a filter solution file */
int read_filters(const char *fname, int nspec_expect, int nchan_expect,
        fftwf_complex **data);

/* Init a FITS-format filter output file */
int init_filter_file(const char *fname, int nspec, int nchan, int nbin,
        int imjd0, double fmjd0);
/* Add a filter to the FITS file */
int append_filter(const char *fname, int ispec, int imjd, double fmjd, 
        struct filter_freq *hf);
/* Add a profile to the FITS file */
int append_profile(const char *fname, int ispec, struct profile_phase *pp);

#endif
