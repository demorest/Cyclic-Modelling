/* cyclic_fileio.h
 *
 * Structs for cyclic spectra I/O
 */
#ifndef _CYCLIC_FILIO_H
#define _CYCLIC_FILIO_H

#include "cyclic_utils.h"

/* Structure that will allow us to abstract out the 
 * actual file reading implementation.  We will use
 * a error status variable to indicate error conds.
 * For the fitsio-based routines this is just the
 * usual status var.
 */
struct cyclic_file {
    void *file_ptr; /* pointer to the actual file object */
    int err_status; /* 0=no error, non-zero=error */
};

/* Open the filename, and put fill results into cf struct */
int cyclic_file_open(struct cyclic_file *cf, const char *fname);

/* Check the error condition and exit with a message if non-zero */
void cyclic_file_error_check_fatal(struct cyclic_file *cf);

/* Load dimension params from datafile								*/
/* MAW modified to include *nspec									*/
int cyclic_load_params(struct cyclic_file *cf, struct cyclic_work *w,
					   int *nspec);

/* Load one periodic spectrum from datafile							*/
int cyclic_load_ps(struct cyclic_file *cf, PS *d, int idx);

#endif

