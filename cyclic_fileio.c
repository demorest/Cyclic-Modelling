/* cyclic_fileio.c
 * P. Demorest, March 2012.
 * File I/O routines for 'old-style' cyclic spectra FITS
 * files.  Split code out of cyclic_utils.c.
 */
#include <fitsio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include "cyclic_utils.h"

/* Load dimension params from fits file								*/
/* MAW added *nspec parameter										*/
int cyclic_load_params(fitsfile *f, struct cyclic_work *w,
					   int *nspec, int *status) {

    int bitpix, naxis; 
    long naxes[4];

    fits_get_img_param(f, 4, &bitpix, &naxis, naxes, status);
    if (naxis!=4) { return(-1); }

    w->nphase = naxes[0];
    w->npol   = naxes[1];
    w->nchan  = naxes[2];
	*nspec    = (int)naxes[3];

    w->nlag   = 0;
    w->nharm  = 0;

    return(*status);
}

/* Load one periodic spectrum from datafile 
 * Space should already be allocated.
 * idx is 1-offset following cfitsio convention.
 */
int cyclic_load_ps(fitsfile *f, PS *d, int idx, int *status) {

    /* Load data */
    long fpixel[4];
    long long nelem;
    fpixel[0] = fpixel[1] = fpixel[2] = 1;
    fpixel[3] = idx;
    nelem = d->nphase * d->npol * d->nchan;
    fits_read_pix(f, TFLOAT, fpixel, nelem, NULL, d->data, NULL, status);

    /* Load header params */
    char key[9];
    sprintf(key, "IMJD%04d", idx-1);
    fits_read_key(f, TINT, key, &d->imjd, NULL, status);
    sprintf(key, "FMJD%04d", idx-1);
    fits_read_key(f, TDOUBLE, key, &d->fmjd, NULL, status);
    sprintf(key, "PHAS%04d", idx-1);
    fits_read_key(f, TDOUBLE, key, &d->ref_phase, NULL, status);
    sprintf(key, "FREQ%04d", idx-1);
    fits_read_key(f, TDOUBLE, key, &d->ref_freq, NULL, status);
    // TODO get these in the file
    d->rf = 428.0;
    d->bw = 4.0;

    return(*status);
}
