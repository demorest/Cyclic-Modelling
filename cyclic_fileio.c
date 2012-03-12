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
#include "cyclic_fileio.h"

void cyclic_file_error_check_fatal(struct cyclic_file *cf) {
    if (cf->err_status) {
        fits_report_error(stderr, cf->err_status);
        exit(1);
    }
}

int cyclic_file_open(struct cyclic_file *cf, const char *fname) {
    fitsfile *f;
    fits_open_file(&f, fname, READONLY, &cf->err_status);
    cf->file_ptr = (void *)f;
    return(cf->err_status);
}

/* Load dimension params from fits file								*/
/* MAW added *nspec parameter										*/
int cyclic_load_params(struct cyclic_file *cf, struct cyclic_work *w,
					   int *nspec) {

    fitsfile *f = (fitsfile *)cf->file_ptr;
    int bitpix, naxis; 
    long naxes[4];

    fits_get_img_param(f, 4, &bitpix, &naxis, naxes, &cf->err_status);
    if (naxis!=4) { return(-1); }

    w->nphase = naxes[0];
    w->npol   = naxes[1];
    w->nchan  = naxes[2];
	*nspec    = (int)naxes[3];

    w->nlag   = 0;
    w->nharm  = 0;

    return(cf->err_status);
}

/* Load one periodic spectrum from datafile 
 * Space should already be allocated.
 * idx is 1-offset following cfitsio convention.
 */
int cyclic_load_ps(struct cyclic_file *cf, PS *d, int idx) {

    /* Load data */
    fitsfile *f = (fitsfile *)cf->file_ptr;
    int *status = &cf->err_status;
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
