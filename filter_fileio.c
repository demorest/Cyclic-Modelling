/* filter_fileio.c
 *
 * Routines to read/write filter funcs and profiles.
 * Split off from cyclic_utils.c, PBD 2012/08/14.
 */
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <fitsio.h>

#include "cyclic_utils.h"
#include "filter_fileio.h"

void write_profile(const char *fname, struct profile_phase *p) {
    FILE *f = fopen(fname, "w");
    int i;
    for (i=0; i<p->nphase; i++) {
        fprintf(f,"%.7e %.7e\n",
                (double)i/(double)p->nphase, p->data[i]);
    }
    fclose(f);
}

void write_fprofile(const char *fname, struct profile_harm *p) {
    FILE *f = fopen(fname, "a");
    int i;
    for (i=0; i<p->nharm; i++) {
        fprintf(f,"%d %.7e %.7e\n",
                i, creal(p->data[i]), cimag(p->data[i]));
    }
    fprintf(f,"\n\n");
    fclose(f);
}

void write_filter(const char *fname, struct filter_time *h) {
    FILE *f = fopen(fname, "a");
    int i;
    for (i=0; i<h->nlag; i++) {
        fprintf(f,"%d %.7e %.7e\n",
                i, creal(h->data[i]), cimag(h->data[i]));
    }
    fprintf(f,"\n\n");
    fclose(f);
}

void write_filter_freq(const char *fname, struct filter_freq *h) {
    FILE *f = fopen(fname, "a");
    int i;
    for (i=0; i<h->nchan; i++) {
        fprintf(f,"%d %.7e %.7e\n",
                i, creal(h->data[i]), cimag(h->data[i]));
    }
    fprintf(f,"\n\n");
    fclose(f);
}

void read_profile(const char *fname, struct profile_phase *pp) {
    FILE *f = fopen(fname, "r");
    int i;
    float ptmp, dtmp;
    for (i=0; i<pp->nphase; i++) {
        fscanf(f,"%f %f", &ptmp, &dtmp);
        pp->data[i] = dtmp;
    }
    fclose(f);
}

int read_filters(const char *fname, int nspec_expect, int nchan_expect,
        fftwf_complex **data) {

    FILE *fpointer = fopen(fname, "r");
    if (fpointer==NULL) {
        fprintf(stderr,"Cannot open filter file '%s'\n", fname);
        return(-1);
    }

    int nspec, nchan;
    fscanf(fpointer,"%d", &nspec);
    fscanf(fpointer,"%d", &nchan);
    if (nchan!=nchan_expect || nspec!=nspec_expect) {
        fprintf(stderr,"File dimensions (%d,%d) don't match expected (%d,%d)\n",
                nspec, nchan, nspec_expect, nchan_expect);
        fclose(fpointer);
        return(-2);
    }

    int is, ic, rv;
    float rprev=0, iprev=0;
    for (is=0; is<nspec; is++) {
        for (ic=0; ic<nchan; ic++) {
            rv = fscanf(fpointer, "%f %f", &rprev, &iprev);
            if (rv != 2) {
                fprintf(stderr, "Not enough data in filter file\n");
                fclose(fpointer);
                return(-3);
            }
            data[is][ic] = rprev + I * iprev;
        }
    }

    fclose(fpointer);
    return(0);

}

int init_filter_file(const char *fname, int nspec, int nchan, int nphase,
        int imjd0, double fmjd0) {
    int status=0;
    fitsfile *fptr;
    fits_create_file(&fptr, fname, &status);
    fits_movabs_hdu(fptr, 1, NULL, &status);
    long naxes[3] = {2, nchan, nspec};
    fits_create_img(fptr, FLOAT_IMG, 3, naxes, &status);
    fits_write_key(fptr, TINT, "IMJD", &imjd0, NULL, &status);
    fits_write_key(fptr, TDOUBLE, "FMJD", &fmjd0, NULL, &status);
    naxes[0] = nphase;
    naxes[1] = nspec;
    fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status);
    fits_movabs_hdu(fptr, 2, NULL, &status);
    fits_write_key(fptr, TSTRING, "EXTNAME", "PROFILE", NULL, &status);
    fits_close_file(fptr, &status);
    if (status) {
        fprintf(stderr, "Error in init_filter_file:\n");
        fits_report_error(stderr, status);
        return(-1);
    }
    return(0);
}

int append_filter(const char *fname, int ispec, int imjd, double fmjd, 
        struct filter_freq *hf) {
    int status=0;
    fitsfile *fptr;
    fits_open_file(&fptr, fname, READWRITE, &status);
    fits_movabs_hdu(fptr, 1, NULL, &status);
    int naxis;
    long naxes[3];
    fits_get_img_dim(fptr, &naxis, &status);
    if (naxis!=3) { 
        fits_close_file(fptr, &status);
        return(-1); 
    }
    fits_get_img_size(fptr, 3, naxes, &status);
    if (naxes[0]!=2 || naxes[1]!=hf->nchan || naxes[2]<ispec) {
        fits_close_file(fptr, &status);
        return(-2);
    }
    int imjd0;
    double fmjd0;
    fits_read_key(fptr, TINT, "IMJD", &imjd0, NULL, &status);
    fits_read_key(fptr, TDOUBLE, "FMJD", &fmjd0, NULL, &status);
    double toffs = (double)(imjd-imjd0) + (fmjd-fmjd0);
    toffs *= 86400.0;
    char key[80];
    sprintf(key, "DT%5.5d", ispec);
    fits_write_key(fptr, TDOUBLE, key, &toffs, "[s]", &status);
    naxes[0] = naxes[1] = 1;
    naxes[2] = ispec;
    fits_write_pix(fptr, TFLOAT, naxes, 2*hf->nchan, 
            (float *)hf->data, &status);
    fits_close_file(fptr, &status);
    if (status) {
        fprintf(stderr, "Error in append_filter:\n");
        fits_report_error(stderr, status);
        return(-1);
    }
    return(0);
}

int append_profile(const char *fname, int ispec, struct profile_phase *pp) {
    int status=0;
    fitsfile *fptr;
    fits_open_file(&fptr, fname, READWRITE, &status);
    fits_movnam_hdu(fptr, IMAGE_HDU, "PROFILE", 0, &status);
    int naxis;
    long naxes[2];
    fits_get_img_dim(fptr, &naxis, &status);
    if (naxis!=2) { 
        fits_close_file(fptr, &status);
        return(-1); 
    }
    fits_get_img_size(fptr, 2, naxes, &status);
    if (naxes[0]!=pp->nphase || naxes[1]<ispec) {
        fits_close_file(fptr, &status);
        return(-2);
    }
    naxes[0] = 1;
    naxes[1] = ispec;
    fits_write_pix(fptr, TFLOAT, naxes, pp->nphase, pp->data, &status);
    fits_close_file(fptr, &status);
    if (status) {
        fprintf(stderr, "Error in append_profile:\n");
        fits_report_error(stderr, status);
        return(-1);
    }
    return(0);
}
