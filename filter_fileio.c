/* filter_fileio.c
 *
 * Routines to read/write filter funcs and profiles.
 * Split off from cyclic_utils.c, PBD 2012/08/14.
 */
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>

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
