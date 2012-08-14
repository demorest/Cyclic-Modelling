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

