/*
 * test_merit.c
 *
 *  Created on: Nov 5, 2012
 *      Author: gjones
 *
 *  This file performs a series of tests using the functions defined in the Cyclic-Modelling package on random data
 *  The random test vectors and the results are written out to simple txt files that can be ingested by pycyc
 *  test routines to compare the two programs and ensure that they are functionally identical.
 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <signal.h>
#include <getopt.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <fitsio.h>
#include <nlopt.h>
#include <string.h>

/*#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h> */

#include "cyclic_utils.h"
#include "model_cyclic.h"
#include "merit_functions.h"

#define randf() ((((float)rand())/((float)RAND_MAX)) - 0.5)

int verbose = 0;
int sample_ncalls = 0;

#define NLAG 1024
int main(int argc, char *argv[]) {
    struct cyclic_work w;
    int fft_threads = 4;
    int nspec;
    int rv;
    int k,m,n;
    int rindex;
    float error;
    w.nphase = 512;
	w.npol   = 1;
	w.nchan  = NLAG;
	nspec     = 1;

	w.nlag   = 0;
	w.nharm  = 0;

	float bw = 1.0;
	float ref_freq = 250.0;
	float rf = 430.0;



    /* Initialise FFTs												*/
    fftwf_init_threads();
    fftwf_plan_with_nthreads(fft_threads);
    if (verbose) { printf("Planning FFTs\n"); fflush(stdout); }
#define WF "cyclic_wisdom.dat"
    FILE *wf = fopen(WF,"r");
    if (wf!=NULL) { fftwf_import_wisdom_from_file(wf); fclose(wf); }

    rv = cyclic_init_ffts(&w); //sets nharm, nlag

    if (rv) {
        fprintf(stderr, "Error planning ffts (rv=%d)\n", rv);
        exit(1);
    }
    wf = fopen(WF,"w");
    if (wf!=NULL) { fftwf_export_wisdom_to_file(wf); fclose(wf); }


    /* Allocate some stuff											*/
    struct periodic_spectrum raw;
    struct cyclic_spectrum cs, cs2;
    struct filter_freq hf, hf_prev;
	struct filter_time ht, ht_grad;
    struct profile_phase pp, pp_ref, pp_int;
    struct profile_harm  ph, ph_ref;

    raw.nphase = pp.nphase = pp_ref.nphase = pp_int.nphase = w.nphase;
    raw.nchan = cs.nchan = cs2.nchan = hf.nchan = hf_prev.nchan = w.nchan;
    cs.nharm = cs2.nharm = ph.nharm = ph_ref.nharm = w.nharm;
	ht.nlag = ht_grad.nlag = w.nlag;
    raw.npol = 1;
    cs.npol = cs2.npol = 1;

    cs.ref_freq = cs2.ref_freq = ref_freq;
    cs.bw = cs2.bw = bw;
    cs.rf = cs2.rf = rf;


    double grad[2*NLAG-1]; // assumes nlag = NLAG!
    double x[2*NLAG-1];
    double merit;

    cyclic_alloc_ps(&raw);
    cyclic_alloc_cs(&cs);
    cyclic_alloc_cs(&cs2);
    filter_alloc_freq(&hf);
    filter_alloc_freq(&hf_prev);
	filter_alloc_time(&ht);
	filter_alloc_time(&ht_grad);
    profile_alloc_phase(&pp);
    profile_alloc_phase(&pp_ref);
    profile_alloc_phase(&pp_int);
    profile_alloc_harm(&ph);
    profile_alloc_harm(&ph_ref);


    /*
     * PHASE 1:
     * First load hf with random data, and write it to hf.tst
     * then load pp with random data and write it to pp.tst
  	 * Compute ph from pp and write to ph.tst
  	 * Recompute pp_ref from ph to make sure the result is identical with original pp. write pp_ref to pp2pp.tst
	 * using ph and hf, compute a model cs
	 * Do the same for hf_prev, pp_ref, ph_ref, and cs2
	 *
	 * PHASE 2:
	 * Now compute best fit ph_ref given hf and cs
	 * Write reconstructed ph_ref to ph_reconstructed.tst
	 * Then compare reconstructed ph_ref to original ph using profile_ms_difference
	 *
	 * PHASE 3:
	 * Next compute ht from hf
	 * Find rindex as max filter time
	 * Transfer ht to parameter array x
	 * Compute lag merit and gradient from ht, cs, and ph. Gradient stored in grad1.tst
	 * Since ht came directly from hf used to generate cs, this should give good (low) merit and small gradient
	 *
	 * Now compute ht from hf_prev (a random transfer function unrelated to cs)
	 * Compute lag merit and gradient from ht, cs, and ph. Gradient stored in grad2.tst
	 * Since ht came from unrelated hf_prev, this should give poor (high) merit and large gradient
	 *
	 * Scalar parameters are saved to params.tst
     */

    FILE *params = fopen("params.tst", "w");
    fprintf(params,"dict(\n");
    fprintf(params,"bw= %f,\n",bw);
    fprintf(params,"ref_freq= %f,\n",ref_freq);
    fprintf(params,"rf= %f,\n",rf);


    /*
     * PHASE 1
     */
    printf("Loading hf\n");
    FILE *hff = fopen("hf.tst","w");
    float re,im;
    for(n=0; n < hf.nchan; n++) {
    	re = randf();
    	im = randf();
    	hf.data[n] = re + I * im;
    	fprintf(hff,"%.15g, %.15g\n",re,im);
    }
    fclose(hff);

    printf("Loading profile\n");

    FILE *ppf = fopen("pp.tst","w");

	for(n=0; n < pp.nphase; n++) {
		pp.data[n] = randf();
		fprintf(ppf,"%.15g\n",pp.data[n]);
	}
	fclose(ppf);

	printf("computing phase2harm\n");
	profile_phase2harm(&pp,&ph,&w);

	FILE *phf = fopen("ph.tst","w");

	for(n=0; n < ph.nharm; n++) {
		fprintf(phf,"%.15g, %.15g\n",creal(ph.data[n]),cimag(ph.data[n]));
	}
	fclose(phf);

	profile_harm2phase(&ph,&pp_ref,&w);
	phf = fopen("pp2pp.tst","w");

	for(n=0; n < pp_ref.nphase; n++) {
		fprintf(phf,"%.15g\n",pp_ref.data[n]);
	}
	fclose(phf);

	printf("computing CS from hf and ph\n");
	rv = make_model_cs(&cs, &hf, &ph, &w);


    printf("Loading hf_prev\n");
    hff = fopen("hf_prev.tst","w");

    for(n=0; n < hf_prev.nchan; n++) {
    	re = randf();
    	im = randf();
    	hf_prev.data[n] = re + I * im;
    	fprintf(hff,"%.15g, %.15g\n",re,im);
    }
    fclose(hff);

    printf("Loading profile ref\n");

    ppf = fopen("pp_ref.tst","w");

	for(n=0; n < pp_ref.nphase; n++) {
		pp_ref.data[n] = randf();
		fprintf(ppf,"%.15g\n",pp_ref.data[n]);
	}
	fclose(ppf);

	printf("computing phase2harm for pp_ref\n");
	profile_phase2harm(&pp_ref,&ph_ref,&w);

	phf = fopen("ph_ref.tst","w");

	for(n=0; n < ph_ref.nharm; n++) {
		fprintf(phf,"%.15g, %.15g\n",creal(ph_ref.data[n]),cimag(ph_ref.data[n]));
	}
	fclose(phf);

	printf("computing CS2\n");
	printf("cs2.nharm: %d, cs2.npol: %d, cs2.nchan %d\n", cs2.nharm, cs2.npol, cs2.nchan);

	rv = make_model_cs(&cs2, &hf_prev, &ph_ref, &w);
	if (rv) {
		printf("error creating model cs2 %d\n", rv);
	}

	/*
	 * PHASE 2
	 */
	printf("now computing best profile\n");

	rv = optimise_profile(&ph_ref, &cs, &hf, &w);
	if (rv) {
		printf("error calculating inverse profile %d\n", rv);
	}

	phf = fopen("ph_reconstructed.tst","w");

	for(n=0; n < ph_ref.nharm; n++) {
		fprintf(phf,"%.15g, %.15g\n",creal(ph_ref.data[n]),cimag(ph_ref.data[n]));
	}
	fclose(phf);

	error = profile_ms_difference(&ph, &ph_ref,ph.nharm);
	printf("Error between initial profile and reconstructed profile: %.15g ... ", error);
	if (error < 1e-10) {
		printf("PASS!\n");
	}
	else {
		printf("FAIL!\n");
	}


	/*
	 * PHASE 3
	 */
	filter_freq2time(&hf,&ht,&w);

	rindex = maximum_filter_time(&ht);

	printf("Max filter time (rindex): %d\t",rindex);
	fprintf(params,"rindex1= %d,\n",rindex);

	struct cyclic_data cdata;
	cdata.cs     = &cs;
	cdata.s0     = &ph;
	cdata.w	     = &w;
	cdata.rindex = &rindex;

	struct2parms_time(&ht,x,rindex);
	merit =  cyclic_merit_nlopt_lag(2*NLAG-1, x,grad, &cdata);
	printf("merit: %.15g\n",merit);

	fprintf(params,"merit1= %.15g,\n",merit);

	phf = fopen("grad1.tst","w");

	for(n=0; n < 2*NLAG-1; n++) {
		fprintf(phf,"%.15g\n",grad[n]);
	}
	fclose(phf);


	filter_freq2time(&hf_prev,&ht,&w);

	rindex = maximum_filter_time(&ht);

	printf("Max filter time (rindex): %d\t",rindex);
	fprintf(params,"rindex2= %d,\n",rindex);

	struct2parms_time(&ht,x,rindex);
	merit =  cyclic_merit_nlopt_lag(2*NLAG-1, x,grad, &cdata);

	phf = fopen("grad2.tst","w");

	for(n=0; n < 2*NLAG-1; n++) {
		fprintf(phf,"%.15g\n",grad[n]);
	}
	fclose(phf);

	printf("merit: %.15g,\n",merit);
	fprintf(params,"merit2= %.15g,\n",merit);

	fprintf(params,")");
	fclose(params);

	cyclic_free_ps(&raw);
	cyclic_free_cs(&cs);
	cyclic_free_cs(&cs2);
	filter_free_freq(&hf);
	filter_free_freq(&hf_prev);
	filter_free_time(&ht);
	filter_free_time(&ht_grad);
	profile_free_phase(&pp);
	profile_free_harm(&ph);
	profile_free_phase(&pp_ref);
	profile_free_harm(&ph_ref);
	profile_free_phase(&pp_int);
	cyclic_free_ffts(&w);
	exit(0);

}
