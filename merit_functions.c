/* merit_functions.c												*/
/*																	*/
/* Created by Mark Walker on 23/08/11								*/
/* Functions for calculating the figure-of-merit, and its			*/
/* gradient with respect to the real- and imaginary-parts			*/
/* of the filter coefficients.										*/
/* Patterned on original code by Paul Demorest (Dec 2010)			*/

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

#include "cyclic_utils.h"
#include "model_cyclic.h"
#include "merit_functions.h"

/* Functions for timing												*/
double cur_time_in_sec() {
    struct timeval tv;
    int rv=0;
    if ((rv=gettimeofday(&tv,NULL))!=0) { return(0); }
    return (tv.tv_sec+tv.tv_usec*1e-6);
}

/* Return the sum of square differences between current model and	*/
/* data, using the functional form that nlopt wants.				*/
/* Compute the gradient of the sum-of-squares if needed.			*/
/* Vector "x" contains the parameter values H(freq).				*/
double cyclic_merit_nlopt_freq(unsigned n, const double *x, 
							   double *grad, void *_data) {
	
	extern int verbose;
	extern int sample_ncalls;
    static int ncalls=0;
    static double tot_time = 0.0;
    double t0 = cur_time_in_sec();
    ncalls++;
	sample_ncalls++;
	
    /* Pointer to input data										*/
    struct cyclic_data *data = (struct cyclic_data *)_data;
	
    /* check dimensions												*/
	if (n != 2*(data->w->nchan)-1) {
        fprintf(stderr, 
		"cyclic_merit_nlopt_freq: error, inconsistent sizes!\n");
        exit(1);
    }
	
	/* Put the parameters x into the struct hf, for ease of use		*/
	int rchan = *(data->rindex);
	struct filter_freq hf;
	hf.nchan = data->w->nchan;
    filter_alloc_freq(&hf);
	parms2struct_freq(x, &hf, rchan);
	
    /* Compute the model cyclic spectrum from hf and s0				*/
	struct cyclic_spectrum cs_model;
	cs_copy_parms(data->cs, &cs_model);
	cyclic_alloc_cs(&cs_model);
	
	int imod = make_model_cs(&cs_model, &hf, data->s0, data->w);
	
    if (imod != 0) { fprintf(stderr, 
		"cyclic_merit_nlopt_freq: error in make_model_cs (%d)\n",imod);
        exit(1);
    }	
	
    /* Evaluate and return the sum of |data-model|^2				*/
    double merit = cyclic_square_difference(data->cs, &cs_model);
		
	/* Evaluate the gradient, if needed								*/
    if (grad != NULL) {		
		struct filter_freq complex_gradient;
		complex_gradient.nchan=data->w->nchan;
		filter_alloc_freq(&complex_gradient);
		
		int mg;
		mg = merit_gradient_freq_via_lag(&complex_gradient,
					data->cs, &hf, data->s0, data->w);
		
		if (mg != 0) { fprintf(stderr, 
		"cyclic_merit_nlopt_freq: error in merit_gradient_freq (%d)\n",
							   mg);
			exit(1);
		}
				
		/* Assign NLOPT-grad values from the complex_gradient		*/
		struct2parms_freq(&complex_gradient, grad, rchan);
				
		filter_free_freq(&complex_gradient);
    }	
	
    /* Print timing info											*/
    double t1 = cur_time_in_sec();
    tot_time += t1-t0;
	if (verbose) {
		printf("Freq Solver: ncalls=%d, avg %.2e sec/call, m = %.9e\n", 
			   ncalls, tot_time/(double)ncalls, merit);
		/* printf("%d %.9e\n", ncalls, merit);						*/
	}	
	
	cyclic_free_cs(&cs_model);
    filter_free_freq(&hf);
	
    return(merit);
}

/* Return the sum of square differences between current model and	*/
/* data, using the functional form that nlopt wants.				*/
/* Compute the gradient of the sum-of-squares if needed.			*/
/* Vector "x" contains the parameter values h(lag).					*/
double cyclic_merit_nlopt_lag(unsigned n, const double *x, 
							   double *grad, void *_data) {
	
	extern int verbose;
	extern int sample_ncalls;
    static int ncalls=0;
    static double tot_time = 0.0;
    double t0 = cur_time_in_sec();
    ncalls++;
	sample_ncalls++;
	
    /* Pointer to input data										*/
    struct cyclic_data *data = (struct cyclic_data *)_data;
	
    /* check dimensions												*/
	if (n != 2*(data->w->nlag)-1) {
        fprintf(stderr, 
			"cyclic_merit_nlopt_lag: error, inconsistent sizes!\n");
        exit(1);
    }
	
	/* Put the parameters x into the struct ht, for ease of use		*/
	int rlag = *(data->rindex);
	struct filter_time ht;
	ht.nlag = data->w->nlag;
    filter_alloc_time(&ht);
	parms2struct_time(x, &ht, rlag);
		
	/* Convert ht=h(lag) to hf=H(freq)								*/
	struct filter_freq hf;
	hf.nchan = data->w->nchan;
    filter_alloc_freq(&hf);
	filter_time2freq(&ht, &hf, data->w);
		
    /* Compute the model cyclic spectrum from hf and s0				*/
	struct cyclic_spectrum cs_model;
	cs_copy_parms(data->cs, &cs_model);
	cyclic_alloc_cs(&cs_model);
	
	int imod = make_model_cs(&cs_model, &hf, data->s0, data->w);
	
    if (imod != 0) { fprintf(stderr, 
	  "cyclic_merit_nlopt_lag: error in make_model_cs (%d)\n",imod);
        exit(1);
    }	
	
    /* Evaluate and return the sum of |data-model|^2				*/
    double merit = cyclic_square_difference(data->cs, &cs_model);
	
	
	/* Evaluate the gradient, if needed								*/
    if (grad != NULL) {		
		struct filter_time complex_gradient_time;
		complex_gradient_time.nlag=data->w->nlag;
		filter_alloc_time(&complex_gradient_time);
		
		int mg;
		mg = merit_gradient_lag(&complex_gradient_time, data->cs,
									&hf, data->s0, data->w);
		
		if (mg != 0) { fprintf(stderr, 
		   "cyclic_merit_nlopt_lag: error in gradient (%d)\n",mg);
			exit(1);
		}
				
		/* Assign NLOPT-grad values from the complex_gradient_time	*/
		struct2parms_time(&complex_gradient_time, grad, rlag);
		
		filter_free_time(&complex_gradient_time);
    }	
	
    /* Print timing info											*/
    double t1 = cur_time_in_sec();
    tot_time += t1-t0;
	if (verbose) {
		printf("Lag Solver: ncalls=%d, m = %.9e\n", ncalls, merit);				
		/* printf("%d %.9e\n", ncalls, merit);						*/
	}	
	
	cyclic_free_cs(&cs_model);
    filter_free_freq(&hf);
    filter_free_time(&ht);
	
    return(merit);
}

int merit_gradient_lag(struct filter_time *gradient, const CS *cs, 
						const struct filter_freq *hf, 
						const struct profile_harm *s0, 
						const struct cyclic_work *w) {
	
	/* This function added by MAW 08/11/2011 for use with NLOPT		*/
	/* Evaluates the gradient of the merit function with respect 	*/
	/* to the coefficients of the filter function h(lag)			*/
	
    /* Only valid for 1-pol data at present							*/
    if (cs->npol!=1 || w->nlag!=gradient->nlag ||
		cs->nchan!=hf->nchan || cs->nharm!=s0->nharm) { 
		printf("merit_gradient_lag : Incompatible dimensions\n");
		exit(1); }
	
	/* Allocate temporary CS and CC structs							*/
	struct cyclic_spectrum cs_res;
	struct cyclic_spectrum cs_tmp1;
	struct cyclic_spectrum cs_tmp2;
	cs_copy_parms(cs, &cs_res);
	cs_copy_parms(cs, &cs_tmp1);
	cs_copy_parms(cs, &cs_tmp2);
	cyclic_alloc_cs(&cs_res);
	cyclic_alloc_cs(&cs_tmp1);
	cyclic_alloc_cs(&cs_tmp2);
	
	struct cyclic_correlation cc1;
	cc1.nharm     = cs->nharm;
    cc1.npol      = cs->npol;
    cc1.nlag      = gradient->nlag;
    cc1.imjd      = cs->imjd;
    cc1.fmjd      = cs->fmjd;
    cc1.ref_phase = cs->ref_phase;
    cc1.ref_freq  = cs->ref_freq;
    cc1.rf        = cs->rf;
    cc1.bw        = cs->bw;
	cyclic_alloc_cc(&cc1);
	
	
	/* Construct the current model cyclic spectrum from hf and s0	*/
	make_model_cs(&cs_res, hf, s0, w);
		
	/* And form the residual ( = model - data )						*/
	cs_subtract(cs, &cs_res);
	
	/* Need to handle +ve alpha and -ve alpha separately, so copy	*/
	/* the residual cyclic spectrum to save re-computing			*/
	cs_copy_data(&cs_res, &cs_tmp1);
	
	/* Positive values of alpha first								*/
	/* Make a filter array											*/
	filter2cs(&cs_tmp2, hf);
	/* And shear the filter array									*/
	double shear = -0.5;
	cyclic_shear_cs(&cs_tmp2, shear, w);
	/* Form the product of residuals with sheared filters			*/
	cs_multiply(&cs_tmp1, &cs_tmp2);
	/* Convert from cyclic-spectrum to cyclic-correlation			*/
	cyclic_cs2cc(&cs_tmp2, &cc1, w);
	
	/* For each lag, sum over all alpha								*/
	/* Only one polarisation at present								*/
	int ih, ilag, ip=0;
	double tau, phs;
	fftwf_complex phasor;
	for (ilag=0; ilag<cc1.nlag; ilag++) {
		gradient->data[ilag] = 0.0 + I * 0.0;
		int lag = (ilag<=cc1.nlag/2) ? ilag : ilag-cc1.nlag;		
		tau = (double)lag * (double)cs->nchan / 
		                  ( (double)cc1.nlag * cc1.bw*1.e6 );
		for (ih=1; ih<cc1.nharm; ih++) {
			phs = M_PI * tau * (double)ih * cc1.ref_freq;
			phasor = cos(phs)+I*sin(phs);
			fftwf_complex *ccval = get_cc(&cc1,ih,ip,ilag);
			gradient->data[ilag] += 4.0 * (*ccval) * phasor
			 * conj(s0->data[ih]) / (float)cs->nchan;
		}
	}
	
	/* Now repeat the above for negative values of alpha			*/		
	/* The CS for -ve alpha is conjugate of that for +ve alpha		*/
	cs_copy_data(&cs_res, &cs_tmp1);
	cs_conjugate(&cs_tmp1);
	/* Make a filter array											*/
	filter2cs(&cs_tmp2, hf);
	/* And shear the filter array									*/
	shear = 0.5;
	cyclic_shear_cs(&cs_tmp2, shear, w);
	/* Form the product of residuals with sheared filters			*/
	cs_multiply(&cs_tmp1, &cs_tmp2);
	/* Convert from cyclic-spectrum to cyclic-correlation			*/
	cyclic_cs2cc(&cs_tmp2, &cc1, w);
	
	/* For each lag, sum over all alpha								*/
	/* Only one polarisation at present								*/
	for (ilag=0; ilag<cc1.nlag; ilag++) {
		int lag = (ilag<=cc1.nlag/2) ? ilag : ilag-cc1.nlag;		
		tau = (double)lag * (double)cs->nchan / 
		( (double)cc1.nlag * cc1.bw*1.e6 );
		for (ih=1; ih<cc1.nharm; ih++) {
			phs = M_PI * tau * (double)ih * cc1.ref_freq;
			phasor = cos(phs)-I*sin(phs);
			fftwf_complex *ccval = get_cc(&cc1,ih,ip,ilag);
			gradient->data[ilag] += 4.0 * (*ccval) * phasor
				* s0->data[ih] / (float)cs->nchan;
		}
	}
	
	/* Free-up the arrays allocated herein							*/
	cyclic_free_cs(&cs_res);
	cyclic_free_cs(&cs_tmp1);
	cyclic_free_cs(&cs_tmp2);	
	cyclic_free_cc(&cc1);	
    return(0);
}



int merit_gradient_freq_via_lag(struct filter_freq *gradient,  
								const CS *cs, 
								const struct filter_freq *hf, 
								const struct profile_harm *s0, 
								const struct cyclic_work *w) {
	
	/* This function added by MAW 1/3/2012 for use with NLOPT		*/
	/* Evaluates the gradient of the merit function with respect 	*/
	/* to the coefficients of the filter function H(freq)			*/
	/* Freq-gradients are computed as the FT of lag-gradients		*/
	
    /* Only valid for 1-pol data at present							*/
    if (cs->npol!=1 || w->nchan!=gradient->nchan ||
		cs->nchan!=hf->nchan || cs->nharm!=s0->nharm) { 
		printf(
		"merit_gradient_freq_via_lag: Incompatible dimensions\n");
		exit(1);
	}
	
	
	struct filter_time lag_gradient;
	lag_gradient.nlag=w->nlag;
	filter_alloc_time(&lag_gradient);
	
	int mg;
	mg = merit_gradient_lag(&lag_gradient, cs, hf, s0, w);
	
	if (mg != 0) { fprintf(stderr, 
	"merit_gradient_freq_via_lag:error in merit_gradient_lag (%d)\n",
						   mg);
		exit(1);
	}
	
	/* Now convert the frequency-gradient to the lag-gradient		*/
	/* filter_freq_renorm(&freq_gradient);							*/
	filter_time2freq(&lag_gradient, gradient, w);
	
	
	/* Free-up the arrays allocated herein							*/
	filter_free_time(&lag_gradient);
	
    return(0);
}

