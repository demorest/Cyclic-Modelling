/* merit_functions.h												*/
/*																	*/
/* Created by Mark Walker on 23/08/11								*/
/*																	*/


#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "cyclic_utils.h"


/* Struct for passing the data and other information to nlopt		*/
struct cyclic_data {
    struct cyclic_spectrum *cs;	/* The data	themselves				*/
    struct profile_harm    *s0;	/* Reference profile harmonics		*/
    struct cyclic_work      *w;	/* FFTW plans, etc					*/
	int		           *rindex;	/* Array index : cimag(H(rindex))=0	*/
};

/* Functions for timing												*/
double cur_time_in_sec(); 


/* Function to compute the sum of square diffs between model and	*/
/* data, using the functional form that nlopt wants.				*/
/* Computes the gradient of the sum-of-squares if needed.			*/
/* Vector "x" contains the parameter values H(freq).				*/
double cyclic_merit_nlopt_freq(unsigned n, const double *x, 
							   double *grad, void *_data);


/* Function to compute the sum of square diffs between model and	*/
/* data, using the functional form that nlopt wants.				*/
/* Computes the gradient of the sum-of-squares if needed.			*/
/* Vector "x" contains the parameter values h(lag).					*/
double cyclic_merit_nlopt_lag(unsigned n, const double *x, 
							   double *grad, void *_data);


/* Function to compute the gradient of the sum-of-squares			*/
/* with respect to the coefficients of h(lag)						*/
/* This function added by MAW 08/11/2011 for use with NLOPT			*/
int merit_gradient_lag(struct filter_time *gradient, const CS *cs, 
						const struct filter_data *fd,
						const struct profile_harm *s0, 
						const struct cyclic_work *w);
	

/* Function to compute the gradient of the sum-of-squares			*/
/* with respect to the coefficients of H(freq)						*/
/* This function added by MAW 1/3/2012 for use with NLOPT			*/
/* Freq-gradients are computed as the FT of lag-gradients			*/
int merit_gradient_freq_via_lag(struct filter_freq *gradient, 
								const CS *cs, 
								const struct filter_data *fd,
								const struct profile_harm *s0, 
								const struct cyclic_work *w);

