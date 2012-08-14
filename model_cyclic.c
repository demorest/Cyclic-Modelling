/* model_cyclic.c													*/

/* M. Walker July 2011 (not much left of PBD's original code)		*/

#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "cyclic_utils.h"
#include "model_cyclic.h"

/* Compute freq-domain filter and sheared CS */
void compute_filter_data_from_time(struct filter_data *fd, 
        const CS *cs_in, const struct cyclic_work *w) {

	/* Convert ht=h(lag) to hf=H(freq)								*/
	fd->hf.nchan = w->nchan;
    filter_alloc_freq(&fd->hf);
	filter_time2freq(&fd->ht, &fd->hf, w);

    /* Positive shear version
     * Note, "positive" refers to sign of the "shear" argument.
     * It looks like this sign is opposite to the sign of alpha
     * that the sheared version applies to.
     */
    cs_copy_parms(cs_in, &fd->cs_pos);
	cyclic_alloc_cs(&fd->cs_pos);
	filter2cs(&fd->cs_pos, &fd->hf);
	double shear = 0.5;
	cyclic_shear_cs(&fd->cs_pos, shear, w);
	
    /* Negative shear version */
    cs_copy_parms(cs_in, &fd->cs_neg);
	cyclic_alloc_cs(&fd->cs_neg);
	filter2cs(&fd->cs_neg, &fd->hf);
	shear = -0.5;
	cyclic_shear_cs(&fd->cs_neg, shear, w);
}

/* Compute time-domain filter and sheared CS */
void compute_filter_data_from_freq(struct filter_data *fd, 
        const CS *cs_in, const struct cyclic_work *w) {

	/* Convert hf=H(freq) to ht=h(lag)							*/
	fd->ht.nlag = w->nlag;
    filter_alloc_time(&fd->ht);
	filter_freq2time(&fd->hf, &fd->ht, w);
		
    /* Positive shear version
     * Note, "positive" refers to sign of the "shear" argument.
     * It looks like this sign is opposite to the sign of alpha
     * that the sheared version applies to.
     */
    cs_copy_parms(cs_in, &fd->cs_pos);
	cyclic_alloc_cs(&fd->cs_pos);
	filter2cs(&fd->cs_pos, &fd->hf);
	double shear = 0.5;
	cyclic_shear_cs(&fd->cs_pos, shear, w);
	
    /* Negative shear version */
    cs_copy_parms(cs_in, &fd->cs_neg);
	cyclic_alloc_cs(&fd->cs_neg);
	filter2cs(&fd->cs_neg, &fd->hf);
	shear = -0.5;
	cyclic_shear_cs(&fd->cs_neg, shear, w);
}

/* Free up filter_data memory */
void free_filter_data(struct filter_data *fd) {
    filter_free_time(&fd->ht);
    filter_free_freq(&fd->hf);
    cyclic_free_cs(&fd->cs_pos);
    cyclic_free_cs(&fd->cs_neg);
}

int make_model_cs(CS *cs_model, const struct filter_data *fd,
				  const struct profile_harm *s0,
				  const struct cyclic_work *w) {

/* Construct a model cyclic spectrum given H(freq) and reference	*/
/* pulse-profile harmonics. Only 1-pol for now.						*/
/* This function returns the model CS for alpha >= 0				*/
/* For -ve alpha, take the complex conjugate of these values		*/
	
    /* Check dimensions												*/
    if (cs_model->npol!=1 || cs_model->nchan!=fd->hf.nchan || 
		cs_model->nharm!=s0->nharm) { 
		printf("make_model_cs : Incompatible dimensions\n");
        printf("cs_model->nchan=%d fd->hf.nchan=%d\n", 
                cs_model->nchan, fd->hf.nchan);
        printf("cs_model->nharm=%d s0->nharm=%d\n", 
                cs_model->nharm, s0->nharm);
		return(-1); }
	
	/* Create the intrinsic cyclic spectrum							*/
	/* Currently this has no dependence on radio frequency			*/
	profile2cs(cs_model, s0);
	/* Propagate this spectrum through the filter in two steps		*/
		
	/* Multiply the intrinsic spectrum by this sheared array		*/
	cs_multiply(&fd->cs_pos, cs_model);

	cs_multiply_conj(&fd->cs_neg, cs_model);
	
	/* Replace the ends of the model CS with zeros					*/
	cyclic_padding(cs_model);
		
	/* Free-up temporary cyclic spectrum and return					*/
    return(0);
}



int optimise_profile(struct profile_harm *s0, const CS *cs,
					 const struct filter_freq *hf,
					 const struct cyclic_work *w) {

/* Construct the optimum (least-squares) estimate of the pulse		*/
/* profile (harmonics), given the data and a model filter.			*/
/* Only 1-pol for now.												*/
		
    /* Check dimensions												*/
    if (cs->npol!=1 || cs->nchan!=hf->nchan || cs->nharm!=s0->nharm) { 
		printf("optimise_profile : Incompatible dimensions\n");
		return(-1); }
			
	/* Allocate two temporary cyclic spectral arrays				*/
	struct cyclic_spectrum cs_tmp1;
	struct cyclic_spectrum cs_tmp2;
	cs_copy_parms(cs, &cs_tmp1);
	cs_copy_parms(cs, &cs_tmp2);
	cyclic_alloc_cs(&cs_tmp1);
	cyclic_alloc_cs(&cs_tmp2);
		
	/* Determine the (conjugate of) the effect of the filter on		*/
	/* the CS. Two steps. First form the +ve shifted filter array.	*/
	/* Create a "cyclic spectrum" array filled with filter hf		*/
	filter2cs(&cs_tmp1, hf);
	/* Shear the filter spectrum by +alpha/2						*/
	double shear = 0.5;
	cyclic_shear_cs(&cs_tmp1, shear, w);
	/* Form the complex conjugate									*/
	cs_conjugate(&cs_tmp1);
	
	/* Now form the -ve shifted filter array						*/
	/* Recreate a "cyclic spectrum" array filled with filter hf		*/
	filter2cs(&cs_tmp2, hf);
	/* Shear the filter spectrum by -alpha/2						*/
	shear = -0.5;
	cyclic_shear_cs(&cs_tmp2, shear, w);
	/* And multiply by the +ve shifted array						*/
	cs_multiply(&cs_tmp1, &cs_tmp2);

	/* Copy the result into cs_tmp1									*/
	cs_copy_data(&cs_tmp2, &cs_tmp1);
	/* Form the complex conjugate									*/
	cs_conjugate(&cs_tmp2);
	
	/* Now form the arrays cs H(-)H(+)*  and |H(-)|^2 |H(+)|^2		*/
	cs_multiply(&cs_tmp1, &cs_tmp2);
	cs_multiply(cs, &cs_tmp1);

	/* Scrunch both CS arrays in the frequency direction			*/
	struct profile_harm s_tmp;
	s_tmp.nharm = s0->nharm;
	profile_alloc_harm(&s_tmp);
	cyclic_fscrunch_cs(s0, &cs_tmp1);
	cyclic_fscrunch_cs(&s_tmp, &cs_tmp2);

	/* Ratio these two arrays to get the optimised profile			*/
	s0->data[0] = 0.0 + I * 0.0;
	float denominator = 0.0;
	int ih;
	for (ih=1; ih<s0->nharm; ih++) {
		denominator = creal(s_tmp.data[ih]);
		if (denominator>0.0) {
			s0->data[ih] /= denominator;
		}
		else {
			s0->data[ih] = 0.0 + I * 0.0;
		}
	}
	
	/* Free-up temporary arrays and return							*/
	profile_free_harm(&s_tmp);
	cyclic_free_cs(&cs_tmp1);
	cyclic_free_cs(&cs_tmp2);	
    return(0);
}
