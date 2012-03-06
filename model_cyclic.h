/* model_cyclic.h													*/

/* M. Walker July 2011 (not much left of PBD's original code)		*/

#ifndef _MODEL_CYCLIC_H
#define _MODEL_CYCLIC_H

#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "cyclic_utils.h"

/* Added by MAW 15/07/2011. Streamlined version of PBD model		*/
int make_model_cs(CS *cs_model, const struct filter_freq *hf,
				  const struct profile_harm *s0,
				  const struct cyclic_work *w);


/* Construct the least-squares fit of the pulse-profile (harmonics)	*/
/* to the data, given a model filter. Output s0 must be allocated	*/
/* before calling this function.  Only 1-pol for now.				*/
int optimise_profile(struct profile_harm *s0, const CS *cs,
					 const struct filter_freq *hf,
					 const struct cyclic_work *w);
	


#endif
