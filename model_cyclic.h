/* model_cyclic.h													*/

/* M. Walker July 2011 (not much left of PBD's original code)		*/

#ifndef _MODEL_CYCLIC_H
#define _MODEL_CYCLIC_H

#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "cyclic_utils.h"

/* Useful struct for holding both time/freq version of the current
 * filter, along with the 'sheared' cyclic spectra computed from 
 * the filter.  Helps avoid recomputing all this stuff more than
 * necessary.
 */
struct filter_data {
    struct filter_time ht;
    struct filter_freq hf;
    struct cyclic_spectrum cs_pos;
    struct cyclic_spectrum cs_neg;
};

/* Populate the filter_data struct assuming the time-domain filter 
 * is already filled in.
 */
void compute_filter_data_from_time(struct filter_data *fd, 
        const CS *cs_in, const struct cyclic_work *w);

/* Populate the filter_data struct assuming the freq-domain filter 
 * is already filled in.
 */
void compute_filter_data_from_freq(struct filter_data *fd, 
        const CS *cs_in, const struct cyclic_work *w);

/* Free all memory associated with filter_data struct */
void free_filter_data(struct filter_data *fd);

/* Added by MAW 15/07/2011. Streamlined version of PBD model		*/
int make_model_cs(CS *cs_model, const struct filter_data *fd,
				  const struct profile_harm *s0,
				  const struct cyclic_work *w);


/* Construct the least-squares fit of the pulse-profile (harmonics)	*/
/* to the data, given a model filter. Output s0 must be allocated	*/
/* before calling this function.  Only 1-pol for now.				*/
int optimise_profile(struct profile_harm *s0, const CS *cs,
					 const struct filter_freq *hf,
					 const struct cyclic_work *w);
	


#endif
