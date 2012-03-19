/* cyclic_utils.c													*/

/* P. Demorest, December 2010 with additions by M. Walker July 2011	*/


#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>

#include "cyclic_utils.h"

/* Allocs / frees */
void cyclic_alloc_ps(PS *d) {
    d->data = (float *)fftwf_malloc(sizeof(float) * 
            d->nphase * d->nchan * d->npol);
}
void cyclic_alloc_cs(CS *d) {
    d->data = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * 
            d->nharm * d->nchan * d->npol);
}
void cyclic_alloc_cc(CC *d) {
    d->data = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * 
            d->nharm * d->nlag * d->npol);
}
void cyclic_alloc_pc(PC *d) {
    d->data = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * 
            d->nphase * d->nlag * d->npol);
}
void cyclic_free_ps(PS *d) { fftwf_free(d->data); }
void cyclic_free_cs(CS *d) { fftwf_free(d->data); }
void cyclic_free_cc(CC *d) { fftwf_free(d->data); }
void cyclic_free_pc(PC *d) { fftwf_free(d->data); }

void filter_alloc_time(struct filter_time *f) {
    f->data = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) *
            f->nlag);
}
void filter_alloc_freq(struct filter_freq *f) {
    f->data = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) *
            f->nchan);
}
void filter_free_time(struct filter_time *f) { fftwf_free(f->data); }
void filter_free_freq(struct filter_freq *f) { fftwf_free(f->data); }

void profile_alloc_phase(struct profile_phase *f) { 
    f->data = (float *)fftwf_malloc(sizeof(float) * f->nphase); 
}
void profile_alloc_harm(struct profile_harm *f) { 
    f->data = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * f->nharm); 
}
void profile_free_phase(struct profile_phase *f) { fftwf_free(f->data); }
void profile_free_harm(struct profile_harm *f) { fftwf_free(f->data); }

/* Set up fft plans.  Need to have npol, nphase, nchan 
 * already filled in struct */
int cyclic_init_ffts(struct cyclic_work *w) {

    /* Infer lag, harmonic sizes from chan/phase */
    w->nlag = w->nchan; // Total number of lags including + and -
    w->nharm = w->nphase/2 + 1; // Only DC and positive harmonics

    /* Alloc temp arrays for planning */
    PS ps;
    CS cs;
    CC cc;
    PC pc;
    struct filter_time ft;
    struct filter_freq ff;
    struct profile_phase pp;
    struct profile_harm ph;

    ps.npol = cs.npol = cc.npol = pc.npol = w->npol;
    ps.nphase = pc.nphase = w->nphase;
    ps.nchan = cs.nchan = w->nchan;
    cs.nharm = cc.nharm = w->nharm;
    cc.nlag = pc.nlag = w->nlag;

    cyclic_alloc_ps(&ps);
    cyclic_alloc_cs(&cs);
    cyclic_alloc_cc(&cc);
    cyclic_alloc_pc(&pc);

    ft.nlag = w->nlag;
    ff.nchan = w->nchan;
    pp.nphase = w->nphase;
    ph.nharm = w->nharm;

    filter_alloc_time(&ft);
    filter_alloc_freq(&ff);
    profile_alloc_phase(&pp);
    profile_alloc_harm(&ph);

    /* FFT plans */
    int rv=0;

    /* ps2cs - r2c fft along phase (fastest) axis */
    w->ps2cs = fftwf_plan_many_dft_r2c(1, &w->nphase, w->npol*w->nchan,
            ps.data, NULL, 1, w->nphase,
            cs.data, NULL, 1, w->nharm,
            FFTW_MEASURE | FFTW_PRESERVE_INPUT);
    if (w->ps2cs == NULL) rv++; 

    /* cs2cc - c2c ifft along channel axis */
    w->cs2cc = fftwf_plan_many_dft(1, &w->nchan, w->npol*w->nharm,
            cs.data, NULL, w->nharm*w->npol, 1,
            cc.data, NULL, w->nharm*w->npol, 1,
            FFTW_BACKWARD, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
    if (w->cs2cc == NULL) rv++; 
    
    /* cc2cs - c2c fft along lag axis */
    w->cc2cs = fftwf_plan_many_dft(1, &w->nlag, w->npol*w->nharm,
            cc.data, NULL, w->nharm*w->npol, 1,
            cs.data, NULL, w->nharm*w->npol, 1,
            FFTW_FORWARD, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
    if (w->cc2cs == NULL) rv++; 

    /* time2freq, freq2time for filters */
    w->time2freq = fftwf_plan_dft_1d(w->nlag, ft.data, ff.data,
            FFTW_FORWARD, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
    if (w->time2freq == NULL) rv++;
    w->freq2time = fftwf_plan_dft_1d(w->nchan, ff.data, ft.data,
            FFTW_BACKWARD, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
    if (w->freq2time == NULL) rv++;

    /* phase2harm, harm2phase for profiles */
    w->phase2harm = fftwf_plan_dft_r2c_1d(w->nphase, pp.data, ph.data, 
            FFTW_MEASURE | FFTW_PRESERVE_INPUT);
    if (w->phase2harm == NULL) rv++;
    w->harm2phase = fftwf_plan_dft_c2r_1d(w->nphase, ph.data, pp.data, 
            FFTW_MEASURE | FFTW_PRESERVE_INPUT);
    if (w->harm2phase == NULL) rv++;

    cyclic_free_ps(&ps);
    cyclic_free_cs(&cs);
    cyclic_free_cc(&cc);
    cyclic_free_pc(&pc);

    filter_free_time(&ft);
    filter_free_freq(&ff);
    profile_free_phase(&pp);
    profile_free_harm(&ph);

    return(rv);
}

void cyclic_free_ffts(struct cyclic_work *w) {
    if (w->ps2cs!=NULL) fftwf_destroy_plan(w->ps2cs);
    if (w->cs2cc!=NULL) fftwf_destroy_plan(w->cs2cc);
    if (w->cc2cs!=NULL) fftwf_destroy_plan(w->cc2cs);
}

int cyclic_pscrunch_ps(PS *d, float xgain, float ygain) {
	
	/* MAW 7/7/2011	Modified version of PBD's code which fixes the	*/
	/* polarisation interleave problem								*/
	
    if (d->npol<2) { return(-1); }
	
    int ichan, iphase, ichan2;
    float *xx, *yy, *xxyy;
    for (ichan=0; ichan<d->nchan; ichan++) {
		ichan2 = ichan / 2;
        for (iphase=0; iphase<d->nphase; iphase++) {
            xx   = get_ps(d, iphase, 0, ichan);
            yy	 = get_ps(d, iphase, 1, ichan);
            xxyy = get_ps(d, iphase, 0, ichan2);
			xxyy += (ichan % 2) * (d->nphase);
            *xxyy = xgain * (*xx) + ygain * (*yy);
        }
    }
	
    d->npol = 1;
    return(0);
}

int cyclic_fscrunch_ps(struct profile_phase *out, PS *in) {

    /* Only 1 pol for now */
    if (in->npol>1) return(-1); 

    int iphase, ichan;
    for (iphase=0; iphase<in->nphase; iphase++) {
        out->data[iphase] = 0.0;
        for (ichan=0; ichan<in->nchan; ichan++) {
            const float *tmp = get_ps(in,iphase,0,ichan);
            out->data[iphase] += *tmp;
        }
        out->data[iphase] /= (float)in->nchan;
    }

    return(0);
}

void cyclic_ps2cs(PS *in, CS *out, const struct cyclic_work *w) {
	
    fftwf_execute_dft_r2c(w->ps2cs, in->data, out->data);
	
    out->imjd = in->imjd;
    out->fmjd = in->fmjd;
    out->ref_phase = in->ref_phase;
    out->ref_freq = in->ref_freq;
    out->rf = in->rf;
    out->bw = in->bw;
	
	cyclic_ps2cs_renorm(out);
	
}
void cyclic_cs2cc(CS *in, CC *out, const struct cyclic_work *w) {
	
    fftwf_execute_dft(w->cs2cc, in->data, out->data);
	
    out->imjd = in->imjd;
    out->fmjd = in->fmjd;
    out->ref_phase = in->ref_phase;
    out->ref_freq = in->ref_freq;
    out->rf = in->rf;
    out->bw = in->bw;
	
}
void cyclic_cc2cs(CC *in, CS *out, const struct cyclic_work *w) {
	
    fftwf_execute_dft(w->cc2cs, in->data, out->data);
	
    out->imjd = in->imjd;
    out->fmjd = in->fmjd;
    out->ref_phase = in->ref_phase;
    out->ref_freq = in->ref_freq;
    out->rf = in->rf;
    out->bw = in->bw;
	
	cyclic_cc2cs_renorm(out);
	
}

void filter_time2freq(struct filter_time *in, struct filter_freq *out, 
        const struct cyclic_work *w) {
	
    fftwf_execute_dft(w->time2freq, in->data, out->data);
	filter_freq_renorm(out);
							
}

void filter_freq2time(struct filter_freq *in, struct filter_time *out, 
        const struct cyclic_work *w) {
	
    fftwf_execute_dft(w->freq2time, in->data, out->data);
	
}

void profile_phase2harm(struct profile_phase *in,
						struct profile_harm *out, 
						const struct cyclic_work *w) {
	
    fftwf_execute_dft_r2c(w->phase2harm, in->data, out->data);
	profile_harm_renorm(out);
}

void profile_harm2phase(struct profile_harm *in,
						struct profile_phase *out, 
						const struct cyclic_work *w) {
	
    fftwf_execute_dft_c2r(w->harm2phase, in->data, out->data);
	
}


double cyclic_square_difference(const CS *cs1, const CS *cs2) {
	/* Excludes ih = 0  from the summation							*/
	/* Computation is for +ve alpha only. 							*/
	/* Negative alpha contribution is same, so x2 at end			*/
	
    /* Direct sum, about 2x faster than old version, maybe 
     * more accurate also?
     * Only handles 1-poln for now.
     */
    double sum = 0.0;
    int ic, ih, chan_min, chan_max;
    for (ih=1; ih<cs1->nharm; ih++) {
        chan_limits_cs(&chan_min, &chan_max, ih, cs1);
        for (ic=chan_min; ic<chan_max; ic++) {
            fftwf_complex *d1 = get_cs(cs1, ih, 0, ic);
            fftwf_complex *d2 = get_cs(cs2, ih, 0, ic);
            fftwf_complex diff = *d1 - *d2;
            sum += (double)creal(diff)*(double)creal(diff)
                + (double)cimag(diff)*(double)cimag(diff);
        }
    }
    sum *= 2.0;

    return(sum);
}

double profile_ms_difference(struct profile_harm *p1,
							 struct profile_harm *p2, int max_harm) {
    double sum = 0.0;
    int i;
    for (i=1; i<max_harm; i++) {
        fftwf_complex diff = p1->data[i] - p2->data[i];
        sum += creal(diff*conj(diff));
    }
    return(sum);
}

double filter_ms_difference(struct filter_time *f1,
							struct filter_time *f2) {
    double sum = 0.0;
    int i;
    for (i=1; i<f1->nlag; i++) {
        fftwf_complex diff = f1->data[i] - f2->data[i];
        sum += creal(diff*conj(diff));
    }
    return(sum);
}

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

int	maximum_cs1(const CS *cs) {
	/* Routine added by MAW 13/07/2011								*/
	/* Returns the channel having the largest absolute value		*/
	/* of the cyclic spectrum at ih=1								*/
	
	int ic, ic_max=0, chan_min, chan_max, ih=1;
	chan_limits_cs(&chan_min, &chan_max, ih, cs);
	float cs_abs_sq_max=0.;
	for (ic=chan_min; ic<chan_max; ic++) {
		fftwf_complex *d1 = get_cs(cs, ih, 0, ic);
		float cs_abs_sq = creal((*d1) * conj(*d1));
		if (cs_abs_sq>cs_abs_sq_max) {
			ic_max = ic;
			cs_abs_sq_max = cs_abs_sq;
		}
	}
			
	return(ic_max);
}

float rms_cs(const CS *cs, const int ih) {
	
	/* Routine added by MAW 14/07/2011								*/
	/* Returns the rms of the cyclic spectrum at harmonic number ih	*/
	
	int ic, chan_min, chan_max;
	chan_limits_cs(&chan_min, &chan_max, ih, cs);
	float csrms=0.;
	for (ic=chan_min; ic<chan_max; ic++) {
		fftwf_complex *d1 = get_cs(cs, ih, 0, ic);
		csrms += creal((*d1) * conj(*d1));
	}
	csrms /= (float)(chan_max-chan_min);
	
	return(sqrt(csrms));
}

int	cs_copy_data(const CS *cs1, CS *cs2) {
	/* Routine added by MAW 13/07/2011								*/
	/* Copies cyclic spectrum cs1 to cs2							*/
	
	/* Check dimensions												*/
	if (cs1->nharm!=cs2->nharm || cs1->nchan!=cs2->nchan 
		|| cs1->npol!=cs2->npol) {
		printf("Error : cs_copy : incompatible dimensions\n");
		exit(1);
	}
	
    int ih, ic, ip;
    for (ic=0; ic<cs1->nchan; ic++) {
        for (ip=0; ip<cs1->npol; ip++) {
            for (ih=0; ih<cs1->nharm; ih++) {
                fftwf_complex *d1 = get_cs(cs1, ih, ip, ic);
                fftwf_complex *d2 = get_cs(cs2, ih, ip, ic);
                *d2 = *d1;
            }
        }
    }
    return(0);
}

int	cs_copy_parms(const CS *cs1, CS *cs2) {
	/* Routine added by MAW 16/07/2011								*/
	/* Copies cyclic spectrum parameters from cs1 to cs2			*/
	
	cs2->nchan     = cs1->nchan;
	cs2->nharm     = cs1->nharm;
	cs2->npol      = cs1->npol;
	cs2->imjd      = cs1->imjd;
    cs2->fmjd      = cs1->fmjd;
    cs2->ref_phase = cs1->ref_phase;
    cs2->ref_freq  = cs1->ref_freq;
    cs2->rf        = cs1->rf;
    cs2->bw        = cs1->bw;
	
    return(0);
}


int	cs_subtract(const CS *cs1, CS *cs2) {
	/* Routine added by MAW 15/07/2011								*/
	/* Subtracts cyclic spectrum cs1 from cs2, with the result		*/
	/* returned in cs2												*/
	
	/* Check dimensions												*/
	if (cs1->nharm!=cs2->nharm || cs1->nchan!=cs2->nchan 
		|| cs1->npol!=cs2->npol) {
		printf("Error : cs_subtract : incompatible dimensions\n");
		exit(1);
	}
	
    int ih, ic, ip;
    for (ic=0; ic<cs1->nchan; ic++) {
        for (ip=0; ip<cs1->npol; ip++) {
            for (ih=0; ih<cs1->nharm; ih++) {
                fftwf_complex *d1 = get_cs(cs1, ih, ip, ic);
                fftwf_complex *d2 = get_cs(cs2, ih, ip, ic);
                *d2 -= *d1;
            }
        }
    }
    return(0);
}

int	cs_multiply(const CS *cs1, CS *cs2) {
	/* Routine added by MAW 15/07/2011								*/
	/* Multiplies cyclic spectrum cs2 by cs1, with the result		*/
	/* returned in cs2												*/
	
	/* Check dimensions												*/
	if (cs1->nharm!=cs2->nharm || cs1->nchan!=cs2->nchan 
		|| cs1->npol!=cs2->npol) {
		printf("Error : cs_multiply : incompatible dimensions\n");
		exit(1);
	}
	
    int ih, ic, ip;
    for (ic=0; ic<cs1->nchan; ic++) {
        for (ip=0; ip<cs1->npol; ip++) {
            for (ih=0; ih<cs1->nharm; ih++) {
                fftwf_complex *d1 = get_cs(cs1, ih, ip, ic);
                fftwf_complex *d2 = get_cs(cs2, ih, ip, ic);
                *d2 *= *d1;
            }
        }
    }
    return(0);
}

int cs_multiply_conj(const CS *cs1, CS *cs2) {
    /* Do conj(cs1) * cs2, put result in cs2 */
	/* Check dimensions												*/
	if (cs1->nharm!=cs2->nharm || cs1->nchan!=cs2->nchan 
		|| cs1->npol!=cs2->npol) {
		printf("Error : cs_multiply : incompatible dimensions\n");
		exit(1);
	}
	
    int ih, ic, ip;
    for (ic=0; ic<cs1->nchan; ic++) {
        for (ip=0; ip<cs1->npol; ip++) {
            for (ih=0; ih<cs1->nharm; ih++) {
                fftwf_complex *d1 = get_cs(cs1, ih, ip, ic);
                fftwf_complex *d2 = get_cs(cs2, ih, ip, ic);
                *d2 *= conj(*d1);
            }
        }
    }
    return(0);
}

int	cs_conjugate(CS *cs1) {
	/* Routine added by MAW 15/07/2011								*/
	/* Forms the complex conjugate of cs1, in place					*/
	
    int ih, ic, ip;
	fftwf_complex tmp = 0.0 + I * 0.0;
    for (ic=0; ic<cs1->nchan; ic++) {
        for (ip=0; ip<cs1->npol; ip++) {
            for (ih=0; ih<cs1->nharm; ih++) {
                fftwf_complex *d1 = get_cs(cs1, ih, ip, ic);
				tmp	= conj(*d1);
                *d1 = tmp;
            }
        }
    }
    return(0);
}

int	filter2cs(CS *cs1, const struct filter_freq *hf) {
	/* Routine added by MAW 15/07/2011								*/
	/* Fills the cyclic spectrum cs1 uniformly with the filter hf	*/
	
	/* Check dimensions												*/
	if (cs1->nchan!=hf->nchan) {
		printf("Error : filter2cs : incompatible dimensions\n");
		exit(1);
	}
	
    int ih, ic, ip;
    for (ic=0; ic<cs1->nchan; ic++) {
        for (ip=0; ip<cs1->npol; ip++) {
            for (ih=0; ih<cs1->nharm; ih++) {
                fftwf_complex *d1 = get_cs(cs1, ih, ip, ic);
                *d1 = hf->data[ic];
            }
        }
    }
    return(0);
}

int	profile2cs(CS *cs1, const struct profile_harm *s0) {
	/* Routine added by MAW 15/07/2011								*/
	/* Fills the cyclic spectrum cs1 uniformly with the profile s0	*/
	
	/* Check dimensions												*/
	if (cs1->nharm!=s0->nharm) {
		printf("Error : profile2cs : incompatible dimensions\n");
		exit(1);
	}
	
    int ih, ic, ip;
    for (ic=0; ic<cs1->nchan; ic++) {
        for (ip=0; ip<cs1->npol; ip++) {
            for (ih=0; ih<cs1->nharm; ih++) {
                fftwf_complex *d1 = get_cs(cs1, ih, ip, ic);
                *d1 = s0->data[ih];
            }
        }
    }
    return(0);
}


int rotate_phase_filter_freq(struct filter_freq *hf, const int rc) {
	/* Routine added by MAW 14/07/2011								*/
	/* Rotates the phase of H(freq) so that H(rc) is real			*/
	
	fftwf_complex phasor = conj(hf->data[rc]);
	float norm = creal(phasor * conj(phasor));
	phasor /= sqrt(norm);
	int ic;
	for (ic=0; ic<hf->nchan; ic++) {
		hf->data[ic] *= phasor;
	}
	
	return(0);
}

int rotate_phase_filter_time(struct filter_time *ht, const int rl) {
	/* Routine added by MAW 24/08/2011								*/
	/* Rotates the phase of h(lag) so that h(rl) is real			*/
	
	fftwf_complex phasor = conj(ht->data[rl]);
	float norm = creal(phasor * conj(phasor));
	phasor /= sqrt(norm);
	int ic;
	for (ic=0; ic<ht->nlag; ic++) {
		ht->data[ic] *= phasor;
	}
	
	return(0);
}


int parms2struct_freq(const double *x, struct filter_freq *hf, 
					  const int rchan) {
	/* Puts the 2*nchan-1 parameters x into the nchan complex		*/
	/* numbers making up hf, with the channel rchan having zero		*/
	/* imaginary component											*/
	/* Added by MAW 14/07/2011										*/
	
    int ic;
	const double *xtmp = x;
	for (ic=0; ic<rchan; ic++) {
		hf->data[ic] = xtmp[0] + I*xtmp[1];
		xtmp += 2;
	}
	hf->data[rchan] = xtmp[0] + I * 0.0;
	xtmp += 1;
	for (ic=rchan+1;ic<hf->nchan; ic++) {
		hf->data[ic] = xtmp[0] + I*xtmp[1];
		xtmp += 2;
	}
	
	return(0);
}

int struct2parms_freq(const struct filter_freq *hf, double *x, 
					  const int rchan) {
	/* Takes the 2*nchan-1 parameters x from the nchan complex		*/
	/* numbers making up hf, with the channel rchan having zero		*/
	/* imaginary component											*/
	/* Added by MAW 14/07/2011										*/
	
    int ic;
	double *xtmp = x;
	for (ic=0; ic<rchan; ic++) {
		xtmp[0] = (double)creal(hf->data[ic]);
		xtmp[1] = (double)cimag(hf->data[ic]);
		xtmp += 2;
	}
	xtmp[0] = (double)creal(hf->data[rchan]);
	xtmp += 1;
	for (ic=rchan+1;ic<hf->nchan; ic++) {
		xtmp[0] = (double)creal(hf->data[ic]);
		xtmp[1] = (double)cimag(hf->data[ic]);
		xtmp += 2;
	}
	
	return(0);
}


int parms2struct_time(const double *x, struct filter_time *ht, 
					  const int rlag) {
	/* Puts the 2*nlag-1 parameters x into the nlag complex			*/
	/* numbers making up ht, with ht(rlag) having zero				*/
	/* imaginary component											*/
	/* Added by MAW 24/08/2011										*/
	
    int ic;
	const double *xtmp = x;
	for (ic=0; ic<rlag; ic++) {
		ht->data[ic] = xtmp[0] + I*xtmp[1];
		xtmp += 2;
	}
	ht->data[rlag] = xtmp[0] + I * 0.0;
	xtmp += 1;
	for (ic=rlag+1;ic<ht->nlag; ic++) {
		ht->data[ic] = xtmp[0] + I*xtmp[1];
		xtmp += 2;
	}
	
	return(0);
}

int struct2parms_time(const struct filter_time *ht, double *x, 
					  const int rlag) {
	/* Takes the 2*nlag-1 parameters x from the nlag complex		*/
	/* numbers making up ht, with ht(rlag) having zero				*/
	/* imaginary component											*/
	/* Added by MAW 24/08/2011										*/
	
    int ic;
	double *xtmp = x;
	for (ic=0; ic<rlag; ic++) {
		xtmp[0] = (double)creal(ht->data[ic]);
		xtmp[1] = (double)cimag(ht->data[ic]);
		xtmp += 2;
	}
	xtmp[0] = (double)creal(ht->data[rlag]);
	xtmp += 1;
	for (ic=rlag+1;ic<ht->nlag; ic++) {
		xtmp[0] = (double)creal(ht->data[ic]);
		xtmp[1] = (double)cimag(ht->data[ic]);
		xtmp += 2;
	}
	
	return(0);
}


int cyclic_shear_cs(CS *d, double shear,
					const struct cyclic_work *w) {
	
	/* Added by MAW 16/7/2011. Modelled on PBD's cyclic_shift_cs	*/
	/* This version replaces the (integer) parameter "sign" with	*/
	/* the (double) parameter "shear", which is the amount by which	*/
	/* the cs should be sheared, in units of alpha. Thus we have	*/
	/*																*/
	/* FOR ALPHA > 0												*/
	/* shift = -0.5 to form CS(nu-alpha/2,alpha) from CS(nu,alpha)	*/
	/* shift = +1.0 to form CS(nu+alpha,  alpha) from CS(nu,alpha)	*/
	/*																*/
	/* FOR ALPHA < 0												*/
	/* shift = +0.5 to form CS(nu-alpha/2,alpha) from CS(nu,alpha)	*/
	/* shift = -1.0 to form CS(nu+alpha,  alpha) from CS(nu,alpha)	*/
	
    CC tmp_cc;
    tmp_cc.nharm = d->nharm;
    tmp_cc.npol  = d->npol;
    tmp_cc.nlag  = d->nchan;
    cyclic_alloc_cc(&tmp_cc);
	
    const double dtau   = 1.0/(d->bw*1e6);/* lag step, in seconds	*/
    const double dalpha = d->ref_freq;    /* harmonic step in Hz	*/
		
    /* Move to cc domain											*/
    cyclic_cs2cc(d, &tmp_cc, w);
	
    /* Multiply by shift function									*/
    int iharm, ilag, ipol;
    for (ilag=0; ilag<tmp_cc.nlag; ilag++) {
        for (iharm=0; iharm<tmp_cc.nharm; iharm++) {
            int lag = (ilag<=tmp_cc.nlag/2) ? ilag : ilag-tmp_cc.nlag;
            double phs = shear * (-2.0*M_PI) *(dalpha*(double)iharm) * 
						(dtau*(double)lag);
            fftwf_complex fac = (cos(phs)+I*sin(phs));
            for (ipol=0; ipol<tmp_cc.npol; ipol++) {
                fftwf_complex *dat = get_cc(&tmp_cc,iharm,ipol,ilag);
                *dat *= fac;
            }
        } 
    }
	
    /* Back to cs domain											*/
    cyclic_cc2cs(&tmp_cc, d, w);
	
	/* Free the storage and return									*/
    cyclic_free_cc(&tmp_cc);
    return(0);
}

int cyclic_fscrunch_cs(struct profile_harm *out, const CS *in) {
	/* Added by MAW 17/7/2011										*/
	/* Sums the input cyclic spectrum over radio frequency			*/
	/* and puts the result into the output harmonic profile			*/

    /* Only 1 pol for now */
    if (in->npol>1 || in->nharm!=out->nharm) {
		printf("cyclic_fscrunch_cs : incompatible dimensions\n");
		return(-1); 
	}
	
    int iharm, ichan, chan_min, chan_max;
    for (iharm=0; iharm<in->nharm; iharm++) {
		chan_limits_cs(&chan_min, &chan_max, iharm, in);
        out->data[iharm] = 0.0 + I * 0.0;
        for (ichan=chan_min; ichan<chan_max; ichan++) {
            const fftwf_complex *tmp = get_cs(in,iharm,0,ichan);
            out->data[iharm] += *tmp;
        }
    }
	
    return(0);
}


int	chan_limits_cs(int *chan_min, int *chan_max, const int iharm, 
				   const CS *cs) {
	
	/* Function added by MAW 17/7/2011								*/
	/* Returns the min and max channel numbers giving a valid CS	*/
	/* estimate at harmonic # iharm									*/
	/* Usage of these limits is as follows							*/
	/* for (ichan=chan_min; ichan<chan_max; ichan++) data are valid	*/
	
	if (iharm<0 || iharm>=cs->nharm) {
		printf("Error : chan_limits_cs : iharm out of range\n");
		exit(1);
	}
		
		
	double inv_aspect = cs->ref_freq * (double)cs->nchan ;
	inv_aspect *= (double)iharm / (cs->bw*1.e6);
	inv_aspect -= 1.0;
	inv_aspect /= 2.0;
	int ichan = (int)inv_aspect + 1;
	
	if (ichan>cs->nchan/2) { ichan = cs->nchan / 2; }

	*chan_min = ichan;
	*chan_max = cs->nchan - ichan;
	
	return(0);
	
}

int normalise_profile(struct profile_harm *s0) {
	
	/* Function added by MAW 22/07/2011								*/
	/* Normalises the reference pulse profile harmonics:|s0(1)| = 1 */

	int ih = 1;
	float normfac = 1./sqrt(creal(s0->data[ih] * conj(s0->data[ih])));
	
	for (ih=0; ih<s0->nharm; ih++) { s0->data[ih] *= normfac; }
		
	return(0);
}

int normalise_cs(CS *cs) {
	
	/* Function added by MAW 05/03/2012								*/
	/* Normalises the cyclic spectrum so that the rms signal power	*/
	/* at ih = 1 is unity											*/
	
	int ic, ip, ih = 1;
	float rms1 = rms_cs(cs,ih); 
	ih = cs->nharm-1;
	float rmsn = rms_cs(cs,ih); 
	float normfac = 1./sqrt(fabs(rms1*rms1 - rmsn*rmsn));; 
	
	for (ih=0; ih<cs->nharm; ih++) {
		for (ip=0; ip<cs->npol; ip++) {
			for (ic=0; ic<cs->nchan; ic++) {
				fftwf_complex *cs1 = get_cs(cs, ih, ip, ic);
				*cs1 *= normfac;
			}
		}
	}
	
	return(0);
}


int normalise_cs_old(CS *cs) {
	
	/* Function added by MAW 24/11/2011								*/
	/* Normalises the cyclic spectrum :	rms(|cs(1,freq)|) = 1		*/
	
	int ic, ip, ih = 1;
	float normfac = 1./rms_cs(cs,ih); 
	
	for (ih=0; ih<cs->nharm; ih++) {
		for (ip=0; ip<cs->npol; ip++) {
			for (ic=0; ic<cs->nchan; ic++) {
				fftwf_complex *cs1 = get_cs(cs, ih, ip, ic);
				*cs1 *= normfac;
			}
		}
	}
	
	return(0);
}


int	maximum_filter_time(const struct filter_time *ht) {
	/* Routine added by MAW 25/08/2011								*/
	/* Returns lag_max: |ht(lag_max)| is maximum					*/
	
	int il, lag_max=0;
	float ht_abs_sq=0., ht_abs_sq_max=0.;
	for (il=0; il<ht->nlag; il++) {
		ht_abs_sq = creal( ht->data[il] * conj( ht->data[il] ) );
		if (ht_abs_sq > ht_abs_sq_max) {
			ht_abs_sq_max = ht_abs_sq;
			lag_max	= il;
		}
	}
	
	return(lag_max);
}


int cyclic_padding(CS *d) {
	
	/* Added by MAW 17/11/2011.										*/
	/* Replaces the meaningless (unsampled) channels at the ends of	*/
	/* a cyclic spectrum with zeros.								*/
	
    int imin = 0;
	int imax = d->nchan;
	int ih, ip, ic;
	
	for (ih=0; ih<d->nharm; ih++) {
		chan_limits_cs(&imin, &imax, ih, d);
		for (ip=0; ip<d->npol; ip++) {
			for (ic=0; ic<imin; ic++) {
				fftwf_complex *d2 = get_cs(d, ih, ip, ic);
				*d2 = 0. + I * 0.;
			}
			for (ic=imax; ic<d->nchan; ic++) {
				fftwf_complex *d2 = get_cs(d, ih, ip, ic);
				*d2 =  0. + I * 0.;
			}
		}
	}
	
	
	return(0);
}


int phase_gradient(const CS *d, const struct profile_harm *s0) {
	
	/* Added by MAW 19/11/2011.										*/
	/* Estimates the mean filter phase gradient, given an input		*/
	/* cyclic spectrum and reference pulse profile harmonics.		*/
	/* The estimated phase gradient is then converted to a delay	*/
	/* and returned as an integer value, corresponding to the		*/
	/* non-zero element in the lag-space filter representation		*/
	
    int imin = 0, imax = d->nchan;
	int ih=1, ip=0, ic;
	int delay;
	float norm, phase_angle;
	
	/* Only one polarisation for now								*/
	
	chan_limits_cs(&imin, &imax, ih, d);
	fftwf_complex grad_sum = 0. + I * 0.;
	/* Sum over the cyclic spectrum values for ih = 1				*/
	for (ic=imin; ic<imax; ic++) {
		fftwf_complex *d2 = get_cs(d, ih, ip, ic);
		grad_sum += *d2;
	}

	/* Divide by the pulse-profile harmonic, to give the phasor		*/
	/* resulting from the mean phase gradient						*/
	grad_sum /= s0->data[ih];
	
	/* Normalise the phasor to unity								*/
	norm	 = sqrt((float)creal(grad_sum)*(float)creal(grad_sum) +
					(float)cimag(grad_sum)*(float)cimag(grad_sum));
	grad_sum /= norm;

	/* Solve for the corresponding phase angle in the range			*/
	/* -Pi < phase_angle < Pi										*/
	if (cimag(grad_sum) >= 0.) {
		phase_angle = acos(creal(grad_sum));
	}
	else {
		phase_angle	= -1. * acos(creal(grad_sum));
	}

	/* Express the phase angle in units of the pulsar's angular		*/
	/* rotation frequency, so that the result is a measure of delay	*/
	phase_angle /= -2. * M_PI * d->ref_freq;
	
	/* Now express this in units of the lag-space resolution		*/
	phase_angle *= 1.e6 * d->bw;
	
	/* And map this onto the actual lag array						*/
	if (phase_angle > d->nchan/2 ) {
		delay =	d->nchan/2;
	}
	else if (phase_angle < -(d->nchan/2) ) {
		delay =	d->nchan/2 + 1;
	}
	else if (phase_angle < 0.) {
		delay = (int)phase_angle + d->nchan - 1;
	}
	else {
		delay = (int)phase_angle;
	}


	if (delay < 0 || delay >= d->nchan) {
		printf("Error in phase_gradient : delay = %d\n", delay);
	}
	
	return(delay);
}


int filter_freq_renorm(struct filter_freq *hf) {
	/* Routine added by MAW 23/11/2011								*/
	/* Divides H(freq) by nchan.									*/
	/* Used by filter_time2freq so that it is the inverse			*/
	/* operation of filter_freq2time								*/
	
	int ic;
	for (ic=0; ic<hf->nchan; ic++) {
		hf->data[ic] /= (float)hf->nchan;
	}
	
	return(0);
}

int profile_harm_renorm(struct profile_harm *ph) {
	/* Routine added by MAW 23/11/2011								*/
	/* Divides ph(iharm) by nharm.									*/
	/* Used by profile_phase2harm so that it is the inverse			*/
	/* operation of profile_harm2phase								*/
	
	int ih;
	for (ih=0; ih<ph->nharm; ih++) {
		ph->data[ih] /= (float)ph->nharm;
	}
	
	return(0);
}

int cyclic_cc2cs_renorm(CS *cs) {
	
	/* Routine added by MAW 23/11/2011								*/
	/* Divides cs by nchan. Used by cyclic_cc2cs so that it is the	*/
	/* inverse operation of cyclic_cs2cc							*/
	
	int ih, ip, ic;
	for (ih=0; ih<cs->nharm; ih++) {
		for (ip=0; ip<cs->npol; ip++) {
			for (ic=0; ic<cs->nchan; ic++) {
				fftwf_complex *d1 = get_cs(cs, ih, ip, ic);
				*d1 /= (float)cs->nchan;
				/* *d1 /= 1.;										*/
			}
		}
	}
	
	return(0);
}


int cyclic_ps2cs_renorm(CS *cs) {
	
	/* Routine added by MAW 23/11/2011								*/
	/* Divides cs by nharm. Used by cyclic_ps2cs so that it is the	*/
	/* inverse operation of cyclic_cs2ps (not yet implemented!)		*/
	
	int ih, ip, ic;
	for (ih=0; ih<cs->nharm; ih++) {
		for (ip=0; ip<cs->npol; ip++) {
			for (ic=0; ic<cs->nchan; ic++) {
				fftwf_complex *d1 = get_cs(cs, ih, ip, ic);
				*d1 /= (float)cs->nharm;
				/* *d1 /= 1.;										*/
			}
		}
	}
	
	return(0);
}

int match_two_filters(const struct filter_freq *hf1,
							struct filter_freq *hf2) {
	
/* Routine added by MAW 30/11/2011									*/
/* Determines the phasor z such that z * H2(freq) is the best match	*/
/* to the filter H1(freq), and modifies H2 accordingly				*/
/* Also normalises H2 so that the r.m.s. value is unity				*/
	
	if (hf1->nchan != hf2->nchan) {
		printf("Error: match_two_filters: incompatible dimensions\n");
		exit(1);
	}
	
	
	/* Determine the phasor which minimises the						*/
	/* difference between H1(freq) and H2(freq)						*/
	int ic;	
	float znorm;
	fftwf_complex z  = 0. + I * 0.;
	fftwf_complex z2 = 0. + I * 0.;
	for (ic=0; ic<hf1->nchan; ic++) {	
		z	 += hf1->data[ic] * conj(hf2->data[ic]);
		z2   += hf2->data[ic] * conj(hf2->data[ic]);						
	}
	znorm = sqrt(creal(z * conj(z)));
	z    /= znorm;
	znorm = sqrt((float)hf1->nchan / creal(z2));
	z    *= znorm;
		
	/* Scale H2(freq) by this factor and return						*/
	for (ic=0; ic<hf2->nchan; ic++) { hf2->data[ic] *= z; }
	 
	return(0);
}


int cyclic_variance(float *variance, float *vpoints, const CS *d) {
	
	/* Routine added by MAW 16/01/2012								*/
	/* Estimates the variance of the cyclic spectrum, using the		*/
	/* highest available harmonic of the pulse (modulation) freq.	*/
	/* Also computes the number of valid points contributing to		*/
	/* the estimate of the merit function.							*/
	/* Routine assumes that the noise is white.						*/
	
	/* Only 1 pol for now											*/
    if (d->npol>1) {
		printf("cyclic_variance : npol > 1 in data\n");
		return(-1); 
	}
	
	float valid_points = 0., cyclic_var = 0.;
	
	/* Estimate the variance using the highest harmonic data		*/
	int ip = 0, ih = d->nharm - 1;
	int minchan=0, maxchan=0;
	chan_limits_cs(&minchan, &maxchan, ih, d);
	
	int ic;
	for (ic=minchan; ic<maxchan; ic++) {
		fftwf_complex *d1 = get_cs(d, ih, ip, ic);
		cyclic_var += (float)creal(*d1 * conj(*d1));
		valid_points += 1.;
	}
	cyclic_var /= valid_points;
	
	for (ih=1; ih<d->nharm-1; ih++) {
		chan_limits_cs(&minchan, &maxchan, ih, d);
		valid_points += (float)maxchan;
		valid_points -= (float)minchan;
	}
	
	*vpoints  = 2. * valid_points;
	*variance = cyclic_var;

	return(0);
}



