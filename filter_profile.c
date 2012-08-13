/* filter_profile.c													*/
/*																	*/
/* Original code by P. Demorest, December 2010						*/
/* Modified by M. Walker July 2011 - March 2012						*/
/*																	*/
/* Given a measured periodic spectrum, and a reference intrinsic	*/
/* pulse profile, this code uses the NLOPT optimization library to	*/
/* determine the best fit (in a least-squares sense) impulse		*/
/* response function. And for this IRF, the intrinsic pulse			*/
/* profile implied by the measured periodic spectrum is also		*/
/* determined.														*/
/*																	*/
/* Filter optimisation is available in lag-space (default),			*/
/* or in frequency-space (-f option)								*/
/*																	*/
/* If no reference pulse profile is given, the code will determine	*/
/* a pulse profile from the data, assuming that the IRF	is a		*/
/* delta-function at lag=0.											*/
/*																	*/

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
#include <nlopt.h>
#include <string.h>

#include "cyclic_utils.h"
#include "cyclic_fileio.h"
#include "model_cyclic.h"
#include "merit_functions.h"

void usage() {printf("filter_profile [Options] filename\n");
    printf("Options:\n");
    printf("  -v verbose\n");
    printf("  -t nthread Number of FFTW threads\n");
    printf("  -S isub    Start at subint # isub\n");
    printf("  -N nsub    Process nsub total subints\n");
    printf("  -I nchan   Ignore this many channels at each band edge\n");
    printf("  -f frequency-space optimisation\n");
    printf("     (default is lag-space optimisation)\n");
	printf("  -i initialise (no filter optimisation)\n");
    printf("  -R fname  Use the Reference pulse profile in fname\n");
    printf("Must invoke with -i or -R options\n");
}

/* Integers which can be accessed by various functions				*/
int verbose = 0;
int sample_ncalls = 0;


/* Catch sigint														*/
int run=1;
void cc(int sig) { run=0; }

int main(int argc, char *argv[]) {
    int isub=1, opt=0, opcheck=0, do_optimisation=0, lagspace=1;
    int nchan_ignore=0;
    int nsub_proc=0;
    int fft_threads = 2;
	extern int verbose;
	extern int sample_ncalls;
	char *ref_prof="";
    while ((opt=getopt(argc,argv,"fivR:t:S:N:I:"))!=-1) {
        switch (opt) {
            case 'f':
				lagspace--;
                break;
            case 'i':
				opcheck++;
                break;
            case 'v':
                verbose++;
                break;
            case 'R':
                ref_prof = optarg;
				do_optimisation++;
				opcheck++;
                break;
            case 't':
                fft_threads = atoi(optarg);
                break;
            case 'S':
                isub = atoi(optarg);
                break;
            case 'N':
                nsub_proc = atoi(optarg);
                break;
            case 'I':
                nchan_ignore = atoi(optarg);
                break;
        }
    }

    if (optind==argc || opcheck!=1) { usage(); exit(1); }
	
	int ic, ih, j, is, rv, nspec=1;
	
    /* Open generic datafile 		            					*/
    struct cyclic_file cf;
    cf.err_status = 0; // Needs to be initialized to 0
    cyclic_file_open(&cf, argv[optind]);
    cyclic_file_error_check_fatal(&cf);

    /* Get basic dimensions											*/
    struct cyclic_work w;
    cyclic_load_params(&cf, &w, &nspec);
    cyclic_file_error_check_fatal(&cf);
    
    int nsub_max = nspec;
    if (nsub_proc > 0) {
        nsub_max = isub + nsub_proc - 1;
        if (nsub_max > nspec) 
            nsub_max = nspec;
    }
	if (verbose) {
       printf("Read nphase = %d, npol = %d, nchan = %d, nspec = %d\n", 
			   w.nphase, w.npol, w.nchan, nspec);
       printf("Processing %d total spectra, starting at subint %d\n",
               nsub_max - isub + 1, isub);
       printf("Ignoring %d channels at each band edge\n", nchan_ignore);
	   fflush(stdout);
    }

    if (nchan_ignore*2 >= w.nchan) {
        fprintf(stderr, "nchan_ingore=%d is too large (nchan=%d)\n",
                nchan_ignore, w.nchan);
        exit(1);
    }
    
    int orig_npol = w.npol;
    w.npol = 1;

    int orig_nchan = w.nchan;
    w.nchan = w.nchan - 2*nchan_ignore;

    /* Initialise FFTs												*/
    fftwf_init_threads();
    fftwf_plan_with_nthreads(fft_threads);
    if (verbose) { printf("Planning FFTs\n"); fflush(stdout); }
#define WF "cyclic_wisdom.dat"
    FILE *wf = fopen(WF,"r");
    if (wf!=NULL) { fftwf_import_wisdom_from_file(wf); fclose(wf); }
    rv = cyclic_init_ffts(&w);
    if (rv) {
        fprintf(stderr, "Error planning ffts (rv=%d)\n", rv);
        exit(1);
    }
    wf = fopen(WF,"w");
    if (wf!=NULL) { fftwf_export_wisdom_to_file(wf); fclose(wf); }

	
    /* Allocate some stuff											*/
    struct periodic_spectrum raw, rawraw;
    struct cyclic_spectrum cs;
    struct filter_freq hf, hf_prev;
	struct filter_time ht;
    struct profile_phase pp, pp_ref, pp_int;
    struct profile_harm  ph, ph_ref;

    raw.nphase = pp.nphase = pp_ref.nphase = pp_int.nphase = w.nphase;
    raw.nchan = cs.nchan = hf.nchan = hf_prev.nchan = w.nchan;
    cs.nharm = ph.nharm = ph_ref.nharm = w.nharm;
	ht.nlag = w.nlag;
    raw.npol = 1;
    cs.npol = 1;

    rawraw.nphase = raw.nphase;
    rawraw.nchan = orig_nchan;
    rawraw.npol = orig_npol;

    cyclic_alloc_ps(&raw);
    cyclic_alloc_ps(&rawraw);
    cyclic_alloc_cs(&cs);
    filter_alloc_freq(&hf);
    filter_alloc_freq(&hf_prev);
	filter_alloc_time(&ht);
    profile_alloc_phase(&pp);
    profile_alloc_phase(&pp_ref);
    profile_alloc_phase(&pp_int);
    profile_alloc_harm(&ph);
    profile_alloc_harm(&ph_ref);

	/* Initialise arrays for dynamic spectrum and optimised filters	*/
	float dynamic_spectrum[nspec][w.nchan];
	fftwf_complex optimised_filters[nspec][w.nchan];
	for (is=0; is<nspec; is++) {
		for (ic=0; ic<w.nchan;ic++) {
			dynamic_spectrum[is][ic]  = 0.;
			optimised_filters[is][ic] = 0. + I * 0.;
		}
	}
	
	/* Initialise arrays to record optimisation stats				*/
	float minima[nspec];
	int outcome[nspec];
	int	callno[nspec];
	for (is=0; is<nspec; is++) {
		minima[is]  = 0.0;
		outcome[is] = 0;
		callno[is]  = 0;
	}
	
	/* Initialise the intrinsic pulse profile to zero				*/
	for (j=0; j<pp_int.nphase; j++) { pp_int.data[j] = 0.0; }
		
	/* Initialise the "previous" filter coefficients to unity		*/
	for (ic=0; ic<hf_prev.nchan; ic++) {
		hf_prev.data[ic] = 1.0 + I * 0.0;
	}
	
	if (do_optimisation) {
		/* Read in the reference pulse profile						*/
		read_profile(ref_prof, &pp_ref);
		/* Convert reference profile to harmonics					*/
		profile_phase2harm(&pp_ref, &ph_ref, &w);
		/* Normalise : |ph_ref(1)| = 1.								*/
		normalise_profile(&ph_ref);
		/* Ensure the zero-frequency term is zero					*/
		ph_ref.data[0] = 0.0 + I * 0.0;
	}

	
	/* Loop over all subintegrations								*/
	int noptimised=0;
	while (isub <= nsub_max) {
	
		if (verbose) {
			printf("Subintegration %d of %d\n", isub, nspec);
		}
		
		/* Load data												*/
		rawraw.npol = orig_npol;
		cyclic_load_ps(&cf, &rawraw, isub);
		cyclic_file_error_check_fatal(&cf);
				
		/* Add polarisations without calibration					*/
		cyclic_pscrunch_ps(&rawraw, 1.0, 1.0);
		/* Only one polarisation from this point in loop			*/

        /* Trim requested number of edge channels                   */
        cyclic_remove_edge_chans(&rawraw, &raw, nchan_ignore);
		
		/* Convert input data to cyclic spectrum					*/
		cyclic_ps2cs(&raw, &cs, &w);
		
		/* Normalise the data to give unit rms signal power in the	*/
		/* cyclic spectrum at modulation frequency = pulsar			*/
		/* rotation frequency										*/
		normalise_cs(&cs);
		
		/* Replace the invalid (unsampled) end channels of the 		*/
		/* cyclic spectrum with zeros.								*/
		cyclic_padding(&cs);
		
		/* The zero modulation-frequency component of the cyclic	*/
		/* spectrum gives us the dynamic spectrum -> output			*/
		for (ic=0; ic<w.nchan; ic++) {
			fftwf_complex *d1 = get_cs(&cs, 0, 0, ic);
			dynamic_spectrum[isub-1][ic]  = creal(*d1);
		}
		
		/* If we're initialising (opt = -i) then set H(freq)=1		*/
		int delay=0;
		if (!do_optimisation) {
			for (ic=0; ic<hf.nchan; ic++) {
				hf.data[ic] = 1.0 + I * 0.0;
			}
		}
		else if (!noptimised) {
		/* Otherwise if this is the first sample then we estimate	*/
		/* the phase gradient in H(freq) and set h(lag)	to be a		*/
		/* delta-function centred on the corresponding lag			*/
			delay = phase_gradient(&cs, &ph_ref);
			printf("Initial filter: delta-function at delay = %d\n",
				   delay);
			for (ic=0; ic<ht.nlag; ic++) {
				ht.data[ic] = 0.0 + I * 0.0;
			}
			ht.data[delay] = (float)ht.nlag + I * 0.0;
			/* Convert model filter to frequency-space				*/
			filter_time2freq(&ht, &hf, &w);						
		}
		
		/* Convert current filter to lag-space						*/
		filter_freq2time(&hf, &ht, &w);
		
		int rindex = 0;
		/* Find the index for which the filter is maximised			*/
		if (lagspace) { rindex = maximum_filter_time(&ht); }
		else          { rindex = maximum_cs1(&cs); }
		
		if (verbose) {
			printf("Real-valued filter coefficient at index = %d\n",
				   rindex);  
			fflush(stdout);
		}

		if (do_optimisation) {
			/* Rotate filter phases so that rindex is real-valued	*/
			if (lagspace) { rotate_phase_filter_time(&ht, rindex); }
			else {          rotate_phase_filter_freq(&hf, rindex); }

			/* Enforce ph_ref as the reference profile				*/
			for (ih=0; ih<ph.nharm; ih++) {
				ph.data[ih] = ph_ref.data[ih];
			}

			/* Fill in data struct for NLOPT						*/
			struct cyclic_data cdata;
			cdata.cs     = &cs;
			cdata.s0     = &ph;
			cdata.w	     = &w;
			cdata.rindex = &rindex;

			sample_ncalls = 0;
			
			/* Set up NLOPT minimizer								*/
			int dim0;				/* no. of free parameters		*/
			if (lagspace) { dim0 = 2*w.nlag-1; }
			else          { dim0 = 2*w.nchan-1;}
			
			const int dim = dim0;
			nlopt_opt op;
			op = nlopt_create(NLOPT_LD_LBFGS, dim);
			
			/* Other NLopt algorithms which use gradients			*/
			/* You shouldn't normally be using these as they		*/
			/* don't perform as well as LBFGS						*/
			/* op = nlopt_create(NLOPT_LD_VAR1, dim);				*/
			/* op = nlopt_create(NLOPT_LD_VAR2, dim);				*/
			/* op = nlopt_create(NLOPT_LD_MMA, dim);				*/
			/*op=nlopt_create(NLOPT_LD_TNEWTON_PRECOND_RESTART,dim);*/
			/* op = nlopt_create(NLOPT_LD_SLSQP, dim);				*/
						
			if (lagspace) {
				nlopt_set_min_objective(op, cyclic_merit_nlopt_lag,
										&cdata);
				if (verbose) {printf("Lag-space  optimisation\n");}		
			}
			else {
				nlopt_set_min_objective(op, cyclic_merit_nlopt_freq,
										&cdata);
				if (verbose) {printf("Freq-space optimisation\n");}				
			}

			if (verbose) {
				printf("Number of fit parameters = %d\n", dim);
			}
			
			float cs_var = 0., merit_samples = 0., dof = 0.;
			/* Determine the noise level in the data				*/
			cyclic_variance(&cs_var, &merit_samples, &cs);
			dof = merit_samples - (float)dim0 - (float)w.nphase;
			
			if (verbose) {
				printf("Variance of data         = %.5e\n",cs_var);
				printf("Number of samples        = %.5e\n",
					   merit_samples);
				printf("Degrees of Freedom       = %.5e\n",dof);
				printf("Expected Minimum Demerit = %.5e\n",
					   dof * cs_var);
			}
			
						
			/* Assert the stopping criterion						*/
			double tolerance = 1.e-1 / (double)dof;
			/* Fractional tolerance on the merit function			*/
			nlopt_set_ftol_rel(op, tolerance);
			if (verbose) {
				printf("NLOPT: set FTOL-rel      = %.5e\n",tolerance);
			}
			
			/* Set initial values of parameters	to current filter	*/
			double *x = (double *)malloc(sizeof(double) * dim);			
			
			if (x==NULL) {
				fprintf(stderr, "malloc(x) : insufficient space\n");
				exit(1);
			}
			
			if (lagspace) struct2parms_time(&ht, x, rindex);
			else          struct2parms_freq(&hf, x, rindex);			
			
			/* Run NLOPT optimization								*/
			double min;
			outcome[isub-1] = nlopt_optimize(op, x, &min);
			minima[isub-1]  = min;
			callno[isub-1]  = sample_ncalls;

			/* Construct H(freq) from the optimum values of "x"		*/
			if (lagspace) {
				parms2struct_time(x, &ht, rindex);
				filter_time2freq(&ht, &hf, &w);
			}
			else parms2struct_freq(x, &hf, rindex);
						
			/* Free-up allocations									*/
			free(x);
			nlopt_destroy(op);
		}
		
		/* Rotate the phases of the current filter so that they		*/
		/* match the previous time-step as closely as possible		*/
		/* And normalise the r.m.s. amplitude to unity				*/
		match_two_filters(&hf_prev, &hf);
		
		/* Put current optimised filter(freq) into output array		*/
		for (ic=0; ic<hf.nchan; ic++) {
			optimised_filters[isub-1][ic]=hf.data[ic];
		}
		
		/* Get optimised profile, given filter,	for this sub-int	*/
		int iprof = 0;
		iprof = optimise_profile(&ph, &cs, &hf, &w);
		if (iprof!=0) {
			fprintf(stderr, "Error in optimise_profile.\n");
			exit(1);
		}
		
        /* Convert profile(harmonic) to profile(phase)				*/
        ph.data[0] = 0.0 + I * 0.0;
        profile_harm2phase(&ph, &pp, &w);
				
		
		if (outcome[isub-1]<0 && verbose) 
            printf("NLOPT failed (code=%d)\n",outcome[isub-1]);
			
		/* Add profile for subint to the intrinsic estimate			*/
		for (j=0; j<pp_int.nphase; j++) { 
			pp_int.data[j] += pp.data[j];
		}
		/* Update the "previous" filter, used to determine the		*/
		/* complex scaling applied to the next optimised filter		*/
		for (j=0; j<hf.nchan; j++) { 
			hf_prev.data[j] = hf.data[j];
		}
			
		isub++;
		noptimised++;

	}

	char *inputptr = argv[optind];
	char dotchar = '.';
	char *dotpos = strrchr(inputptr,dotchar);
	size_t fnamelength = dotpos - inputptr + 1 ;
	char *outputptr = inputptr;	
	strcpy(outputptr,inputptr);
		
	/* Output intrinsic profile (= scattered profile if opt = -i)	*/
	outputptr[fnamelength] = '\0';
	strcat(outputptr,"profile.txt");
	write_profile(outputptr, &pp_int);
		
	/* Output the optimised filters									*/
	outputptr[fnamelength] = '\0';
	strcat(outputptr,"filters.txt");
	FILE *fpointer = fopen(outputptr, "w");
	fprintf(fpointer, "%d\n", nspec);
	fprintf(fpointer, "%d\n", w.nchan);
	for (is=0; is<nspec; is++) {
		for (ic=0; ic<w.nchan; ic++) {
			fprintf(fpointer, "%.7e %.7e\n",
					(float)creal(optimised_filters[is][ic]),
					(float)cimag(optimised_filters[is][ic]));
		}
	}
    fprintf(fpointer,"\n\n");
    fclose(fpointer);
	
	/* Output the dynamic spectrum									*/
	outputptr[fnamelength] = '\0';
	strcat(outputptr,"dynamic_spectrum.txt");
	fpointer = fopen(outputptr, "w");
	fprintf(fpointer, "%d\n", nspec);
	fprintf(fpointer, "%d\n", w.nchan);
	for (is=0; is<nspec; is++) {
		for (ic=0; ic<w.nchan; ic++) {
			fprintf(fpointer, "%.7e\n", dynamic_spectrum[is][ic]);
		}
	}
    fprintf(fpointer,"\n\n");
    fclose(fpointer);
	
	/* Output the optimisation stats								*/
	outputptr[fnamelength] = '\0';
	strcat(outputptr,"optimisation_stats.txt");
	fpointer = fopen(outputptr, "w");
	fprintf(fpointer, "%d\n", nspec);
	for (is=0; is<nspec; is++) {
		fprintf(fpointer, "%d %d %d %.9e\n", is+1, outcome[is],
				callno[is], minima[is]);
	}
    fprintf(fpointer,"\n\n");
    fclose(fpointer);
	
	/* Free-up the structs											*/
    cyclic_free_ps(&raw);
	cyclic_free_cs(&cs);
    filter_free_freq(&hf);
    filter_free_freq(&hf_prev);
    filter_free_time(&ht);
    profile_free_phase(&pp);
    profile_free_harm(&ph);
    profile_free_phase(&pp_ref);
    profile_free_harm(&ph_ref);
    profile_free_phase(&pp_int);
	cyclic_free_ffts(&w);
	
    exit(0);
}

