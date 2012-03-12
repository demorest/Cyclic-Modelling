/* cyclic_fileio_psrchive.C
 * P. Demorest, March 2012.
 * File I/O routines using PSRCHIVE.
 */
#include <fstream>
#include <iostream>

#include "Pulsar/Archive.h"
#include "Pulsar/Integration.h"
#include "Pulsar/Profile.h"
#include "MJD.h"

extern "C" {
#include "cyclic_fileio.h"
#include "cyclic_utils.h"
}

using namespace Pulsar;

void cyclic_file_error_check_fatal(struct cyclic_file *cf) {
    if (cf->err_status) {
        // TODO how to access the PSRCHIVE error stack??
        fprintf(stderr, "cyclic_file error\n");
        exit(1);
    }
}

int cyclic_file_open(struct cyclic_file *cf, const char *fname) try {

    // Note, need to be careful about Reference::To objects 
    // going out-of-scope prematurely in all this... 

    Archive *arch = Archive::load(fname);
    cf->file_ptr = (void *)arch;

    return(cf->err_status);
}
catch (Error &err) {
    // How to store the error message for future use...
    cf->err_status = 1;
    return (cf->err_status);
}


/* Load dimension params from a Pulsar::Archive                     */
int cyclic_load_params(struct cyclic_file *cf, struct cyclic_work *w,
					   int *nspec) try {

    Archive *arch = reinterpret_cast<Archive *>(cf->file_ptr);

    w->nphase = arch->get_nbin();
    w->npol   = arch->get_npol();
    w->nchan  = arch->get_nchan();
    *nspec    = arch->get_nsubint();

    w->nlag   = 0;
    w->nharm  = 0;

    return(cf->err_status);
}
catch (Error &err) {
    cf->err_status = 1;
    return (cf->err_status);
}

/* Load one periodic spectrum from a Pulsar::Archive 
 * Space should already be allocated.
 * idx is 1-offset following cfitsio convention.
 */
int cyclic_load_ps(struct cyclic_file *cf, PS *d, int idx) try {

    Archive *arch = reinterpret_cast<Archive *>(cf->file_ptr);

    // This integration object will get destroyed after this function
    // exits.
    Reference::To<Integration> subint = arch->get_Integration(idx-1);

    // We need to fill the psrchive-format data into the expected
    // array.
    const unsigned nchan = subint->get_nchan();
    const unsigned npol = subint->get_npol();
    const unsigned nbin = subint->get_nbin();
    for (unsigned ichan=0; ichan<nchan; ichan++) {
        for (unsigned ipol=0; ipol<npol; ipol++) {
            Reference::To<Profile> prof = subint->get_Profile(ipol,ichan);
            for (unsigned ibin=0; ibin<nbin; ibin++) {
                float *dat = get_ps(d, ibin, ipol, ichan);
                *dat = prof->get_amps()[ibin];
            }
        }
    }

    // Load timestamp, by psrchive convention this should always
    // correspond to a time when phase=0.
    MJD epoch = subint->get_epoch();
    d->imjd = epoch.intday();
    d->fmjd = epoch.fracday();
    d->ref_phase = 0.0;
    d->ref_freq = 1.0/subint->get_folding_period();

    // Freq, BW params
    d->rf = subint->get_centre_frequency();
    d->bw = subint->get_bandwidth();

    return(cf->err_status);
}
catch (Error &err) {
    cf->err_status = 1;
    return (cf->err_status);
}

