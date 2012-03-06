README file for Cyclic Modelling code 

This file should be read in conjunction with the paper
"Cyclic spectral analysis of radio pulsars"
P. Demorest 2011, Monthly Notices of the Royal Astronomical Society 416, 2821

The major part of this code was written by Paul Demorest in December 2010.
Subsequent contributions were made by Mark Walker over the period July 2011 - March 2012.

The main body of the code is least-squares optimisation software for fitting filter functions and intrinsic pulse profiles to measured cyclic spectra (or periodic spectra). The compiled optimiser is named filter_profile. Execute the makefile to compile.

Suppose you have a fits file (YourDataFile.fits) in which your periodic spectra are stored. To run the optimisation code for the first time, execute the command

filter_profile -v -i YourDataFile.fits

The -v option is the verbose option, so you can see what's going on.

The -i option is the initialise option. With this option the code doesn't attempt to solve for the filters. All it does is to construct the mean pulse profile over all the data you supply in the file YourDataFile.fits (this file may include more than one periodic spectrum). The mean pulse profile is written to the file YourDataFile.profile.txt and can be used as input to the next stage. With the initialise option, the output mean pulse profile is just the scattered profile.

Once you have a pulse profile as a reference, the software will generate a least-squares estimate of the complex filter function for each periodic spectrum you supply. Following initialisation, as above, you can execute

filter_profile -v -R YourDataFile.profile.txt YourDataFile.fits

which specifies the scattered pulse profile as the reference pulse profile. When invoked in this way the code uses the NLopt library
http://ab-initio.mit.edu/wiki/index.php/NLopt
to form a least-squares fit of the filter function to each of your cyclic spectra, given a reference pulse profile (which in this example is just the scattered pulse profile). These optimised filters are written to the file YourDataFile.filters.txt.

By default the optimisation is undertaken in lag-space, i.e. it is the impulse response function which is optimised. But by using the -f option it is possible to optimise the frequency-space coefficients of the filter. In the tests undertaken to date, lag-space optimisation usually gives slightly better results (i.e. lower residuals), and does so slightly faster, so you should not normally employ frequency-space optimisation.

Once the solver has determined an optimum filter, with respect to the input reference pulse profile, for a given cyclic spectrum, it uses that filter to form a least-squares estimate of the intrinsic pulse profile for each cyclic spectrum. This step is trivial and does not utilise NLopt. The average of these optimised pulse profiles, taken over all the input cyclic spectra, is written out, as described above, to YourDataFile.profile.txt.

When using the -i (initialise) option, the filter function is taken to be unity, so the output mean pulse profile is the scattered pulse profile. But when a reference pulse profile is input, so that the software generates least-squares estimates of the filter functions, the output mean pulse profile is a much better approximation to the intrinsic pulse profile. To the extent that the model filter correctly captures the structure of the actual impulse response function, the resulting mean pulse profile will correctly describe the intrinsic pulse profile.

Because the filter and profile models are optimised separately, a single execution of the code will not generate the best possible models if you use the scattered pulse profile as the reference profile. Rather you should go through several more iterations of 

filter_profile -v -R YourDataFile.profile.txt YourDataFile.fits

using, at each iteration, the YourDataFile.profile.txt generated at the previous iteration. These iterations can, of course, be skipped if you've already used the optimiser to determine a good model  of the intrinsic pulse profile and you now want to apply that model to some new data.

Although the separation of profile and filter optimisations does require one to approach the best possible models iteratively, the use of a reference pulse profile enforces a single timing fiducial on the resulting models.

One final note: the number of threads is currently set (in filter_profile.c) to nthreads=2, which might not be the best choice for your machine. 


M. Walker, Manly Astrophysics, 5th March 2012
