; NAME:
;       Fitline MCMC
;
; AUTHOR:
;       Ryan Cybulski
;
; PURPOSE: 
;        Identify spectral lines (in emission), fit them with a Gaussian,
;        and measure the integrated line flux of the identified line, all 
;        with a Markov Chain Monte Carlo (MCMC) framework. The Gaussian 
;        has the form y = a[0]*exp(-z^2/2) + a[3], where z=(x-a[1])/a[2],
;        and so "a" is a 4-element array consisting of:
;        a[0]: normalization
;        a[1]: central frequency
;        a[2]: standard deviation
;        a[3]: vertical offset
;
;        The code is designed to read in a Redshift Search Receiver (RSR) 
;        spectrum, but in principle it could work for any ASCII files 
;        having the format (with no header lines) of: 
;        Freq(GHz)     T_A^*(K)
;        
;
; CALLING SEQUENCE:
; IDL> fitline_mcmc, specfile, f0, [fwindow=fwindow, g0=g0, $
;      normrange=normrange, frange=frange, sigrange=sigrange, $
;      offsetrange=offsetrange, nsteps=nsteps, nburnin=nburnin, 
;      savefigs=savefigs, specplotname=specplotname, $
;      paramplotname=paramplotname, stepplotname=steppolotname, $
;      showsteps=showsteps, movie=movie, stepsizes=stepsizes, $
;      base_poly_ord=base_poly_ord, baseplot=baseplot, $
;      baseplotname=baseplotname, local_baseline=local_baseline, $
;      nobasesubtract=nobasesubtract, $
;      first_guess_gaussian=first_guess_gaussian, $
;      savsetname=savsetname, save_mcmc_steps=save_mcmc_steps, $
;      noplots=noplots, gain=gain]
;
; INPUTS:
;       specfile - full path to the ASCII file with the spectrum
;
;       f0 - central frequency (GHz). The might be your expected central 
;            freq. of the line, but it need not be if your goal is to 
;            sample a wide range of parameter space. It can just provide a 
;            starting point for the MCMC chain, if your "frange" is set  
;            to span a wide range of frequencies.
;
; OPTIONAL INPUTS:
;
;       fwindow - 2-element array with: [f_lower, f_upper], where:
;        f_lower is the lower lim of the frequency range being examined (GHz)
;        f_upper is the upper lim of the frequency range being examined (GHz)
;        default = [f0 - 1.0, f0 + 1.0]
;        *NOTE* that this is different from "frange" (see below), which forces 
;        the Gaussian to be within a user-specified frequency range near f0. 
;        Even if one is only fitting a Gaussian on a very limited frequency 
;        range, it's necessary to have a larger window (i.e. around +/- 1 GHz) 
;        in which to measure the RMS of the spectrum.
;
;            g0 - initial guess at Gaussian parameter array with:
;             g0[0] = normalization (mK) - default=0.3
;             g0[1] = central frequency (GHz) - default=f0
;             g0[2] = standard deviation (GHz) - default=0.015
;             g0[3] = DC offset (mK) - default=0.01
;        NOTE: if these parameters are initially given values that lie outside 
;        of the ranges restricted elsewhere (see frange, normrange, etc.), then 
;        the appropriate g0 entries will be automatically set as the mean of the 
;        specified restricted range
;
;     normrange - 2-element array giving the limits one wishes to apply to 
;        restrict the MCMC chains in their normalization - default=[0.0, 100.0]
;
;        frange - 2-element array giving the limits one wishes to apply to 
;        restrict the MCMC chains in their central frequency - default=[f0-0.08, f0+0.08]
;
;      sigrange - 2-element array giving the limits one wishes to apply to 
;        restrict the MCMC chains in their std. deviation - default=[0.0, 0.15]
;
;   offsetrange - 2-element array giving the limits one wishes to apply to 
;        restrict the MCMC chains in their DC offset - default=[-100.0, 100.0]
;
;        nsteps - number of MCMC steps to take - default = 100000d
;
;       nburnin - number of initial "burn-in" steps to take before officially 
;        recording the range of parameter space being explored as reflecting the 
;        statistical distribution of the parameters themselves - default = 40000d
;        (NOTE that the length of the burn-in period will be inherently dependent 
;        on the size of parameter space one is exploring, and the step sizes 
;        being taken [see "paramsigma" below]. The best way to be absolutely sure 
;        that you have an appropriate number of MCMC steps and burn-in period is 
;        to look at the chains of parameter values themselves, which you can do 
;        by setting the /showsteps flag [see below])
;
;      savefigs - save EPS figures of the spectrum and best-fit Gaussian, as 
;        well as the histograms of statistical distributions of parameter values 
;        and MCMC chains (NOTE, this will suppress graphics from being displayed, 
;        except for "movie" plots)
;
;  specplotname - name of plot showing the spectrum and best-fit Gaussian
;
; paramplotname - name of plot showing the histograms of parameter values
;
;  stepplotname - name of plot showing the chains of MCMC parameter values
;
;      showsteps - display the chain of MCMC parameter values explored by the code
;
;        movie - if set, it will display a plot of the spectrum and show periodic snapshots, 
;        every 1000 MCMC steps, of the current Gaussian being explored by the MCMC chain
;
;     stepsizes - four-element array containing multiplicative values to modify the 
;        relative step sizes for MCMC to take for each parameter. NOTE that the default
;        step sizes being used have been calibrated to the application of finding 
;        the highest S/N line in a frequency window of ~ +/- 1 GHz in an RSR spectrum, 
;        and for measuring the relevant parameters after finding that line. Therefore,
;        if one wishes to search for a line in the full 74-111 GHz range, then 
;        it might be a good idea to increase the step size taken for the frequency 
;        parameter, i.e. by specifying stepsizes=[1.0, 20.0, 1.0, 1.0] in addition 
;        to setting fwindow=[74.0, 111.0] and frange=[74.0, 111.0]. However, note that 
;        while it might be a good idea to identify the location of your line, relative 
;        to the whole spectrum, by setting the step size in frequency space to a 
;        larger number, once you identify the location of your highest S/N line you'll 
;        want to then use a more focused call to fitline_mcmc, with a smaller fwindow, 
;        frange, and smaller step sizes to adequately sample the statistical distribution 
;        of frequency space and get a more accurate estimate of the line.
;
; base_poly_ord - order of the polynomial to fit to the spectrum to remove a low-
;        frequency baseline prior to analysis - default=2
;
;  baseplotname - file name of the EPS figure showing the aforementioned "baseplot"
;
; local_baseline - calculate the baseline of the spectrum over only the part of the 
;                  spectrum within fwindow (the default is to fit the baseline 
;                  over the entire spectrum). 
;
; nobasesubtract - boolean to NOT subtract off a baseline from the spectrum. NOTE 
;        that it will still calculate a baseline, but it just won't actually 
;        apply the baseline subtraction to the spectrum.
;
; first_guess_gaussian - 4-element array containing the parameters of a Gaussian 
;        to subtract from the spectrum prior to removing the baseline. By default,
;        the code will attempt to fit a Gaussian, near the specified frequency 
;        using its MCMC framework, and then repeat the whole fitting analysis by 
;        first removing that initial best-fit Gaussian to measure a more accurate 
;        baseline, and then running its MCMC stuff on the baseline-subtracted 
;        spectrum. However, if you have prior knowledge of what your line SHOULD
;        look like, you can skip the initial run-through by specifying the param-
;        eters of the Gaussian describing the emission line.
;
; savsetname - name of an IDL saveset to output all the relevant results, as well 
;        as the spectrum
;
; save_mcmc_steps - if flagged, the output saveset "savsetname" will also include 
;        the MCMC steps taken for each parameter (note that this can add a great 
;        deal to the size of the output savesets!)
;
;  noplots - suppress all graphics (to window or otherwise)
;
;  gain - the telescope/instrument gain (in mJy/mK) used to convert the units of the 
;         spectrum. For the RSR on the LMT, it's 7.0, which is the default used here.
;
; MODIFICATION HISTORY:
; June 25, 2015 - version 1 finished
;  July 7, 2015 - Fixed a typo in the header, and added /noplots option
; July 21, 2015 - modified the integration over the Gaussian line to round to the
;                 nearest pixel, rather than always rounding down. Also added ffp_fwhm
;                 array to the saveset output, which gives the resulting Gaussian fit 
;                 of the width of the line in terms of velocity.
; July 17, 2016 - Added a "gain" parameter so that the user can potentially apply this 
;                 to spectra from other instruments/facilities
;
; EXAMPLE:
; 
;      ; Locate the peak S/N line in a wide frequency window centered on 98.605 GHz:
;
; IDL> fitline_mcmc, 'spectra/A963_noHI_04_nofilt.txt', 98.605, frange=[97.6, 99.6]
;
;      ; This search turns up a line centered on f0=98.601 GHz, so you can do a more
;      ; focused fitting of that line (the default frange will use f0 +/- 0.08 GHz):
;
; IDL> fitline_mcmc, 'spectra/A963_noHI_04_nofilt.txt', 98.601
;
;      ; The previous call to fitline_mcmc fit a Gaussian with parameters:
;      ; [1.005, 98.6006, 0.0448, -0.0622] and a S/N of 8.178. However, we can 
;      ; iterate on this once more by using those best-fit parameters as the 
;      ; starting Gaussian that gets removed from the spectrum to fit a better 
;      ; baseline on the next run-through.
;
; IDL> ffgau = [1.005, 98.6006, 0.0448, -0.0622]
; IDL> fitline_mcmc, 'spectra/A963_noHI_04_nofilt.txt', 98.601, first_guess_gaussian=ffgau
;
;      ; Now we have a slightly better fit to the line, with S/N=8.906. 
;  
