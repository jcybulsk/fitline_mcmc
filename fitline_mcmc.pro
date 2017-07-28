; Here's our Gaussian
FUNCTION thegaussian, x, a
	z = (x-a[1])/a[2]
	y = a[0]*exp(-z^2/2) + a[3]
	return, y
END

; This function handles the probability (likelihood calculation)
FUNCTION pi, x, y, a, thesig
	npts = n_elements(x)
	gau = thegaussian(x,a)
	tot=0.0
	for i=0d,npts-1 do begin
		tot = tot + alog10(exp(-((y[i]-gau[i])^2)/0.001*thesig^2))
	endfor
	return, tot
END

; With this function, we decide what the next step for each parameter (if taken) will be
FUNCTION stepcalc, a, b0sig, b1sig, b2sig, b3sig, a0lims, a1lims, a2lims, a3lims
	b = dblarr(4)
	b0test = randomn(seed)*b0sig+a[0]
	b1test = randomn(seed)*b1sig+a[1]
	b2test = randomn(seed)*b2sig+a[2]
	b3test = randomn(seed)*b3sig+a[3]
	
	good0 = -1.0
	good1 = -1.0
	good2 = -1.0
	good3 = -1.0
		
	while good0 lt 0 do begin
		if b0test gt a0lims[0] AND b0test lt a0lims[1] then begin
			good0 = 1.0
		endif else begin
			b0test = randomn(seed)*b0sig+a[0]
		endelse
	endwhile

	while good1 lt 0 do begin
		if b1test gt a1lims[0] AND b1test lt a1lims[1] then begin
			good1 = 1.0
		endif else begin
			b1test = randomn(seed)*b1sig+a[1]
		endelse
	endwhile
	
	while good2 lt 0 do begin
		if b2test gt a2lims[0] AND b2test lt a2lims[1] then begin
			good2 = 1.0
		endif else begin
			b2test = randomn(seed)*b2sig+a[2]
		endelse
	endwhile

	while good3 lt 0 do begin
		if b3test gt a3lims[0] AND b3test lt a3lims[1] then begin
			good3 = 1.0
		endif else begin
			b3test = randomn(seed)*b3sig+a[3]
		endelse
	endwhile

	b[0] = b0test
	b[1] = b1test
	b[2] = b2test
	b[3] = b3test
	return, b
END

PRO fitline_mcmc, filename, f0, g0=g0, frange=frange, normrange=normrange, sigrange=sigrange, offsetrange=offsetrange, fwindow=fwindow, nsteps=nsteps, nburnin=nburnin, savefigs=savefigs, specplotname=specplotname, paramplotname=paramplotname, showsteps=showsteps, stepplotname=stepplotname, cplot=cplot, contplotname=contplotname, movie=movie, stepsizes=stepsizes, base_poly_ord=base_poly_ord, baseplotname=baseplotname, local_baseline=local_baseline, nobasesubtract=nobasesubtract, first_guess_gaussian=first_guess_gaussian, savsetname=savsetname, save_mcmc_steps=save_mcmc_steps, noplots=noplots, gain=gain

if(not keyword_set(filename)) then message, 'ERROR: No spectrum ASCII file specified!'
if(not keyword_set(f0)) then message, 'ERROR: No frequency specified!'
if(not keyword_set(fwindow)) then fwindow=[f0-1.0, f0+1.0]
if fwindow[0] ge fwindow[1] then message, 'ERROR: fwindow array must be [f_lower, f_upper]!'

if(keyword_set(specplotname) OR keyword_set(paramplotname) OR keyword_set(stepplotname) OR keyword_set(contplotname) OR keyword_set(baseplotname)) then savefigs=1
if(keyword_set(contplotname)) then cplot=1
if(keyword_set(stepplotname)) then showsteps=1

loadct, 13, ncolors=255
if(keyword_set(savefigs)) then set_plot, 'ps'

; Read in the spectrum
FMT = 'F, F'
readcol,filename,F=FMT,sfreq,svalue,/silent

; Replace nans with 0.0, and convert from K to mK
inan = where(finite(svalue, /nan))
if n_elements(inan) gt 1 then svalue[inan] = 0.0
freq_full = sfreq
y_full = svalue*1000.0

; We've built in an approach here to calculate a first-time fitting to improve on the 
; baseline subtraction. So basically the fitting process runs twice: The first time, it 
; just tries to find a Gaussian that fits the line, and the second time it subtracts that 
; prior best-fitting Gaussian, THEN fits the baseline with the line subtracted, AND THEN 
; it finally runs the line fitting again, but the second time it should be more reliable 
; because it's done a more thorough job of fitting the baseline. Note that the user can 
; specify their own Gaussian parameters to use to skip that first time around.
first_time_done = -1.0
for z=0,1 do begin

if(keyword_set(first_guess_gaussian)) then z=1

; Trim the data to the frequency range of interest
x = sfreq[where(sfreq ge fwindow[0] AND sfreq le fwindow[1])]
ytrim = svalue[where(sfreq ge fwindow[0] AND sfreq le fwindow[1])]*1000.0    ; Put spectrum in mK
x_raw = x
y_raw = ytrim
ndata = double(n_elements(ytrim))
ndata_full = double(n_elements(y_full))

if(not keyword_set(base_poly_ord)) then base_poly_ord=2

if(z eq 1) then begin
	if(keyword_set(first_guess_gaussian)) then begin
		print, 'I hear you and I will use the user-specified first guess Gaussian with parameters: ', +string(first_guess_gaussian)
		if(keyword_set(local_baseline)) then begin
			first_bestfit_gaussian = thegaussian(x,first_guess_gaussian)
		endif else begin
			first_bestfit_gaussian = thegaussian(freq_full,first_guess_gaussian)
		endelse
	endif else begin
		print, 'Initial run finds a first guess Gaussian with parameters: ', +string([ffp_amp[1], ffp_nu0[1], ffp_sig[1], ffp_off[1]])
		if(keyword_set(local_baseline)) then begin
			first_bestfit_gaussian = thegaussian(x,[ffp_amp[1], ffp_nu0[1], ffp_sig[1], ffp_off[1]])
		endif else begin
			first_bestfit_gaussian = thegaussian(freq_full,[ffp_amp[1], ffp_nu0[1], ffp_sig[1], ffp_off[1]])
		endelse		
	endelse

	; The baseline is fit using Savitzky-Golay filter, which is a rolling polynomial fit (default order of 2), 
	; which is accomplished using the IDL routine poly_smooth. It's very important not to "over-fit", though, 
	; and so the baseline should be calculated over a pretty wide range in frequency (approximately 1GHz should 
	; be used for the width, which is why I have for a "local" calculation (+/- 1 GHz) ndata/2.0, and for a 
	; baseline fit over the entire frequency range (~40 GHz) ndata_full/40.0.
	if(keyword_set(local_baseline)) then begin
		ysub = ytrim-first_bestfit_gaussian
		x_baseline = x
		baseline = poly_smooth(ysub, ndata/2.0, degree=base_poly_ord)
		y = ytrim-baseline
	
		sig_1 = stdev(ytrim)
		sig_2 = stdev(y)
	endif else begin
		ysub = y_full-first_bestfit_gaussian
		x_baseline = freq_full
		baseline = poly_smooth(ysub, ndata_full/40.0, degree=base_poly_ord)
		y = y_full-baseline
	
		sig_1_full = stdev(y_full)
		sig_2_full = stdev(y)
		y = y[where(sfreq ge fwindow[0] AND sfreq le fwindow[1])]
		sig_1 = stdev(ytrim)
		sig_2 = stdev(y)		
	endelse

	!P.Multi = [0,1,2]
	if(keyword_set(savefigs)) then begin
		if(not keyword_set(baseplotname)) then baseplotname='baseline_RSR_plot.eps'
		device, filename=baseplotname, /encapsulated, /color
		plot, x, ytrim, xst=1, yst=1, psym=10, xr=fwindow, yr=[min(ytrim)-0.1, max(ytrim)+0.1], xtitle='Frequency (GHz)', ytitle=textoidl('T_A^* (mK)'), charsize=1.3, charthick=1.3, title='Spectrum + Baseline'
		plots, fwindow[0], 0
		plots, fwindow[1], 0, /continue, thick=2
		plots, f0, min(ytrim)-0.1
		plots, f0, max(ytrim)+0.1, /continue, color=254, linestyle=1
		oplot, x_baseline, baseline, color=254
		xyouts, fwindow[0]+0.1, max(ytrim)-0.1, strcompress('rms = '+string(format='(F8.3)',sig_1)+' mK')

		plot, x, y, xst=1, yst=1, psym=10, xr=fwindow, yr=[min(ytrim)-0.1, max(ytrim)+0.1], xtitle='Frequency (GHz)', ytitle=textoidl('T_A^* (mK)'), charsize=1.3, charthick=1.3, title='Baseline-Subtracted Spectrum'
		plots, fwindow[0], 0
		plots, fwindow[1], 0, /continue, thick=2
		plots, f0, min(ytrim)-0.1
		plots, f0, max(ytrim)+0.1, /continue, color=254, linestyle=1
		xyouts, fwindow[0]+0.1, max(ytrim)-0.1, strcompress('rms = '+string(format='(F8.3)',sig_2)+' mK')
		device, /close
	endif else begin
	if not(keyword_set(noplots)) then begin
		window, 6, retain=2
		plot, x, ytrim, xst=1, yst=1, psym=10, xr=fwindow, yr=[min(ytrim)-0.1, max(ytrim)+0.1], xtitle='Frequency (GHz)', ytitle=textoidl('T_A^* (mK)'), charsize=1.3, charthick=1.3, title='Spectrum + Baseline'
		plots, fwindow[0], 0
		plots, fwindow[1], 0, /continue, thick=2
		plots, f0, min(ytrim)-0.1
		plots, f0, max(ytrim)+0.1, /continue, color=254, linestyle=1
		oplot, x_baseline, baseline, color=254
		xyouts, fwindow[0]+0.1, max(ytrim)-0.1, strcompress('rms = '+string(format='(F8.3)',sig_1)+' mK')

		plot, x, y, xst=1, yst=1, psym=10, xr=fwindow, yr=[min(ytrim)-0.1, max(ytrim)+0.1], xtitle='Frequency (GHz)', ytitle=textoidl('T_A^* (mK)'), charsize=1.3, charthick=1.3, title='Baseline-Subtracted Spectrum'
		plots, fwindow[0], 0
		plots, fwindow[1], 0, /continue, thick=2
		plots, f0, min(ytrim)-0.1
		plots, f0, max(ytrim)+0.1, /continue, color=254, linestyle=1
		xyouts, fwindow[0]+0.1, max(ytrim)-0.1, strcompress('rms = '+string(format='(F8.3)',sig_2)+' mK')
	endif
	endelse
	!P.Multi = 0
endif else begin
	y = ytrim
endelse

if(keyword_set(nobasesubtract)) then y=ytrim

; Here's where most of the default parameter values are set, if they're not specified by the user:
if(not keyword_set(g0)) then g0=[0.3, f0, 0.015, 0.01]
if(not keyword_set(frange)) then frange=[f0-0.08, f0+0.08]
if(not keyword_set(normrange)) then normrange=[0.0, 100.0]
if(not keyword_set(sigrange)) then sigrange=[0.0, 0.15]
if(not keyword_set(offsetrange)) then offsetrange=[-100.0,100.0]
if(not keyword_set(stepsizes)) then stepsizes=[1.0, 1.0, 1.0, 1.0]
if(not keyword_set(gain)) then gain=7.0

; Here we catch a few common mistakes that might be made, like specifying starting 
; parameter values that are outside of the specified ranges that can be explored.
if g0[0] lt normrange[0] OR g0[0] gt normrange[1] then g0[0]=mean(normrange)
if g0[1] lt frange[0] OR g0[1] gt frange[1] then g0[1]=mean(frange)
if g0[2] lt sigrange[0] OR g0[2] gt sigrange[1] then g0[2]=mean(sigrange)
if g0[3] lt offsetrange[0] OR g0[3] gt offsetrange[1] then g0[3]=mean(offsetrange)

;;; Apply certain limits to the parameter values?
a0lims = normrange
a1lims = frange
a2lims = sigrange
a3lims = offsetrange

if f0 lt a1lims[0] OR f0 gt a1lims[1] then message, 'ERROR: f0 is outside of the specified frequency range in "frange"!'

if(not keyword_set(nsteps)) then nsteps=100000d
if(not keyword_set(nburnin)) then nburnin=40000d
if nburnin ge nsteps then message, 'ERROR: # of burn-in steps must be < # of MCMC steps!"

ndata = n_elements(x)
nmc = nsteps

a0mc = dblarr(nmc)  ; parameter storage
a1mc = dblarr(nmc)
a2mc = dblarr(nmc)
a3mc = dblarr(nmc)

; starting guesses for the chain
a0mc[0] = g0[0]
a1mc[0] = g0[1]
a2mc[0] = g0[2]
a3mc[0] = g0[3]

a = [a0mc[0], a1mc[0], a2mc[0], a3mc[0]]

;;; The following standard deviations are used to set the scale of step sizes in all the 
;;; parameters. The user can modify these using the stepsizes array, but the defaults have 
;;; been tuned carefully to give a pretty reliable recovery of lines in RSR spectra.
b0sig = 0.05*stepsizes[0]     ; std deviation of the transition probability
b1sig = 0.4*stepsizes[1]      ; std deviation of the transition probability
b2sig = 0.005*stepsizes[2]    ; std deviation of the transition probability
b3sig = 0.01*stepsizes[3]     ; std deviation of the transition probability
thesig = 0.2

b0mc = dblarr(nmc)
b1mc = dblarr(nmc)
b2mc = dblarr(nmc)
b3mc = dblarr(nmc)
alpha = dblarr(nmc)

plotct = 0.0
firsttime = -1

for i=0d,nmc-2 do begin
	;define the a vector
	a[0] = a0mc[i]
	a[1] = a1mc[i]
	a[2] = a2mc[i]
	a[3] = a3mc[i]

	if(keyword_set(movie) AND z eq 1) then begin
		if plotct eq 1000.0 then begin
			current_gaussian = thegaussian(x,a)
			if i lt nburnin then begin
				thetitle = strcompress('Step number '+string(i,format='(I7)'))
			endif else begin
				thetitle = strcompress('Step number '+string(i,format='(I7)') + ' burn-in complete')
			endelse
			if firsttime lt 0 then begin
				set_plot, 'x'
				window, 2, retain=2
				firsttime=1
			endif
			plot, x, y, xst=1, yst=1, psym=10, xtitle='Frequency (GHz)', ytitle=textoidl('T_{A}^* (mK)'), charsize=1.3, charthick=1.3, title=thetitle, yr=[min(y)-0.1, max(y)+0.1]
			oplot, x, current_gaussian, color=254, thick=2
			wait, 0.75
			plotct = 0.0
		endif
	
		plotct = plotct + 1.0
	endif

	;take a potential step
	b = stepcalc(a,b0sig,b1sig,b2sig,b3sig,a0lims,a1lims,a2lims,a3lims)
	b0mc[i] = b[0]
	b1mc[i] = b[1]
	b2mc[i] = b[2]
	b3mc[i] = b[3]

	;calculate acceptance probability
	pib = pi(x,y,b,thesig)

	;likelihood of candidate step
	pia = pi(x,y,a,thesig)

	;likelihood of current position
	alpha[i] = min([0,pib-pia])

	;decide whether to accept or not
	u = alog10(randomu(seed))
	if(u le alpha[i]) then begin
		;looks good, make the step
		a0mc[i+1]=b[0]
		a1mc[i+1]=b[1]
		a2mc[i+1]=b[2]
		a3mc[i+1]=b[3]
	endif else begin
		;no good, try again from same position
		a0mc[i+1]=a0mc[i]
		a1mc[i+1]=a1mc[i]
		a2mc[i+1]=a2mc[i]
		a3mc[i+1]=a3mc[i]
	endelse
endfor

; Now trim off the pre-burn-in stuff
a0mc_save = a0mc[(nburnin-1.0):nmc-1]
a1mc_save = a1mc[(nburnin-1.0):nmc-1]
a2mc_save = a2mc[(nburnin-1.0):nmc-1]
a3mc_save = a3mc[(nburnin-1.0):nmc-1]
alpha_save = alpha[(nburnin-1.0):nmc-1]

; Convert the amplitude and DC offsets from mK to mJy, and the std deviation of the line to a FWHM in km/s
a0mc_save_mJy = a0mc_save*gain
a2mc_save_fwhm = 2.3548*(a2mc_save/a1mc_save)*299792.5
a3mc_save_mJy = a3mc_save*gain

;;; histogram up the amplitude, frequency, sigma, and dc offset distributions, and then fit gaussians to them!
; First identify the overall maximum-likelihood parameters, and center the distributions around those:
best_fit = where(alpha_save eq min(alpha_save))
a0best = a0mc_save[best_fit[0]]
a1best = a1mc_save[best_fit[0]]
a2best = a2mc_save[best_fit[0]]
a3best = a3mc_save[best_fit[0]]
abest = [a0best, a1best, a2best, a3best]
best_gaussian = thegaussian(x,abest)

fwhm_best = a2mc_save_fwhm[best_fit[0]]
a0best_mJy = a0mc_save_mJy[best_fit[0]]

ymax = max([y, best_gaussian])
ymin = min([y, best_gaussian])

xhist_amp = (dindgen(201)/20.0)+a0best
xhist_nu0 = (dindgen(201)/500.0)+a1best
xhist_sig = (dindgen(201)/500.0)+a2best
xhist_off = (dindgen(201)/250.0)+a3best
xhist_ampmJy = (dindgen(201)/10.0)+a0best_mJy
xhist_fwhm = (dindgen(201)*4.0)+fwhm_best

xhistrange_amp = max(xhist_amp)-min(xhist_amp)
xhistrange_nu0 = max(xhist_nu0)-min(xhist_nu0)
xhistrange_sig = max(xhist_sig)-min(xhist_sig)
xhistrange_off = max(xhist_off)-min(xhist_off)
xhistrange_ampmJy = max(xhist_ampmJy)-min(xhist_ampmJy)
xhistrange_fwhm = max(xhist_fwhm)-min(xhist_fwhm)

xhist_amp = xhist_amp-(xhistrange_amp/2.0)
xhist_nu0 = xhist_nu0-(xhistrange_nu0/2.0)
xhist_sig = xhist_sig-(xhistrange_sig/2.0)
xhist_off = xhist_off-(xhistrange_off/2.0)
xhist_ampmJy = xhist_ampmJy-(xhistrange_ampmJy/2.0)
xhist_fwhm = xhist_fwhm-(xhistrange_fwhm/2.0)

nhist_amp = n_elements(xhist_amp)
nhist_nu0 = n_elements(xhist_nu0)
nhist_sig = n_elements(xhist_sig)
nhist_off = n_elements(xhist_off)
nhist_ampmJy = n_elements(xhist_ampmJy)
nhist_fwhm = n_elements(xhist_fwhm)

yhist_amp = dblarr(nhist_amp)
yhist_nu0 = dblarr(nhist_nu0)
yhist_sig = dblarr(nhist_sig)
yhist_off = dblarr(nhist_off)
yhist_ampmJy = dblarr(nhist_ampmJy)
yhist_fwhm = dblarr(nhist_fwhm)

xsize_amp = abs(xhist_amp[1]-xhist_amp[0])
xsize_nu0 = abs(xhist_nu0[1]-xhist_nu0[0])
xsize_sig = abs(xhist_sig[1]-xhist_sig[0])
xsize_off = abs(xhist_off[1]-xhist_off[0])
xsize_ampmJy = abs(xhist_ampmJy[1]-xhist_ampmJy[0])
xsize_fwhm = abs(xhist_fwhm[1]-xhist_fwhm[0])

for i=0,nhist_amp-1 do begin
	test = where(a0mc_save ge (xhist_amp[i]-(xsize_amp/2.0)) AND a0mc_save le (xhist_amp[i]+(xsize_amp/2.0)),num_amp)
	test = where(a1mc_save ge (xhist_nu0[i]-(xsize_nu0/2.0)) AND a1mc_save le (xhist_nu0[i]+(xsize_nu0/2.0)),num_nu0)
	test = where(a2mc_save ge (xhist_sig[i]-(xsize_sig/2.0)) AND a2mc_save le (xhist_sig[i]+(xsize_sig/2.0)),num_sig)
	test = where(a3mc_save ge (xhist_off[i]-(xsize_off/2.0)) AND a3mc_save le (xhist_off[i]+(xsize_off/2.0)),num_off)
	test = where(a0mc_save_mJy ge (xhist_ampmJy[i]-(xsize_ampmJy/2.0)) AND a0mc_save_mJy le (xhist_ampmJy[i]+(xsize_ampmJy/2.0)),num_ampmJy)
	test = where(a2mc_save_fwhm ge (xhist_fwhm[i]-(xsize_fwhm/2.0)) AND a2mc_save_fwhm le (xhist_fwhm[i]+(xsize_fwhm/2.0)),num_fwhm)

	yhist_amp[i] = double(num_amp)
	yhist_nu0[i] = double(num_nu0)
	yhist_sig[i] = double(num_sig)
	yhist_off[i] = double(num_off)
	yhist_ampmJy[i] = double(num_ampmJy)
	yhist_fwhm[i] = double(num_fwhm)
endfor

;;; Fit gaussians to the histograms of parameters to get the central value, and the standard deviations:
xpeak_amp = where(yhist_amp eq max(yhist_amp))
xpeak_nu0 = where(yhist_nu0 eq max(yhist_nu0))
xpeak_sig = where(yhist_sig eq max(yhist_sig))
xpeak_off = where(yhist_off eq max(yhist_off))
xpeak_ampmJy = where(yhist_ampmJy eq max(yhist_ampmJy))
xpeak_fwhm = where(yhist_fwhm eq max(yhist_fwhm))

sp_amp = [max(yhist_amp), xhist_amp[xpeak_amp[0]], 0.2, 0.0]
parinfo_amp = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.d, 0]},4)
parinfo_amp[*].value = sp_amp
parinfo_amp[3].fixed = 1
parinfo_amp[2].limited = [1,0]
parinfo_amp[2].limits[0] = 0.0
err_amp = dblarr(nhist_amp)+1.0
ffp_amp = mpfitfun('thegaussian', xhist_amp, yhist_amp, err_amp, parinfo=parinfo_amp, bestnorm=chi2, dof=dof, perror=fperr, /quiet)
fitvals_amp = thegaussian(xhist_amp, ffp_amp)

sp_nu0 = [max(yhist_nu0), xhist_nu0[xpeak_nu0[0]], 0.01, 0.0]
parinfo_nu0 = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.d, 0]},4)
parinfo_nu0[*].value = sp_nu0
parinfo_nu0[3].fixed = 1
parinfo_nu0[2].limited = [1,0]
parinfo_nu0[2].limits[0] = 0.0
err_nu0 = dblarr(nhist_nu0)+1.0
ffp_nu0 = mpfitfun('thegaussian', xhist_nu0, yhist_nu0, err_nu0, parinfo=parinfo_nu0, bestnorm=chi2, dof=dof, perror=fperr, /quiet)
fitvals_nu0 = thegaussian(xhist_nu0, ffp_nu0)

sp_sig = [max(yhist_sig), xhist_sig[xpeak_sig[0]], 0.01, 0.0]
parinfo_sig = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.d, 0]},4)
parinfo_sig[*].value = sp_sig
parinfo_sig[3].fixed = 1
parinfo_sig[2].limited = [1,0]
parinfo_sig[2].limits[0] = 0.0
err_sig = dblarr(nhist_sig)+1.0
ffp_sig = mpfitfun('thegaussian', xhist_sig, yhist_sig, err_sig, parinfo=parinfo_sig, bestnorm=chi2, dof=dof, perror=fperr, /quiet)
fitvals_sig = thegaussian(xhist_sig, ffp_sig)

sp_off = [max(yhist_off), xhist_off[xpeak_off[0]], 0.1, 0.0]
parinfo_off = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.d, 0]},4)
parinfo_off[*].value = sp_off
parinfo_off[3].fixed = 1
parinfo_off[2].limited = [1,0]
parinfo_off[2].limits[0] = 0.0
err_off = dblarr(nhist_off)+1.0
ffp_off = mpfitfun('thegaussian', xhist_off, yhist_off, err_off, parinfo=parinfo_off, bestnorm=chi2, dof=dof, perror=fperr, /quiet)
fitvals_off = thegaussian(xhist_off, ffp_off)

sp_ampmJy = [max(yhist_ampmJy), xhist_ampmJy[xpeak_ampmJy[0]], 1.0, 0.0]
parinfo_ampmJy = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.d, 0]},4)
parinfo_ampmJy[*].value = sp_ampmJy
parinfo_ampmJy[3].fixed = 1
parinfo_ampmJy[2].limited = [1,0]
parinfo_ampmJy[2].limits[0] = 0.0
err_ampmJy = dblarr(nhist_ampmJy)+1.0
ffp_ampmJy = mpfitfun('thegaussian', xhist_ampmJy, yhist_ampmJy, err_ampmJy, parinfo=parinfo_ampmJy, bestnorm=chi2, dof=dof, perror=fperr, /quiet)
fitvals_ampmJy = thegaussian(xhist_ampmJy, ffp_ampmJy)

sp_fwhm = [max(yhist_fwhm), xhist_fwhm[xpeak_fwhm[0]], 100.0, 0.0]
parinfo_fwhm = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.d, 0]},4)
parinfo_fwhm[*].value = sp_fwhm
parinfo_fwhm[3].fixed = 1
parinfo_fwhm[2].limited = [1,0]
parinfo_fwhm[2].limits[0] = 0.0
err_fwhm = dblarr(nhist_fwhm)+1.0
ffp_fwhm = mpfitfun('thegaussian', xhist_fwhm, yhist_fwhm, err_fwhm, parinfo=parinfo_fwhm, bestnorm=chi2, dof=dof, perror=fperr, /quiet)
fitvals_fwhm = thegaussian(xhist_fwhm, ffp_fwhm)

;;; Print the results, and make plots!
if(z eq 1) then begin
	!P.Multi = [0,3,2]
	if(keyword_set(savefigs)) then begin
		set_plot, 'ps'
		if(not keyword_set(paramplotname)) then paramplotname='MCMC_parameter_histograms.eps'
		device, filename=paramplotname, /encapsulated, /color
		plot, xhist_amp, yhist_amp, xst=1, yst=1, psym=10, yr=[0,max(yhist_amp)+5.0], xtitle='Amplitude (mK)', ytitle='Num', charsize=1.3, charthick=1.3
		oplot, xhist_amp, fitvals_amp, color=254, thick=2

		plot, xhist_nu0, yhist_nu0, xst=1, yst=1, psym=10, yr=[0,max(yhist_nu0)+5.0], xtitle=textoidl('\nu_{CO} (GHz)'), ytitle='Num', charsize=1.3, charthick=1.3
		oplot, xhist_nu0, fitvals_nu0, color=254, thick=2

		plot, xhist_sig, yhist_sig, xst=1, yst=1, psym=10, yr=[0,max(yhist_sig)+5.0], xtitle='Std. Deviation (GHz)', ytitle='Num', charsize=1.3, charthick=1.3
		oplot, xhist_sig, fitvals_sig, color=254, thick=2

		plot, xhist_off, yhist_off, xst=1, yst=1, psym=10, yr=[0,max(yhist_off)+5.0], xtitle='Offset (mK)', ytitle='Num', charsize=1.3, charthick=1.3
		oplot, xhist_off, fitvals_off, color=254, thick=2

		plot, xhist_ampmJy, yhist_ampmJy, xst=1, yst=1, psym=10, yr=[0,max(yhist_ampmJy)+5.0], xtitle='Amplitude (mJy)', ytitle='Num', charsize=1.3, charthick=1.3
		oplot, xhist_ampmJy, fitvals_ampmJy, color=254, thick=2

		plot, xhist_fwhm, yhist_fwhm, xst=1, yst=1, psym=10, yr=[0,max(yhist_fwhm)+5.0], xtitle='FWHM (km/s)', ytitle='Num', charsize=1.3, charthick=1.3
		oplot, xhist_fwhm, fitvals_fwhm, color=254, thick=2
		device, /close
	endif else begin
	if not(keyword_set(noplots)) then begin
		window, 1, retain=2, xsize=900, ysize=600
		plot, xhist_amp, yhist_amp, xst=1, yst=1, psym=10, yr=[0,max(yhist_amp)+5.0], xtitle='Amplitude (mK)', ytitle='Num', charsize=1.3, charthick=1.3
		oplot, xhist_amp, fitvals_amp, color=254, thick=2

		plot, xhist_nu0, yhist_nu0, xst=1, yst=1, psym=10, yr=[0,max(yhist_nu0)+5.0], xtitle=textoidl('\nu_{CO} (GHz)'), ytitle='Num', charsize=1.3, charthick=1.3
		oplot, xhist_nu0, fitvals_nu0, color=254, thick=2

		plot, xhist_sig, yhist_sig, xst=1, yst=1, psym=10, yr=[0,max(yhist_sig)+5.0], xtitle='Std. Deviation (GHz)', ytitle='Num', charsize=1.3, charthick=1.3
		oplot, xhist_sig, fitvals_sig, color=254, thick=2

		plot, xhist_off, yhist_off, xst=1, yst=1, psym=10, yr=[0,max(yhist_off)+5.0], xtitle='Offset (mK)', ytitle='Num', charsize=1.3, charthick=1.3
		oplot, xhist_off, fitvals_off, color=254, thick=2

		plot, xhist_ampmJy, yhist_ampmJy, xst=1, yst=1, psym=10, yr=[0,max(yhist_ampmJy)+5.0], xtitle='Amplitude (mJy)', ytitle='Num', charsize=1.3, charthick=1.3
		oplot, xhist_ampmJy, fitvals_ampmJy, color=254, thick=2

		plot, xhist_fwhm, yhist_fwhm, xst=1, yst=1, psym=10, yr=[0,max(yhist_fwhm)+5.0], xtitle='FWHM (km/s)', ytitle='Num', charsize=1.3, charthick=1.3
		oplot, xhist_fwhm, fitvals_fwhm, color=254, thick=2
	endif
	endelse
	!P.Multi = 0
endif


str_amp = strcompress('Overall best-fit Gaussian has amp: ' + string(ffp_amp[1], format='(F8.3)') + ' +/- ' + string(ffp_amp[2], format='(F8.3)') + ' mK')
str_ampmJy = strcompress('Overall best-fit Gaussian has amp: ' + string(ffp_ampmJy[1], format='(F8.3)') + ' +/- ' + string(ffp_ampmJy[2], format='(F8.3)') + ' mJy')
str_nu0 = strcompress('Overall best-fit Gaussian has f0: ' + string(ffp_nu0[1], format='(F9.4)') + ' +/- ' + string(ffp_nu0[2], format='(F9.4)') + ' GHz')
str_sig = strcompress('Overall best-fit Gaussian has sigma: ' + string(ffp_sig[1], format='(F8.4)') + ' +/- ' + string(ffp_sig[2], format='(F8.4)') + ' GHz')
str_off = strcompress('Overall best-fit Gaussian has DC offset: ' + string(ffp_off[1], format='(F8.4)') + ' +/- ' + string(ffp_off[2], format='(F8.4)') + ' mK')
str_fwhm = strcompress('Overall best-fit Gaussian has FWHM: ' + string(ffp_fwhm[1], format='(F9.3)') + ' +/- ' + string(ffp_fwhm[2], format='(F9.3)') + ' km/s')

;;;; Now get the integrated line flux, and associated error
npeak = value_locate(x, ffp_nu0[1]) ; find the peak pixel
nline = round(abs(2.3548*ffp_sig[1]/31.25e-3))
if nline eq 0 then nline=1
dv=299792.5*31.25e-3/ffp_nu0[1]
linetotal=(dv*(total(y[npeak-nline:npeak+nline])-ffp_off[1]*2*nline*1.e-3))/1000.0
y_noline1 = y[0:(npeak-nline)]
y_noline2 = y[(npeak+nline):(ndata-1)]
y_noline = [y_noline1,y_noline2]
; To calculate the rms leave off the line itself?
rms = stdev(y_noline)
rms_string = strcompress('RMS (excluding the line): '+string(rms, format='(F9.3)')+' mK over the frequency range ' + string(fwindow[0],format='(F9.3)') + ' to ' + string(fwindow[1], format='(F9.3)') + ' GHz')
sigmaline=dv*rms*sqrt(float(2.0*nline))*1.e-3
SNR = linetotal/sigmaline
intstring = strcompress('line integral (K km/s) = ' + string(linetotal, format='(F12.3)') + ' +/- ' + string(sigmaline, format='(F12.3)'))
intstring2 = strcompress('line integral (Jy km/s) = ' + string(linetotal*gain, format='(F14.3)') + ' +/- ' + string(sigmaline*gain, format='(F14.3)'))
intstring3 = strcompress('S/N = '+string(SNR,format='(F12.3)'))

if(z eq 1) then begin
	print, str_amp
	print, str_ampmJy
	print, str_nu0
	print, str_sig
	print, str_fwhm
	print, str_off
	print, rms_string
	print,intstring
	print,intstring2
	print,intstring3
endif

overall_best_gaussian = thegaussian(x, [ffp_amp[1],ffp_nu0[1],ffp_sig[1],ffp_off[1]])
if(z eq 1) then begin
	if(keyword_set(savefigs)) then begin
		if(not keyword_set(specplotname)) then specplotname='RSR_spectrum.eps'
		device, filename=specplotname, /encapsulated, /color
		plot, x, y, xst=1, yst=1, xr=fwindow, yr=[ymin-0.1, ymax+0.1], psym=10, xtitle='Frequency (GHz)', ytitle=textoidl('T_{A}^* (mK)'), charsize=1.2, charthick=1.2
		oplot, x, overall_best_gaussian, color=254, thick=2, linestyle=2
		plots, f0, ymin-0.1
		plots, f0, ymax+0.1, /continue, color=254
		plots, frange[0], ymin-0.1
		plots, frange[0], ymax+0.1, /continue, linestyle=1, color=254
		plots, frange[1], ymin-0.1
		plots, frange[1], ymax+0.1, /continue, linestyle=1, color=254
		device, /close
	endif else begin
	if not(keyword_set(noplots)) then begin
		window, 0, retain=2
		plot, x, y, xst=1, yst=1, xr=fwindow, yr=[ymin-0.1, ymax+0.1], psym=10, xtitle='Frequency (GHz)', ytitle=textoidl('T_{A}^* (mK)'), charsize=1.2, charthick=1.2
		oplot, x, overall_best_gaussian, color=254, thick=2, linestyle=2
		plots, f0, ymin-0.1
		plots, f0, ymax+0.1, /continue, color=254
		plots, frange[0], ymin-0.1
		plots, frange[0], ymax+0.1, /continue, linestyle=1, color=254
		plots, frange[1], ymin-0.1
		plots, frange[1], ymax+0.1, /continue, linestyle=1, color=254
	endif
	endelse
endif

if(z eq 1) then begin
	!P.Multi = [0,3,2]
	if(keyword_set(showsteps)) then begin
		if(keyword_set(savefigs)) then begin
			if(not keyword_set(stepplotname)) then stepplotname = 'mcmc_steps_plot.eps'
			device, filename=stepplotname, /encapsulated, /color
			plot, a1mc, a0mc, xst=1, yst=1, xr=[min(a1mc)-0.01, max(a1mc)+0.01], yr=[min(a0mc)-0.1, max(a0mc)+0.1], xtitle=textoidl('\nu_{CO} (GHz)'), ytitle='Amplitude (mK)', charsize=1.3, charthick=1.3, /nodata, title='Red - Post Burn-In'
			oplot, a1mc, a0mc
			oplot, a1mc_save, a0mc_save, color=254
			plotsym, 3, 1.0, /fill
			plots, ffp_nu0[1], ffp_amp[1], psym=8, color=70

			plot, a2mc, a0mc, xst=1, yst=1, xr=[min(a2mc)-0.01, max(a2mc)+0.01], yr=[min(a0mc)-0.1, max(a0mc)+0.1], xtitle=textoidl('Std. Deviation (GHz)'), ytitle='Amplitude (mK)', charsize=1.3, charthick=1.3, /nodata, title='Red - Post Burn-In'
			oplot, a2mc, a0mc
			oplot, a2mc_save, a0mc_save, color=254
			plotsym, 3, 1.0, /fill
			plots, ffp_sig[1], ffp_amp[1], psym=8, color=70

			plot, a3mc, a0mc, xst=1, yst=1, xr=[min(a3mc)-0.05, max(a3mc)+0.05], yr=[min(a0mc)-0.1, max(a0mc)+0.1], xtitle=textoidl('DC Offset (mK)'), ytitle='Amplitude (mK)', charsize=1.3, charthick=1.3, /nodata, title='Red - Post Burn-In'
			oplot, a3mc, a0mc
			oplot, a3mc_save, a0mc_save, color=254
			plotsym, 3, 1.0, /fill
			plots, ffp_off[1], ffp_amp[1], psym=8, color=70

			plot, a2mc, a1mc, xst=1, yst=1, xr=[min(a2mc)-0.01, max(a2mc)+0.01], yr=[min(a1mc)-0.01, max(a1mc)+0.01], xtitle=textoidl('Std. Deviation (GHz)'), ytitle=textoidl('\nu_{CO} (GHz)'), charsize=1.3, charthick=1.3, /nodata, title='Red - Post Burn-In'
			oplot, a2mc, a1mc
			oplot, a2mc_save, a1mc_save, color=254
			plotsym, 3, 1.0, /fill
			plots, ffp_sig[1], ffp_nu0[1], psym=8, color=70

			plot, a3mc, a1mc, xst=1, yst=1, xr=[min(a3mc)-0.05, max(a3mc)+0.05], yr=[min(a1mc)-0.01, max(a1mc)+0.01], xtitle=textoidl('DC Offset (mK)'), ytitle=textoidl('\nu_{CO} (GHz)'), charsize=1.3, charthick=1.3, /nodata, title='Red - Post Burn-In'
			oplot, a3mc, a1mc
			oplot, a3mc_save, a1mc_save, color=254
			plotsym, 3, 1.0, /fill
			plots, ffp_off[1], ffp_nu0[1], psym=8, color=70

			plot, a3mc, a2mc, xst=1, yst=1, xr=[min(a3mc)-0.05, max(a3mc)+0.05], yr=[min(a2mc)-0.01, max(a2mc)+0.01], xtitle=textoidl('DC Offset (mK)'), ytitle=textoidl('Std. Deviation (GHz)'), charsize=1.3, charthick=1.3, /nodata, title='Red - Post Burn-In'
			oplot, a3mc, a2mc
			oplot, a3mc_save, a2mc_save, color=254
			plotsym, 3, 1.0, /fill
			plots, ffp_off[1], ffp_sig[1], psym=8, color=70
			device, /close
		endif else begin
		if not(keyword_set(noplots)) then begin
			window, 3, retain=2, xsize=900, ysize=600
			plot, a1mc, a0mc, xst=1, yst=1, xr=[min(a1mc)-0.01, max(a1mc)+0.01], yr=[min(a0mc)-0.1, max(a0mc)+0.1], xtitle=textoidl('\nu_{CO} (GHz)'), ytitle='Amplitude (mK)', charsize=1.3, charthick=1.3, /nodata, title='Red - Post Burn-In'
			oplot, a1mc, a0mc
			oplot, a1mc_save, a0mc_save, color=254
			plotsym, 3, 1.0, /fill
			plots, ffp_nu0[1], ffp_amp[1], psym=8, color=70

			plot, a2mc, a0mc, xst=1, yst=1, xr=[min(a2mc)-0.01, max(a2mc)+0.01], yr=[min(a0mc)-0.1, max(a0mc)+0.1], xtitle=textoidl('Std. Deviation (GHz)'), ytitle='Amplitude (mK)', charsize=1.3, charthick=1.3, /nodata, title='Red - Post Burn-In'
			oplot, a2mc, a0mc
			oplot, a2mc_save, a0mc_save, color=254
			plotsym, 3, 1.0, /fill
			plots, ffp_sig[1], ffp_amp[1], psym=8, color=70

			plot, a3mc, a0mc, xst=1, yst=1, xr=[min(a3mc)-0.05, max(a3mc)+0.05], yr=[min(a0mc)-0.1, max(a0mc)+0.1], xtitle=textoidl('DC Offset (mK)'), ytitle='Amplitude (mK)', charsize=1.3, charthick=1.3, /nodata, title='Red - Post Burn-In'
			oplot, a3mc, a0mc
			oplot, a3mc_save, a0mc_save, color=254
			plotsym, 3, 1.0, /fill
			plots, ffp_off[1], ffp_amp[1], psym=8, color=70

			plot, a2mc, a1mc, xst=1, yst=1, xr=[min(a2mc)-0.01, max(a2mc)+0.01], yr=[min(a1mc)-0.01, max(a1mc)+0.01], xtitle=textoidl('Std. Deviation (GHz)'), ytitle=textoidl('\nu_{CO} (GHz)'), charsize=1.3, charthick=1.3, /nodata, title='Red - Post Burn-In'
			oplot, a2mc, a1mc
			oplot, a2mc_save, a1mc_save, color=254
			plotsym, 3, 1.0, /fill
			plots, ffp_sig[1], ffp_nu0[1], psym=8, color=70

			plot, a3mc, a1mc, xst=1, yst=1, xr=[min(a3mc)-0.05, max(a3mc)+0.05], yr=[min(a1mc)-0.01, max(a1mc)+0.01], xtitle=textoidl('DC Offset (mK)'), ytitle=textoidl('\nu_{CO} (GHz)'), charsize=1.3, charthick=1.3, /nodata, title='Red - Post Burn-In'
			oplot, a3mc, a1mc
			oplot, a3mc_save, a1mc_save, color=254
			plotsym, 3, 1.0, /fill
			plots, ffp_off[1], ffp_nu0[1], psym=8, color=70

			plot, a3mc, a2mc, xst=1, yst=1, xr=[min(a3mc)-0.05, max(a3mc)+0.05], yr=[min(a2mc)-0.01, max(a2mc)+0.01], xtitle=textoidl('DC Offset (mK)'), ytitle=textoidl('Std. Deviation (GHz)'), charsize=1.3, charthick=1.3, /nodata, title='Red - Post Burn-In'
			oplot, a3mc, a2mc
			oplot, a3mc_save, a2mc_save, color=254
			plotsym, 3, 1.0, /fill
			plots, ffp_off[1], ffp_sig[1], psym=8, color=70
		endif
		endelse
	endif
	!P.Multi = 0
endif



if(z eq 1) then begin
	!P.Multi = [0,3,2]
	if(keyword_set(cplot)) then begin
		if(keyword_set(savefigs)) then begin
			if(not keyword_set(contplotname)) then contplotname = 'mcmc_contour_plot.eps'
			device, filename=contplotname, /encapsulated, /color
			amp_hist_keep = where(xhist_amp ge (ffp_amp[1]-3.0*ffp_amp[2]) AND xhist_amp le (ffp_amp[1]+3.0*ffp_amp[2]),nkeep_amp)
			nu0_hist_keep = where(xhist_nu0 ge (ffp_nu0[1]-3.0*ffp_nu0[2]) AND xhist_nu0 le (ffp_nu0[1]+3.0*ffp_nu0[2]),nkeep_nu0)
			sig_hist_keep = where(xhist_sig ge (ffp_sig[1]-3.0*ffp_sig[2]) AND xhist_sig le (ffp_sig[1]+3.0*ffp_sig[2]),nkeep_sig)
			off_hist_keep = where(xhist_off ge (ffp_off[1]-3.0*ffp_off[2]) AND xhist_off le (ffp_off[1]+3.0*ffp_off[2]),nkeep_off)

			xhist_amp_trim = xhist_amp[amp_hist_keep]
			yhist_amp_trim = yhist_amp[amp_hist_keep]
			xhist_nu0_trim = xhist_nu0[nu0_hist_keep]
			yhist_nu0_trim = yhist_nu0[nu0_hist_keep]
			xhist_sig_trim = xhist_sig[sig_hist_keep]
			yhist_sig_trim = yhist_sig[sig_hist_keep]
			xhist_off_trim = xhist_off[off_hist_keep]
			yhist_off_trim = yhist_off[off_hist_keep]

;;;; First, amplitude vs nu0
			cont1z = dblarr(nkeep_amp, nkeep_nu0)
			peakrow = -1.0
			peakcol = -1.0
			maxval = 0.0
			for i=0d,nkeep_amp-1 do begin
				for j=0d,nkeep_nu0-1 do begin
					test = where(a0mc_save ge (xhist_amp_trim[i]-xsize_amp/2.0) AND a0mc_save le (xhist_amp_trim[i]+xsize_amp/2.0) AND a1mc_save ge (xhist_nu0_trim[j]-xsize_nu0/2.0) AND a1mc_save le (xhist_nu0_trim[j]+xsize_nu0/2.0), ntest)
					cont1z[i,j] = double(ntest)
					if ntest gt maxval then begin
						maxval = double(ntest)
						peakrow = j
						peakcol = i
					endif
				endfor
			endfor
			interp_amp_1sig_min = interpol(transpose(cont1z[*,peakrow]), xhist_amp_trim, ffp_amp[1]-ffp_amp[2])
			interp_amp_1sig_plus = interpol(transpose(cont1z[*,peakrow]), xhist_amp_trim, ffp_amp[1]+ffp_amp[2])
			interp_nu0_1sig_min = interpol(transpose(cont1z[peakcol,*]), xhist_nu0_trim, ffp_nu0[1]-ffp_nu0[2])
			interp_nu0_1sig_plus = interpol(transpose(cont1z[peakcol,*]), xhist_nu0_trim, ffp_nu0[1]+ffp_nu0[2])
			thresh_1sig = mean([interp_amp_1sig_min, interp_amp_1sig_plus, interp_nu0_1sig_min, interp_nu0_1sig_plus])

			interp_amp_2sig_min = interpol(transpose(cont1z[*,peakrow]), xhist_amp_trim, ffp_amp[1]-2.0*ffp_amp[2])
			interp_amp_2sig_plus = interpol(transpose(cont1z[*,peakrow]), xhist_amp_trim, ffp_amp[1]+2.0*ffp_amp[2])
			interp_nu0_2sig_min = interpol(transpose(cont1z[peakcol,*]), xhist_nu0_trim, ffp_nu0[1]-2.0*ffp_nu0[2])
			interp_nu0_2sig_plus = interpol(transpose(cont1z[peakcol,*]), xhist_nu0_trim, ffp_nu0[1]+2.0*ffp_nu0[2])
			thresh_2sig = mean([interp_amp_2sig_min, interp_amp_2sig_plus, interp_nu0_2sig_min, interp_nu0_2sig_plus])

			contour, transpose(cont1z), xhist_nu0_trim, xhist_amp_trim, xst=1, yst=1, xr=[min(xhist_nu0_trim), max(xhist_nu0_trim)], yr=[min(xhist_amp_trim), max(xhist_amp_trim)], xtitle=textoidl('\nu_{CO} (GHz)'), ytitle='Amplitude (mK)', charsize=1.3, charthick=1.3, nlevels=2, levels=[thresh_2sig, thresh_1sig]
			plotsym, 3, 1.0, /fill
			plots, ffp_nu0[1], ffp_amp[1], psym=8, color=70

;;;; next comes amplitude vs sigma
			cont2z = dblarr(nkeep_amp, nkeep_sig)
			peakrow = -1.0
			peakcol = -1.0
			maxval = 0.0
			for i=0d,nkeep_amp-1 do begin
				for j=0d,nkeep_sig-1 do begin
					test = where(a0mc_save ge (xhist_amp_trim[i]-xsize_amp/2.0) AND a0mc_save le (xhist_amp_trim[i]+xsize_amp/2.0) AND a2mc_save ge (xhist_sig_trim[j]-xsize_sig/2.0) AND a2mc_save le (xhist_sig_trim[j]+xsize_sig/2.0), ntest)
					cont2z[i,j] = double(ntest)
					if ntest gt maxval then begin
						maxval = double(ntest)
						peakrow = j
						peakcol = i
					endif
				endfor
			endfor
			interp_amp_1sig_min = interpol(transpose(cont2z[*,peakrow]), xhist_amp_trim, ffp_amp[1]-ffp_amp[2])
			interp_amp_1sig_plus = interpol(transpose(cont2z[*,peakrow]), xhist_amp_trim, ffp_amp[1]+ffp_amp[2])
			interp_sig_1sig_min = interpol(transpose(cont2z[peakcol,*]), xhist_sig_trim, ffp_sig[1]-ffp_sig[2])
			interp_sig_1sig_plus = interpol(transpose(cont2z[peakcol,*]), xhist_sig_trim, ffp_sig[1]+ffp_sig[2])
			thresh_1sig = mean([interp_amp_1sig_min, interp_amp_1sig_plus, interp_sig_1sig_min, interp_sig_1sig_plus])

			interp_amp_2sig_min = interpol(transpose(cont2z[*,peakrow]), xhist_amp_trim, ffp_amp[1]-2.0*ffp_amp[2])
			interp_amp_2sig_plus = interpol(transpose(cont2z[*,peakrow]), xhist_amp_trim, ffp_amp[1]+2.0*ffp_amp[2])
			interp_sig_2sig_min = interpol(transpose(cont2z[peakcol,*]), xhist_sig_trim, ffp_sig[1]-2.0*ffp_sig[2])
			interp_sig_2sig_plus = interpol(transpose(cont2z[peakcol,*]), xhist_sig_trim, ffp_sig[1]+2.0*ffp_sig[2])
			thresh_2sig = mean([interp_amp_2sig_min, interp_amp_2sig_plus, interp_sig_2sig_min, interp_sig_2sig_plus])

			contour, transpose(cont2z), xhist_sig_trim, xhist_amp_trim, xst=1, yst=1, xr=[min(xhist_sig_trim), max(xhist_sig_trim)], yr=[min(xhist_amp_trim), max(xhist_amp_trim)], xtitle=textoidl('Std. Deviation (GHz)'), ytitle='Amplitude (mK)', charsize=1.3, charthick=1.3, nlevels=2, levels=[thresh_2sig, thresh_1sig]
			plotsym, 3, 1.0, /fill
			plots, ffp_sig[1], ffp_amp[1], psym=8, color=70

;;;; next comes amplitude vs DC Offset
			cont3z = dblarr(nkeep_amp, nkeep_off)
			peakrow = -1.0
			peakcol = -1.0
			maxval = 0.0
			for i=0d,nkeep_amp-1 do begin
				for j=0d,nkeep_off-1 do begin
					test = where(a0mc_save ge (xhist_amp_trim[i]-xsize_amp/2.0) AND a0mc_save le (xhist_amp_trim[i]+xsize_amp/2.0) AND a3mc_save ge (xhist_off_trim[j]-xsize_off/2.0) AND a3mc_save le (xhist_off_trim[j]+xsize_off/2.0), ntest)
					cont3z[i,j] = double(ntest)
					if ntest gt maxval then begin
						maxval = double(ntest)
						peakrow = j
						peakcol = i
					endif
				endfor
			endfor
			interp_amp_1sig_min = interpol(transpose(cont3z[*,peakrow]), xhist_amp_trim, ffp_amp[1]-ffp_amp[2])
			interp_amp_1sig_plus = interpol(transpose(cont3z[*,peakrow]), xhist_amp_trim, ffp_amp[1]+ffp_amp[2])
			interp_off_1sig_min = interpol(transpose(cont3z[peakcol,*]), xhist_off_trim, ffp_off[1]-ffp_off[2])
			interp_off_1sig_plus = interpol(transpose(cont3z[peakcol,*]), xhist_off_trim, ffp_off[1]+ffp_off[2])
			thresh_1sig = mean([interp_amp_1sig_min, interp_amp_1sig_plus, interp_off_1sig_min, interp_off_1sig_plus])

			interp_amp_2sig_min = interpol(transpose(cont3z[*,peakrow]), xhist_amp_trim, ffp_amp[1]-2.0*ffp_amp[2])
			interp_amp_2sig_plus = interpol(transpose(cont3z[*,peakrow]), xhist_amp_trim, ffp_amp[1]+2.0*ffp_amp[2])
			interp_off_2sig_min = interpol(transpose(cont3z[peakcol,*]), xhist_off_trim, ffp_off[1]-2.0*ffp_off[2])
			interp_off_2sig_plus = interpol(transpose(cont3z[peakcol,*]), xhist_off_trim, ffp_off[1]+2.0*ffp_off[2])
			thresh_2sig = mean([interp_amp_2sig_min, interp_amp_2sig_plus, interp_off_2sig_min, interp_off_2sig_plus])

			contour, transpose(cont3z), xhist_off_trim, xhist_amp_trim, xst=1, yst=1, xr=[min(xhist_off_trim), max(xhist_off_trim)], yr=[min(xhist_amp_trim), max(xhist_amp_trim)], xtitle=textoidl('DC Offset (mK)'), ytitle='Amplitude (mK)', charsize=1.3, charthick=1.3, nlevels=2, levels=[thresh_2sig, thresh_1sig]
			plotsym, 3, 1.0, /fill
			plots, ffp_off[1], ffp_amp[1], psym=8, color=70

;;;; next comes nu0 vs sigma
			cont4z = dblarr(nkeep_nu0, nkeep_sig)
			peakrow = -1.0
			peakcol = -1.0
			maxval = 0.0
			for i=0d,nkeep_nu0-1 do begin
				for j=0d,nkeep_sig-1 do begin
					test = where(a1mc_save ge (xhist_nu0_trim[i]-xsize_nu0/2.0) AND a1mc_save le (xhist_nu0_trim[i]+xsize_nu0/2.0) AND a2mc_save ge (xhist_sig_trim[j]-xsize_sig/2.0) AND a2mc_save le (xhist_sig_trim[j]+xsize_sig/2.0), ntest)
					cont4z[i,j] = double(ntest)
					if ntest gt maxval then begin
						maxval = double(ntest)
						peakrow = j
						peakcol = i
					endif
				endfor
			endfor
			interp_nu0_1sig_min = interpol(transpose(cont4z[*,peakrow]), xhist_nu0_trim, ffp_nu0[1]-ffp_nu0[2])
			interp_nu0_1sig_plus = interpol(transpose(cont4z[*,peakrow]), xhist_nu0_trim, ffp_nu0[1]+ffp_nu0[2])
			interp_sig_1sig_min = interpol(transpose(cont4z[peakcol,*]), xhist_sig_trim, ffp_sig[1]-ffp_sig[2])
			interp_sig_1sig_plus = interpol(transpose(cont4z[peakcol,*]), xhist_sig_trim, ffp_sig[1]+ffp_sig[2])
			thresh_1sig = mean([interp_nu0_1sig_min, interp_nu0_1sig_plus, interp_sig_1sig_min, interp_sig_1sig_plus])

			interp_nu0_2sig_min = interpol(transpose(cont4z[*,peakrow]), xhist_nu0_trim, ffp_nu0[1]-2.0*ffp_nu0[2])
			interp_nu0_2sig_plus = interpol(transpose(cont4z[*,peakrow]), xhist_nu0_trim, ffp_nu0[1]+2.0*ffp_nu0[2])
			interp_sig_2sig_min = interpol(transpose(cont4z[peakcol,*]), xhist_sig_trim, ffp_sig[1]-2.0*ffp_sig[2])
			interp_sig_2sig_plus = interpol(transpose(cont4z[peakcol,*]), xhist_sig_trim, ffp_sig[1]+2.0*ffp_sig[2])
			thresh_2sig = mean([interp_nu0_2sig_min, interp_nu0_2sig_plus, interp_sig_2sig_min, interp_sig_2sig_plus])

			contour, transpose(cont4z), xhist_sig_trim, xhist_nu0_trim, xst=1, yst=1, xr=[min(xhist_sig_trim), max(xhist_sig_trim)], yr=[min(xhist_nu0_trim), max(xhist_nu0_trim)], xtitle=textoidl('Std. Deviation (GHz)'), ytitle=textoidl('\nu_{CO} (GHz)'), charsize=1.3, charthick=1.3, nlevels=2, levels=[thresh_2sig, thresh_1sig]
			plotsym, 3, 1.0, /fill
			plots, ffp_sig[1], ffp_nu0[1], psym=8, color=70

;;;; next comes nu0 vs DC Offset
			cont5z = dblarr(nkeep_nu0, nkeep_off)
			peakrow = -1.0
			peakcol = -1.0
			maxval = 0.0
			for i=0d,nkeep_nu0-1 do begin
				for j=0d,nkeep_off-1 do begin
					test = where(a1mc_save ge (xhist_nu0_trim[i]-xsize_nu0/2.0) AND a1mc_save le (xhist_nu0_trim[i]+xsize_nu0/2.0) AND a3mc_save ge (xhist_off_trim[j]-xsize_off/2.0) AND a3mc_save le (xhist_off_trim[j]+xsize_off/2.0), ntest)
					cont5z[i,j] = double(ntest)
					if ntest gt maxval then begin
						maxval = double(ntest)
						peakrow = j
						peakcol = i
					endif
				endfor
			endfor
			interp_nu0_1sig_min = interpol(transpose(cont5z[*,peakrow]), xhist_nu0_trim, ffp_nu0[1]-ffp_nu0[2])
			interp_nu0_1sig_plus = interpol(transpose(cont5z[*,peakrow]), xhist_nu0_trim, ffp_nu0[1]+ffp_nu0[2])
			interp_off_1sig_min = interpol(transpose(cont5z[peakcol,*]), xhist_off_trim, ffp_off[1]-ffp_off[2])
			interp_off_1sig_plus = interpol(transpose(cont5z[peakcol,*]), xhist_off_trim, ffp_off[1]+ffp_off[2])
			thresh_1sig = mean([interp_nu0_1sig_min, interp_nu0_1sig_plus, interp_off_1sig_min, interp_off_1sig_plus])

			interp_nu0_2sig_min = interpol(transpose(cont5z[*,peakrow]), xhist_nu0_trim, ffp_nu0[1]-2.0*ffp_nu0[2])
			interp_nu0_2sig_plus = interpol(transpose(cont5z[*,peakrow]), xhist_nu0_trim, ffp_nu0[1]+2.0*ffp_nu0[2])
			interp_off_2sig_min = interpol(transpose(cont5z[peakcol,*]), xhist_off_trim, ffp_off[1]-2.0*ffp_off[2])
			interp_off_2sig_plus = interpol(transpose(cont5z[peakcol,*]), xhist_off_trim, ffp_off[1]+2.0*ffp_off[2])
			thresh_2sig = mean([interp_nu0_2sig_min, interp_nu0_2sig_plus, interp_off_2sig_min, interp_off_2sig_plus])

			contour, transpose(cont5z), xhist_off_trim, xhist_nu0_trim, xst=1, yst=1, xr=[min(xhist_off_trim), max(xhist_off_trim)], yr=[min(xhist_nu0_trim), max(xhist_nu0_trim)], xtitle=textoidl('DC Offset (mK)'), ytitle=textoidl('\nu_{CO} (GHz)'), charsize=1.3, charthick=1.3, nlevels=2, levels=[thresh_2sig, thresh_1sig]
			plotsym, 3, 1.0, /fill
			plots, ffp_off[1], ffp_nu0[1], psym=8, color=70

;;;; Finally, it's sigma vs DC Offset
			cont6z = dblarr(nkeep_sig, nkeep_off)
			peakrow = -1.0
			peakcol = -1.0
			maxval = 0.0
			for i=0d,nkeep_sig-1 do begin
				for j=0d,nkeep_off-1 do begin
					test = where(a2mc_save ge (xhist_sig_trim[i]-xsize_sig/2.0) AND a2mc_save le (xhist_sig_trim[i]+xsize_sig/2.0) AND a3mc_save ge (xhist_off_trim[j]-xsize_off/2.0) AND a3mc_save le (xhist_off_trim[j]+xsize_off/2.0), ntest)
					cont6z[i,j] = double(ntest)
					if ntest gt maxval then begin
						maxval = double(ntest)
						peakrow = j
						peakcol = i
					endif
				endfor
			endfor
			interp_sig_1sig_min = interpol(transpose(cont6z[*,peakrow]), xhist_sig_trim, ffp_sig[1]-ffp_sig[2])
			interp_sig_1sig_plus = interpol(transpose(cont6z[*,peakrow]), xhist_sig_trim, ffp_sig[1]+ffp_sig[2])
			interp_off_1sig_min = interpol(transpose(cont6z[peakcol,*]), xhist_off_trim, ffp_off[1]-ffp_off[2])
			interp_off_1sig_plus = interpol(transpose(cont6z[peakcol,*]), xhist_off_trim, ffp_off[1]+ffp_off[2])
			thresh_1sig = mean([interp_sig_1sig_min, interp_sig_1sig_plus, interp_off_1sig_min, interp_off_1sig_plus])

			interp_sig_2sig_min = interpol(transpose(cont6z[*,peakrow]), xhist_sig_trim, ffp_sig[1]-2.0*ffp_sig[2])
			interp_sig_2sig_plus = interpol(transpose(cont6z[*,peakrow]), xhist_sig_trim, ffp_sig[1]+2.0*ffp_sig[2])
			interp_off_2sig_min = interpol(transpose(cont6z[peakcol,*]), xhist_off_trim, ffp_off[1]-2.0*ffp_off[2])
			interp_off_2sig_plus = interpol(transpose(cont6z[peakcol,*]), xhist_off_trim, ffp_off[1]+2.0*ffp_off[2])
			thresh_2sig = mean([interp_sig_2sig_min, interp_sig_2sig_plus, interp_off_2sig_min, interp_off_2sig_plus])

			contour, transpose(cont6z), xhist_off_trim, xhist_sig_trim, xst=1, yst=1, xr=[min(xhist_off_trim), max(xhist_off_trim)], yr=[min(xhist_sig_trim), max(xhist_sig_trim)], xtitle=textoidl('DC Offset (mK)'), ytitle=textoidl('Std. Deviation (GHz)'), charsize=1.3, charthick=1.3, nlevels=2, levels=[thresh_2sig, thresh_1sig]
			plotsym, 3, 1.0, /fill
			plots, ffp_off[1], ffp_sig[1], psym=8, color=70
			device, /close
		endif else begin
		if not(keyword_set(noplots)) then begin
			window, 5, retain=2, xsize=900, ysize=600
			amp_hist_keep = where(xhist_amp ge (ffp_amp[1]-3.0*ffp_amp[2]) AND xhist_amp le (ffp_amp[1]+3.0*ffp_amp[2]),nkeep_amp)
			nu0_hist_keep = where(xhist_nu0 ge (ffp_nu0[1]-3.0*ffp_nu0[2]) AND xhist_nu0 le (ffp_nu0[1]+3.0*ffp_nu0[2]),nkeep_nu0)
			sig_hist_keep = where(xhist_sig ge (ffp_sig[1]-3.0*ffp_sig[2]) AND xhist_sig le (ffp_sig[1]+3.0*ffp_sig[2]),nkeep_sig)
			off_hist_keep = where(xhist_off ge (ffp_off[1]-3.0*ffp_off[2]) AND xhist_off le (ffp_off[1]+3.0*ffp_off[2]),nkeep_off)

			xhist_amp_trim = xhist_amp[amp_hist_keep]
			yhist_amp_trim = yhist_amp[amp_hist_keep]
			xhist_nu0_trim = xhist_nu0[nu0_hist_keep]
			yhist_nu0_trim = yhist_nu0[nu0_hist_keep]
			xhist_sig_trim = xhist_sig[sig_hist_keep]
			yhist_sig_trim = yhist_sig[sig_hist_keep]
			xhist_off_trim = xhist_off[off_hist_keep]
			yhist_off_trim = yhist_off[off_hist_keep]

;;;; First, amplitude vs nu0
			cont1z = dblarr(nkeep_amp, nkeep_nu0)
			peakrow = -1.0
			peakcol = -1.0
			maxval = 0.0
			for i=0d,nkeep_amp-1 do begin
				for j=0d,nkeep_nu0-1 do begin
					test = where(a0mc_save ge (xhist_amp_trim[i]-xsize_amp/2.0) AND a0mc_save le (xhist_amp_trim[i]+xsize_amp/2.0) AND a1mc_save ge (xhist_nu0_trim[j]-xsize_nu0/2.0) AND a1mc_save le (xhist_nu0_trim[j]+xsize_nu0/2.0), ntest)
					cont1z[i,j] = double(ntest)
					if ntest gt maxval then begin
						maxval = double(ntest)
						peakrow = j
						peakcol = i
					endif
				endfor
			endfor
			interp_amp_1sig_min = interpol(transpose(cont1z[*,peakrow]), xhist_amp_trim, ffp_amp[1]-ffp_amp[2])
			interp_amp_1sig_plus = interpol(transpose(cont1z[*,peakrow]), xhist_amp_trim, ffp_amp[1]+ffp_amp[2])
			interp_nu0_1sig_min = interpol(transpose(cont1z[peakcol,*]), xhist_nu0_trim, ffp_nu0[1]-ffp_nu0[2])
			interp_nu0_1sig_plus = interpol(transpose(cont1z[peakcol,*]), xhist_nu0_trim, ffp_nu0[1]+ffp_nu0[2])
			thresh_1sig = mean([interp_amp_1sig_min, interp_amp_1sig_plus, interp_nu0_1sig_min, interp_nu0_1sig_plus])

			interp_amp_2sig_min = interpol(transpose(cont1z[*,peakrow]), xhist_amp_trim, ffp_amp[1]-2.0*ffp_amp[2])
			interp_amp_2sig_plus = interpol(transpose(cont1z[*,peakrow]), xhist_amp_trim, ffp_amp[1]+2.0*ffp_amp[2])
			interp_nu0_2sig_min = interpol(transpose(cont1z[peakcol,*]), xhist_nu0_trim, ffp_nu0[1]-2.0*ffp_nu0[2])
			interp_nu0_2sig_plus = interpol(transpose(cont1z[peakcol,*]), xhist_nu0_trim, ffp_nu0[1]+2.0*ffp_nu0[2])
			thresh_2sig = mean([interp_amp_2sig_min, interp_amp_2sig_plus, interp_nu0_2sig_min, interp_nu0_2sig_plus])

			contour, transpose(cont1z), xhist_nu0_trim, xhist_amp_trim, xst=1, yst=1, xr=[min(xhist_nu0_trim), max(xhist_nu0_trim)], yr=[min(xhist_amp_trim), max(xhist_amp_trim)], xtitle=textoidl('\nu_{CO} (GHz)'), ytitle='Amplitude (mK)', charsize=1.3, charthick=1.3, nlevels=2, levels=[thresh_2sig, thresh_1sig]
			plotsym, 3, 1.0, /fill
			plots, ffp_nu0[1], ffp_amp[1], psym=8, color=70

;;;; next comes amplitude vs sigma
			cont2z = dblarr(nkeep_amp, nkeep_sig)
			peakrow = -1.0
			peakcol = -1.0
			maxval = 0.0
			for i=0d,nkeep_amp-1 do begin
				for j=0d,nkeep_sig-1 do begin
					test = where(a0mc_save ge (xhist_amp_trim[i]-xsize_amp/2.0) AND a0mc_save le (xhist_amp_trim[i]+xsize_amp/2.0) AND a2mc_save ge (xhist_sig_trim[j]-xsize_sig/2.0) AND a2mc_save le (xhist_sig_trim[j]+xsize_sig/2.0), ntest)
					cont2z[i,j] = double(ntest)
					if ntest gt maxval then begin
						maxval = double(ntest)
						peakrow = j
						peakcol = i
					endif
				endfor
			endfor
			interp_amp_1sig_min = interpol(transpose(cont2z[*,peakrow]), xhist_amp_trim, ffp_amp[1]-ffp_amp[2])
			interp_amp_1sig_plus = interpol(transpose(cont2z[*,peakrow]), xhist_amp_trim, ffp_amp[1]+ffp_amp[2])
			interp_sig_1sig_min = interpol(transpose(cont2z[peakcol,*]), xhist_sig_trim, ffp_sig[1]-ffp_sig[2])
			interp_sig_1sig_plus = interpol(transpose(cont2z[peakcol,*]), xhist_sig_trim, ffp_sig[1]+ffp_sig[2])
			thresh_1sig = mean([interp_amp_1sig_min, interp_amp_1sig_plus, interp_sig_1sig_min, interp_sig_1sig_plus])

			interp_amp_2sig_min = interpol(transpose(cont2z[*,peakrow]), xhist_amp_trim, ffp_amp[1]-2.0*ffp_amp[2])
			interp_amp_2sig_plus = interpol(transpose(cont2z[*,peakrow]), xhist_amp_trim, ffp_amp[1]+2.0*ffp_amp[2])
			interp_sig_2sig_min = interpol(transpose(cont2z[peakcol,*]), xhist_sig_trim, ffp_sig[1]-2.0*ffp_sig[2])
			interp_sig_2sig_plus = interpol(transpose(cont2z[peakcol,*]), xhist_sig_trim, ffp_sig[1]+2.0*ffp_sig[2])
			thresh_2sig = mean([interp_amp_2sig_min, interp_amp_2sig_plus, interp_sig_2sig_min, interp_sig_2sig_plus])

			contour, transpose(cont2z), xhist_sig_trim, xhist_amp_trim, xst=1, yst=1, xr=[min(xhist_sig_trim), max(xhist_sig_trim)], yr=[min(xhist_amp_trim), max(xhist_amp_trim)], xtitle=textoidl('Std. Deviation (GHz)'), ytitle='Amplitude (mK)', charsize=1.3, charthick=1.3, nlevels=2, levels=[thresh_2sig, thresh_1sig]
			plotsym, 3, 1.0, /fill
			plots, ffp_sig[1], ffp_amp[1], psym=8, color=70

;;;; next comes amplitude vs DC Offset
			cont3z = dblarr(nkeep_amp, nkeep_off)
			peakrow = -1.0
			peakcol = -1.0
			maxval = 0.0
			for i=0d,nkeep_amp-1 do begin
				for j=0d,nkeep_off-1 do begin
					test = where(a0mc_save ge (xhist_amp_trim[i]-xsize_amp/2.0) AND a0mc_save le (xhist_amp_trim[i]+xsize_amp/2.0) AND a3mc_save ge (xhist_off_trim[j]-xsize_off/2.0) AND a3mc_save le (xhist_off_trim[j]+xsize_off/2.0), ntest)
					cont3z[i,j] = double(ntest)
					if ntest gt maxval then begin
						maxval = double(ntest)
						peakrow = j
						peakcol = i
					endif
				endfor
			endfor
			interp_amp_1sig_min = interpol(transpose(cont3z[*,peakrow]), xhist_amp_trim, ffp_amp[1]-ffp_amp[2])
			interp_amp_1sig_plus = interpol(transpose(cont3z[*,peakrow]), xhist_amp_trim, ffp_amp[1]+ffp_amp[2])
			interp_off_1sig_min = interpol(transpose(cont3z[peakcol,*]), xhist_off_trim, ffp_off[1]-ffp_off[2])
			interp_off_1sig_plus = interpol(transpose(cont3z[peakcol,*]), xhist_off_trim, ffp_off[1]+ffp_off[2])
			thresh_1sig = mean([interp_amp_1sig_min, interp_amp_1sig_plus, interp_off_1sig_min, interp_off_1sig_plus])

			interp_amp_2sig_min = interpol(transpose(cont3z[*,peakrow]), xhist_amp_trim, ffp_amp[1]-2.0*ffp_amp[2])
			interp_amp_2sig_plus = interpol(transpose(cont3z[*,peakrow]), xhist_amp_trim, ffp_amp[1]+2.0*ffp_amp[2])
			interp_off_2sig_min = interpol(transpose(cont3z[peakcol,*]), xhist_off_trim, ffp_off[1]-2.0*ffp_off[2])
			interp_off_2sig_plus = interpol(transpose(cont3z[peakcol,*]), xhist_off_trim, ffp_off[1]+2.0*ffp_off[2])
			thresh_2sig = mean([interp_amp_2sig_min, interp_amp_2sig_plus, interp_off_2sig_min, interp_off_2sig_plus])

			contour, transpose(cont3z), xhist_off_trim, xhist_amp_trim, xst=1, yst=1, xr=[min(xhist_off_trim), max(xhist_off_trim)], yr=[min(xhist_amp_trim), max(xhist_amp_trim)], xtitle=textoidl('DC Offset (mK)'), ytitle='Amplitude (mK)', charsize=1.3, charthick=1.3, nlevels=2, levels=[thresh_2sig, thresh_1sig]
			plotsym, 3, 1.0, /fill
			plots, ffp_off[1], ffp_amp[1], psym=8, color=70

;;;; next comes nu0 vs sigma
			cont4z = dblarr(nkeep_nu0, nkeep_sig)
			peakrow = -1.0
			peakcol = -1.0
			maxval = 0.0
			for i=0d,nkeep_nu0-1 do begin
				for j=0d,nkeep_sig-1 do begin
					test = where(a1mc_save ge (xhist_nu0_trim[i]-xsize_nu0/2.0) AND a1mc_save le (xhist_nu0_trim[i]+xsize_nu0/2.0) AND a2mc_save ge (xhist_sig_trim[j]-xsize_sig/2.0) AND a2mc_save le (xhist_sig_trim[j]+xsize_sig/2.0), ntest)
					cont4z[i,j] = double(ntest)
					if ntest gt maxval then begin
						maxval = double(ntest)
						peakrow = j
						peakcol = i
					endif
				endfor
			endfor
			interp_nu0_1sig_min = interpol(transpose(cont4z[*,peakrow]), xhist_nu0_trim, ffp_nu0[1]-ffp_nu0[2])
			interp_nu0_1sig_plus = interpol(transpose(cont4z[*,peakrow]), xhist_nu0_trim, ffp_nu0[1]+ffp_nu0[2])
			interp_sig_1sig_min = interpol(transpose(cont4z[peakcol,*]), xhist_sig_trim, ffp_sig[1]-ffp_sig[2])
			interp_sig_1sig_plus = interpol(transpose(cont4z[peakcol,*]), xhist_sig_trim, ffp_sig[1]+ffp_sig[2])
			thresh_1sig = mean([interp_nu0_1sig_min, interp_nu0_1sig_plus, interp_sig_1sig_min, interp_sig_1sig_plus])

			interp_nu0_2sig_min = interpol(transpose(cont4z[*,peakrow]), xhist_nu0_trim, ffp_nu0[1]-2.0*ffp_nu0[2])
			interp_nu0_2sig_plus = interpol(transpose(cont4z[*,peakrow]), xhist_nu0_trim, ffp_nu0[1]+2.0*ffp_nu0[2])
			interp_sig_2sig_min = interpol(transpose(cont4z[peakcol,*]), xhist_sig_trim, ffp_sig[1]-2.0*ffp_sig[2])
			interp_sig_2sig_plus = interpol(transpose(cont4z[peakcol,*]), xhist_sig_trim, ffp_sig[1]+2.0*ffp_sig[2])
			thresh_2sig = mean([interp_nu0_2sig_min, interp_nu0_2sig_plus, interp_sig_2sig_min, interp_sig_2sig_plus])

			contour, transpose(cont4z), xhist_sig_trim, xhist_nu0_trim, xst=1, yst=1, xr=[min(xhist_sig_trim), max(xhist_sig_trim)], yr=[min(xhist_nu0_trim), max(xhist_nu0_trim)], xtitle=textoidl('Std. Deviation (GHz)'), ytitle=textoidl('\nu_{CO} (GHz)'), charsize=1.3, charthick=1.3, nlevels=2, levels=[thresh_2sig, thresh_1sig]
			plotsym, 3, 1.0, /fill
			plots, ffp_sig[1], ffp_nu0[1], psym=8, color=70

;;;; next comes nu0 vs DC Offset
			cont5z = dblarr(nkeep_nu0, nkeep_off)
			peakrow = -1.0
			peakcol = -1.0
			maxval = 0.0
			for i=0d,nkeep_nu0-1 do begin
				for j=0d,nkeep_off-1 do begin
					test = where(a1mc_save ge (xhist_nu0_trim[i]-xsize_nu0/2.0) AND a1mc_save le (xhist_nu0_trim[i]+xsize_nu0/2.0) AND a3mc_save ge (xhist_off_trim[j]-xsize_off/2.0) AND a3mc_save le (xhist_off_trim[j]+xsize_off/2.0), ntest)
					cont5z[i,j] = double(ntest)
					if ntest gt maxval then begin
						maxval = double(ntest)
						peakrow = j
						peakcol = i
					endif
				endfor
			endfor
			interp_nu0_1sig_min = interpol(transpose(cont5z[*,peakrow]), xhist_nu0_trim, ffp_nu0[1]-ffp_nu0[2])
			interp_nu0_1sig_plus = interpol(transpose(cont5z[*,peakrow]), xhist_nu0_trim, ffp_nu0[1]+ffp_nu0[2])
			interp_off_1sig_min = interpol(transpose(cont5z[peakcol,*]), xhist_off_trim, ffp_off[1]-ffp_off[2])
			interp_off_1sig_plus = interpol(transpose(cont5z[peakcol,*]), xhist_off_trim, ffp_off[1]+ffp_off[2])
			thresh_1sig = mean([interp_nu0_1sig_min, interp_nu0_1sig_plus, interp_off_1sig_min, interp_off_1sig_plus])

			interp_nu0_2sig_min = interpol(transpose(cont5z[*,peakrow]), xhist_nu0_trim, ffp_nu0[1]-2.0*ffp_nu0[2])
			interp_nu0_2sig_plus = interpol(transpose(cont5z[*,peakrow]), xhist_nu0_trim, ffp_nu0[1]+2.0*ffp_nu0[2])
			interp_off_2sig_min = interpol(transpose(cont5z[peakcol,*]), xhist_off_trim, ffp_off[1]-2.0*ffp_off[2])
			interp_off_2sig_plus = interpol(transpose(cont5z[peakcol,*]), xhist_off_trim, ffp_off[1]+2.0*ffp_off[2])
			thresh_2sig = mean([interp_nu0_2sig_min, interp_nu0_2sig_plus, interp_off_2sig_min, interp_off_2sig_plus])

			contour, transpose(cont5z), xhist_off_trim, xhist_nu0_trim, xst=1, yst=1, xr=[min(xhist_off_trim), max(xhist_off_trim)], yr=[min(xhist_nu0_trim), max(xhist_nu0_trim)], xtitle=textoidl('DC Offset (mK)'), ytitle=textoidl('\nu_{CO} (GHz)'), charsize=1.3, charthick=1.3, nlevels=2, levels=[thresh_2sig, thresh_1sig]
			plotsym, 3, 1.0, /fill
			plots, ffp_off[1], ffp_nu0[1], psym=8, color=70

;;;; Finally, it's sigma vs DC Offset
			cont6z = dblarr(nkeep_sig, nkeep_off)
			peakrow = -1.0
			peakcol = -1.0
			maxval = 0.0
			for i=0d,nkeep_sig-1 do begin
				for j=0d,nkeep_off-1 do begin
					test = where(a2mc_save ge (xhist_sig_trim[i]-xsize_sig/2.0) AND a2mc_save le (xhist_sig_trim[i]+xsize_sig/2.0) AND a3mc_save ge (xhist_off_trim[j]-xsize_off/2.0) AND a3mc_save le (xhist_off_trim[j]+xsize_off/2.0), ntest)
					cont6z[i,j] = double(ntest)
					if ntest gt maxval then begin
						maxval = double(ntest)
						peakrow = j
						peakcol = i
					endif
				endfor
			endfor
			interp_sig_1sig_min = interpol(transpose(cont6z[*,peakrow]), xhist_sig_trim, ffp_sig[1]-ffp_sig[2])
			interp_sig_1sig_plus = interpol(transpose(cont6z[*,peakrow]), xhist_sig_trim, ffp_sig[1]+ffp_sig[2])
			interp_off_1sig_min = interpol(transpose(cont6z[peakcol,*]), xhist_off_trim, ffp_off[1]-ffp_off[2])
			interp_off_1sig_plus = interpol(transpose(cont6z[peakcol,*]), xhist_off_trim, ffp_off[1]+ffp_off[2])
			thresh_1sig = mean([interp_sig_1sig_min, interp_sig_1sig_plus, interp_off_1sig_min, interp_off_1sig_plus])

			interp_sig_2sig_min = interpol(transpose(cont6z[*,peakrow]), xhist_sig_trim, ffp_sig[1]-2.0*ffp_sig[2])
			interp_sig_2sig_plus = interpol(transpose(cont6z[*,peakrow]), xhist_sig_trim, ffp_sig[1]+2.0*ffp_sig[2])
			interp_off_2sig_min = interpol(transpose(cont6z[peakcol,*]), xhist_off_trim, ffp_off[1]-2.0*ffp_off[2])
			interp_off_2sig_plus = interpol(transpose(cont6z[peakcol,*]), xhist_off_trim, ffp_off[1]+2.0*ffp_off[2])
			thresh_2sig = mean([interp_sig_2sig_min, interp_sig_2sig_plus, interp_off_2sig_min, interp_off_2sig_plus])

			contour, transpose(cont6z), xhist_off_trim, xhist_sig_trim, xst=1, yst=1, xr=[min(xhist_off_trim), max(xhist_off_trim)], yr=[min(xhist_sig_trim), max(xhist_sig_trim)], xtitle=textoidl('DC Offset (mK)'), ytitle=textoidl('Std. Deviation (GHz)'), charsize=1.3, charthick=1.3, nlevels=2, levels=[thresh_2sig, thresh_1sig]
			plotsym, 3, 1.0, /fill
			plots, ffp_off[1], ffp_sig[1], psym=8, color=70
		endif
		endelse
	endif
	!P.Multi = 0
endif

if(keyword_set(savsetname) AND z eq 1) then begin
	if(keyword_set(save_mcmc_steps)) then begin
		save, filename=savsetname, x, y, x_raw, y_raw, first_guess_gaussian, first_bestfit_gaussian, base_poly_ord, x_baseline, baseline, ffp_amp, ffp_nu0, ffp_sig, ffp_fwhm, ffp_off, fwindow, normrange, frange, sigrange, offsetrange, g0, nsteps, nburnin, a0mc, a1mc, a2mc, a3mc, a0mc_save, a1mc_save, a2mc_save, a3mc_save, overall_best_gaussian, rms, linetotal, sigmaline, SNR, freq_full, y_full
	endif else begin
		save, filename=savsetname, x, y, x_raw, y_raw, first_guess_gaussian, first_bestfit_gaussian, base_poly_ord, x_baseline, baseline, ffp_amp, ffp_nu0, ffp_sig, ffp_fwhm, ffp_off, fwindow, normrange, frange, sigrange, offsetrange, g0, nsteps, nburnin, overall_best_gaussian, rms, linetotal, sigmaline, SNR, freq_full, y_full
	endelse
endif

endfor


END
