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
