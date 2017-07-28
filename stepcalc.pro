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
