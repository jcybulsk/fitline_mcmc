; Here's our Gaussian
FUNCTION thegaussian, x, a
	z = (x-a[1])/a[2]
	y = a[0]*exp(-z^2/2) + a[3]
	return, y
END
