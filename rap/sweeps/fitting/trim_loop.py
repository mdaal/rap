def trim_loop(self,N = 20,Verbose = True,):
	import numpy.ma as ma
	f = f1 = ma.array(self.loop.freq)
	z = z1 = ma.array(self.loop.z)
	# estimate resonant freq using resonance dip
	zr_mag_est = np.abs(z).min()
	zr_est_index = np.where(np.abs(z)==zr_mag_est)[0][0]

	# estimate max transmission mag using max valuse of abs(z)
	z_max_mag = np.abs(z).max()

	#Depth of resonance in dB
	depth_est = 20.0*np.log10(zr_mag_est/z_max_mag)

	#Magnitude of resonance dip at half max
	res_half_max_mag = (z_max_mag+zr_mag_est)/2

	#find the indices of the closest points to this magnitude along the loop, one below zr_mag_est and one above zr_mag_est
	a = np.square(np.abs(z[:zr_est_index+1]) - res_half_max_mag)
	lower_index = np.argmin(a)
	a = np.square(np.abs(z[zr_est_index:]) - res_half_max_mag)
	upper_index = np.argmin(a) + zr_est_index

	#estimate the FWHM bandwidth of the resonance
	f_upper_FWHM = f[upper_index]
	f_lower_FWHM = f[lower_index]
	FWHM_est = np.abs(f_upper_FWHM - f_lower_FWHM)
	fr_est = f[zr_est_index]

	#Bandwidth Cut: cut data that is more than N * FWHM_est away from zr_mag_est
	z = z2 = ma.masked_where((f > fr_est + N*FWHM_est) | (fr_est - N*FWHM_est > f),z)
	f = f2 = ma.array(f,mask = z.mask)

	self.loop.z = ma.compressed(z)
	self.loop.freq = ma.compressed(f)

	if Verbose: 
		print('Bandwidth cut:\n\t{1} points outside of fr_est +/- {0}*FWHM_est removed, {2} remaining data points'.format(N, *self._points_removed(z1,z2)))
