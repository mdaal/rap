def remove_cable_delay(self, Show_Plot = True, Verbose = True, center_freq = None, Force_Recalculate = False):
	'''
	If self.metadate.Electrical_Delay is not None, then use this value as cable delay and remove from data 

	If self.metadate.Electrical_Delay is None:
	- Determine cable delay by finding delay value, tau, which minimizes distance between adjacent S21 points. 
	Then cancel out tau in S21 data and save corrected loop in self.loop.z. Set self.metadate.Electrical_Delay = tau.
	- If S21 is large array, down sample it first before performing minimization

	If self.metadate.Electrical_Delay is None and center_freq is given:
	-If center_freq is given, this function computes the electrical delay by determining the bandwidth over which the S21 
	circle completes a full revolution starting at center_freq and ending at ~ center_freq + tau^-1. Where tau is approximated
	as the vaule deterined by minimum  distance above. 
	-center_freq should only be used when S21 is is sufficiently broadband to generously cover center_freq and ~center_freq + tau^-1.
	center_freq is in Hertz.

	Return tau in any case.

	If Force_Recalculate == False Electrical delay will be recalculated and reset in metadata
	'''

	S21 = self.loop.z
	f= self.loop.freq

	j = np.complex(0,1)
	n = 1

	if (self.metadata.Electrical_Delay == None) or (Force_Recalculate == True):
		cable_delay_max = 200e-9 # Seconds - guess as to maximum value of cable delay
		cable_delay_guess  = 80e-9 # Seconds
		freq_spacing = np.abs(f[0] - f[1])
		# phase change between adjacent frequency points is 360 * tau * freq_spacing --> want tau * freq_spacing < 1 to see loop
		if (3*freq_spacing * cable_delay_max < 1) & (f.size > 3200):
			n1 = int(np.floor( 1./ (3*freq_spacing * cable_delay_max) )) #want as least 3 points per circle
			n2 = int(np.floor(f.size/3200))
			n = min(n1,n2)


		def obj(t):
			'''This objective fuction yields the sum of squared distance between adjacent (if n = 1) or n-separated S21 points.
			'''
			
			S21a = np.exp(2*np.pi*j*f[1::n]*t)*S21[1::n] # go in steps of n 
			S21b = np.exp(2*np.pi*j*f[:-1:n]*t)*S21[:-1:n] # go in steps of n 
			diff = S21a-S21b
			return (diff*diff.conjugate()).real.sum()

		

		# # Could use Nelder-Mead
		# out = minimize(obj,cable_delay_guess, method='Nelder-Mead',tol=1e-20,options={'disp':False})
		# cable_delay = out.x[0] # in seconds
		
		out = minimize(obj,cable_delay_guess, method='Powell',tol=1e-20,options={'disp':False, 'ftol':1e-14,'xtol':1e-14})
		cable_delay_min_distance = out.x.item() #in Seconds 
		cable_delay  = cable_delay_min_distance 

		if center_freq is not None:
			cable_delay_bandwidth = 1/cable_delay_min_distance #Hz - Estimate of winding bandwith using tau

			closest_index_to_center_freq = np.where(np.abs(f-center_freq) == np.abs(f-center_freq).min()) 
			s21 = S21*np.exp(np.complex(0,-np.angle(S21[closest_index_to_center_freq]))) #rotate circle so that S21[center_freq] is close to positive x axis, and angle(S21[center_freq]) ~ 0
			
			
			condition = ((center_freq - .30*cable_delay_bandwidth) < f) & (f<center_freq+.30*cable_delay_bandwidth)
			f_lower_band =np.extract(condition,f)
			s21_lower_band = np.extract(condition,s21)
			ang_lower_band = np.extract(condition,np.angle(s21)) #np.angle has range [+pi,-pi]
			interp_lower_band = interp1d(ang_lower_band, f_lower_band,kind='linear')
			lower_x_axis_crossing_freq = interp_lower_band(0).item()

			center_freq = center_freq + cable_delay_bandwidth #shift to upper band
			condition = ((center_freq - .30*cable_delay_bandwidth) < f) & (f<center_freq+.30*cable_delay_bandwidth)
			f_upper_band =np.extract(condition,f)
			s21_upper_band = np.extract(condition,s21)
			ang_upper_band = np.extract(condition,np.angle(s21)) #np.angle has range [+pi,-pi]
			interp_upper_band = interp1d(ang_upper_band, f_upper_band,kind='linear')
			upper_x_axis_crossing_freq = interp_upper_band(0).item()

			winding_bandwidth = upper_x_axis_crossing_freq - lower_x_axis_crossing_freq

			cable_delay_winding = 1/winding_bandwidth
			cable_delay = cable_delay_winding #override cable_delay_min_distance 
	else:
		cable_delay = self.metadata.Electrical_Delay
		center_freq = None

	S21_Corrected = np.exp(2*np.pi*f*j*cable_delay)*S21
	
	if Verbose == True:
		if n>1:
			print('S21 downsampled by factor n = {}.'.format(n))
		if (self.metadata.Electrical_Delay == None) or (Force_Recalculate == True):
			print('cable delay is {} seconds by minimum distance method'.format(cable_delay_min_distance))	
		else: 
			print('cable delay is {} seconds as found in metadata'.format(self.metadata.Electrical_Delay))
		if center_freq is not None:
			print('cable delay is {} seconds by loop winding method'.format(cable_delay_winding))
		
			
	if Show_Plot:
		fig = plt.figure( figsize=(9,6))#, dpi=150)
		ax = {}	
		def plot_loops(ax):
			from matplotlib.ticker import MaxNLocator
			majormaxnlocator    = MaxNLocator(nbins = 5)
			minormaxnlocator    = MaxNLocator(nbins = 5*5)
			#ax2 = fig.add_subplot(111, aspect='equal')
			line2 = ax.plot(S21.real,S21.imag, color='blue', linestyle='solid', linewidth = 3, label = 'Measured') 
			line1 = ax.plot(S21_Corrected.real, S21_Corrected.imag, 'g-',linewidth = 3, label = 'Corrected')
			ax.grid()
			ax.set_title('Resonance Loop', fontsize=9)
			plt.setp(ax.get_xticklabels(),rotation = 45)
			ax.legend(loc = 'best')

		if center_freq is None:
			gs = gridspec.GridSpec(1, 1)
			ax[1] = plt.subplot(gs[0, 0],aspect='equal')
			plot_loops(ax[1])

		else:
			gs = gridspec.GridSpec(2, 3)#,width_ratios=[2,2,1])

			ax[1] = plt.subplot(gs[:,:2],aspect='equal')
			ax[2] = plt.subplot(gs[0, 2])
			ax[3] = plt.subplot(gs[1, 2], aspect='equal' )
			plot_loops(ax[1])
			curve = ax[2].plot(f_lower_band,ang_lower_band, linestyle = '-')
			curve = ax[2].plot(f_upper_band,ang_upper_band, linestyle = '-')
			curve = ax[3].plot(s21_lower_band.real,s21_lower_band.imag, linestyle = '-')
			curve = ax[3].plot(s21_upper_band.real,s21_upper_band.imag, linestyle = '-')
			plt.setp(ax[2].get_xticklabels(),rotation = 45)	
			plt.setp(ax[3].get_xticklabels(),rotation = 45)				


		#fig.subplots_adjust(wspace = 0.6,bottom = 0.09, top = 0.1)
		#plt.setp(fig, tight_layout = True)
		plt.show()

	self.metadata.Electrical_Delay = cable_delay
	self.loop.z = S21_Corrected
	return cable_delay
