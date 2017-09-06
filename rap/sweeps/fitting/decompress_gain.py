from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

def decompress_gain(Sweep_Array, loop, metadata,Compression_Calibration_Index = -1, Show_Plot = True, Verbose = True):
	''' Assumes the two lowest input powers of the power sweep are not gain compressed, thus
	cannot be used if the two lowest powers are gain compressed. '''


	Sweep_Array_Record_Index = loop.index 
	V = Sweep_Array['Heater_Voltage'][Sweep_Array_Record_Index]
	Fs = Sweep_Array['Fstart'][Sweep_Array_Record_Index]
	P = Sweep_Array['Pinput_dB'][Sweep_Array_Record_Index]

	Sweep_Array = np.extract((Sweep_Array['Heater_Voltage'] == V) & ( Sweep_Array['Fstart']==Fs) , Sweep_Array)


	num_sweep_powers = Sweep_Array['Pinput_dB'].shape[0]

	if num_sweep_powers <= 4:
		print('Number of sweep powers, {0}, is insufficient to perform gain decompression.'.format(num_sweep_powers))
		return
	#else:
	#	print('Performing gain decompression on {0} sweep powers.'.format(num_sweep_powers))

	Pin = np.power(10, Sweep_Array['Pinput_dB']/10.0) #mW, Probe Power

	#ChooseCompression calobration data from Power Sweep Data. 
	#It is the S21(Compression_Calibration_Index) for every sweep power 
	compression_calibration_data = np.power(np.abs(Sweep_Array['S21'][:,Compression_Calibration_Index]),2) #Pout/Pin,  
	# alternatively : np.average(Sweep_Array['S21'][:,Compression_Calibration_Index:Compression_Calibration_Index+n],axis = 1) #average over  n freq points.
	Pout = compression_calibration_data*Pin 



	# calculated_power_gain is power gain calculated from the slope of the two smallest input powers in Pin
	values, indices = np.unique(Pin, return_index=True)
	min_index,min_plus_index =  indices[:2]   
	# When Pin = 0, 0 != Pout = Pin*gaain. There is an offset, i.e. a y-intercept, b, such at y = m*x+b. Next, we find m.  
	calculated_power_gain = (Pout[min_plus_index] - Pout[min_index])/(Pin[min_plus_index ]-Pin[min_index]) 

	#Pout_ideal is the output power assuming linear gain
	Pout_ideal = lambda p_in: calculated_power_gain*(p_in-Pin[0]) + Pout[0]

	Probe_Power_Mag = np.power(10,Sweep_Array[Sweep_Array_Record_Index]['Pinput_dB']/10) #-- Substitute for input power
	S21 = Sweep_Array[Sweep_Array_Record_Index]['S21']
	S21_Pout = np.power(np.abs(S21),2)*Probe_Power_Mag

	# create interpolation funcation to what Pin would be at an arbitrary Pout
	decompression_function = interp1d(Pout,Pin,kind = 'linear')

	# for polynomial to Pout vs Pin curve and use this to extrapolate values where Pout in not in interpolation domain
	def decompression_function_fit(pout, a,b,c):
		return a*np.power(pout,2)+b*pout+c
	popt,pcov = curve_fit(decompression_function_fit, Pout, Pin)
	decompression_function_extrap = lambda pout : decompression_function_fit(pout,popt[0],popt[1],popt[2])

	
	def decompress_element(z):
		z_Pout = np.power(np.abs(z),2)*Probe_Power_Mag
		if z_Pout <= Pout.min(): #Do nothinge when z_Pout is less than the interpolation range, Pout.min() to Pout.max()
			return z
		elif Pout.min() < z_Pout < Pout.max(): # Interpolate to find ideal Pout (assuming linear gain) when z_Pout is in interpolation domain 
			return z*np.sqrt(Pout_ideal(decompression_function(z_Pout))/Probe_Power_Mag)/np.abs(z)
		else: # Pout.max() <= z_Pout --  Extrapolate to find ideal Pout when z_Pout is above interpolation domain
			return z*np.sqrt(Pout_ideal(decompression_function_extrap(z_Pout))/Probe_Power_Mag)/np.abs(z)

	decompress_array = np.vectorize(decompress_element) # Vectorize for speed

	loop.z = S21_Decompressed = decompress_array(S21)

	if Verbose == True:
		print('Gain decompression calculation is based on {0} sweep powers.'.format(num_sweep_powers))
		print('Power out at zero input power is {0} mW'.format(calculated_power_gain*(0-Pin[0]) + Pout[0]))

	if Show_Plot:
		fig1           = plt.figure(figsize = (15,5))
		Freq           = Sweep_Array[Sweep_Array_Record_Index]['Frequencies']
		#majorFormatter = FormatStrFormatter('%d')
		majormaxnlocator    = MaxNLocator(nbins = 5)
		minormaxnlocator    = MaxNLocator(nbins = 5*5)
		#minorLocator   = MultipleLocator((Freq.max() - Freq.min())/25)
		

		ax1 = fig1.add_subplot(131)
		ax1.set_xlabel('Power In [mW]')
		line1 = ax1.plot(Pin,Pout, 'b-', label = 'Measured')
		line2 = ax1.plot(Pin,Pout_ideal(Pin), 'r-', label = 'Ideal')
		ax1.set_ylabel('Power Out [mW]', color='b')
		ax1.set_title('Gain Compression', fontsize=9)
		ax1.legend(loc = 'best', fontsize=9)
		plt.setp(ax1.get_xticklabels(),rotation = 45, fontsize=9)
		ax1.grid()
		#fig1.canvas.manager.resize(800,800)

		
		ax2 = fig1.add_subplot(132, aspect='equal')
		line2 = ax2.plot(S21.real,S21.imag, color='blue', linestyle='solid', linewidth = 3, label = 'Measured') 
		line1 = ax2.plot(S21_Decompressed.real, S21_Decompressed.imag, 'g-',linewidth = 3, label = 'Corrected')
		ax2.grid()
		ax2.set_title('Resonance Loop', fontsize=9)
		plt.setp(ax2.get_xticklabels(),rotation = 45)
		#ax2.legend(loc = 'best')

		
		ax3 = fig1.add_subplot(133)
		ax3.set_xlabel('Freq [Hz]')
		line1 = ax3.plot(Freq,10*np.log10(np.abs(S21)), 'b-',label = 'Measured',linewidth = 3)
		line2 = ax3.plot(Freq,10*np.log10(np.abs(S21_Decompressed)), 'g-', label = 'Corrected',linewidth = 3)
		ax3.set_ylabel('$|S_{21}|$ [dB]', color='k')
		ax3.legend(loc = 'best', fontsize=9)
		ax3.xaxis.set_major_locator(majormaxnlocator)
		#ax3.tick_params( axis='both', labelsize=9)
		plt.setp(ax3.get_xticklabels(),rotation = 45, fontsize=9)
		#ax3.xaxis.set_major_formatter(majorFormatter)
		ax3.xaxis.set_minor_locator(minormaxnlocator)
		ax3.set_title('Resonance Dip', fontsize=9)
		ax3.grid()

		fig1.subplots_adjust(wspace = 0.6,bottom = 0.09, top = 0.1)
		fig1.suptitle('Run: {0}, Sensor: {1}, Ground Plane: {2}, Readout Power: {3} dBm, Date: {4}'.format(metadata.Run, metadata.Sensor,metadata.Ground_Plane,Sweep_Array[Sweep_Array_Record_Index]['Pinput_dB'],metadata.Time_Created), fontsize=10)
		#plt.tight_layout()
		plt.setp(fig1, tight_layout = True)
		plt.show()