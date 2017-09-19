import io
from scipy.signal import gaussian,wiener, filtfilt, butter,  freqz
from scipy.ndimage import filters
from scipy.interpolate import UnivariateSpline


class thermometry:
	def __init__(self):
		pass
	
	def load_MonitoringVI_file(self, filename, temp_list = None, process_therm = 1):
		'''Reads in thermometer data text file created by  MonitoringVI, and plots the temperature as a function of time.
		
		temp_list is a list of tuples, [(heater_voltage,temperature), ...] which are plotted on top of the temperature
		versus time points. This allows one to visually check the calibration, temp_list.

		process_therm is the column number of the thermometer whose data is processed by several filtering algorithms and
		plotted.
		'''
   

		pos = filename.rfind(os.sep)

		try:
			with io.open(filename[:pos+1]+ 'Make_ScanData.m',mode='r') as f:
				while 1:
					line  = f.readline()
					if line == '': # End of file is reached
						break
					elif line.find('ScanData.Heater_Voltage') >= 0:
						Voltages = line[line.find('['):line.find(']')+1]
						break
		except:
			print('Unable to find or read Make_ScanData.m for list of heater voltages')
			Voltages = 'Unknown'

		with io.open(filename,mode='r') as f:
			
			temp_data_header = ''
			while temp_data_header.strip() =='':
				temp_data_header = f.readline()

			therm_list = [t for t in temp_data_header.strip().split('\t')[1:] if (t.strip() != 'None') & (t.strip() != '')]
			

		temp_data = np.loadtxt(filename, dtype=np.float, comments='#', delimiter=None, converters=None, skiprows=3, usecols=None, unpack=False, ndmin=0)
		
		num_col  = temp_data.shape[1]
		start_col = 1 #index of first column in data that has thermometer data
		if process_therm > num_col - start_col:
			print('process_therm = {} exceeds number of thermometers in data. Choose an lower number. Aborting...'.format(process_therm))
			return

		# Gaussian Filter
		num_pts_in_gaussian_window = 20
		b = gaussian(num_pts_in_gaussian_window, 10)
		ga = filters.convolve1d(temp_data[:,process_therm], b/b.sum())

		# buterworth Filter
		npts = temp_data[:,process_therm].size
		end = temp_data[-1,0]
		dt = end/float(npts)
		nyf = 0.5/dt	
		b, a = butter(4, .1)#1.5/nyf)
		fl = filtfilt(b, a, temp_data[:,process_therm])
    	
    	#Spline Fit
		sp = UnivariateSpline(temp_data[:,0], temp_data[:,process_therm])

		#weiner filter
		wi = wiener(temp_data[:,process_therm], mysize=40, noise=10)

		fig1 = plt.figure( facecolor = 'w',figsize = (10,10))
		ax = fig1.add_subplot(1,1,1)		
		
		if isinstance(temp_list, list):
			for temp_tuple in temp_list:
				hline = ax.axhline(y = temp_tuple[1],linewidth=1, color='g', alpha = 0.3 ,linestyle = ':',   label = None)


		color_incr = 1.0/(num_col-start_col)
		for therm_num in xrange(start_col, num_col): # plot all thermometer data present
			line = ax.plot(temp_data[:,0], temp_data[:,therm_num],color=(0,color_incr*therm_num,0), alpha = 0.4 if therm_num != 1 else 1, linewidth = 3,label = therm_list.pop(0) if therm_list[0] != None else 'Therm{0}'.format(therm_num))

		#plot filter outputs for THE FIRST thermometer only 
		line2 = ax.plot(temp_data[:,0], ga, 'y', linewidth = 3, label = 'Gaussian Conv') # Gaussian Convolution
		line3 = ax.plot(temp_data[:,0], fl, 'c', linewidth = 3, label = 'Butterworth') # butterworth 
		line4 = ax.plot(temp_data[:,0], sp(temp_data[:,0]), 'k', linewidth = 3, label = 'Spline') # bspline
		line5 = ax.plot(temp_data[:,0], wi, 'r', linewidth = 3, label = 'Weiner') # weiner

		ax.grid(b=True, which='major', color='b', alpha = 0.2, linestyle='-')
		ax.grid(b=True, which='minor', color='b', alpha = 0.2,linestyle='--')
		ax.set_title('Heater Voltages = {}'.format(Voltages), fontsize=12)
		ax.set_ylabel('Temperature [Kelvin]')
		ax.set_xlabel('Seconds')
		ax.legend(loc = 'best', fontsize=10,scatterpoints =1, numpoints = 1, labelspacing = .1)
		plt.show()