from ..metadata import metadata

from data_management.utils import _read_scandata_from_file
from data_management.utils import _download_data

def load_scandata(metadata, _download_data, _read_scandata_from_file, file_location):
	''' file_location is the locaiton of the scandata.mat file. It can be a URL, filename or /path/filename.
	assumes that self.data is in the form of matlab ScanData Structure'''

	#delete previous metadata object
	del(metadata)
	metadata = metadata()

	if file_location.startswith('http'):
		_download_data(file_location)
	else:
		_read_scandata_from_file(file_location)

	ScanData = self.data['ScanData']
	
	# These tags specify the data to pull out of self.data['ScanData']. syntax is 
	# (field of self.data['ScanData'] to extract, self.metadata name to save to ('key:sub-key' ifself.metadata.key is a dict), 
	#			 type of value (arrays are None),optional sub-field of self.data['ScanData'] to extract)
	tags = [('Run','Run', str), ('Start_Date','Fridge_Run_Start_Date',str), ('Location','Test_Location', str), 
			('Sensor','Sensor',str), ('Ground_Plane','Ground_Plane',str), ('Box','Box',str), ('Press','Press',str), 
			('Notes','Notes',str),('Time','Time_Created',str),('Temperature','Fridge_Base_Temp',float),
			('Powers','Powers', None), ('Resolution','Min_Freq_Resolution', np.float), ('IFBW','IFBW', np.float),
			('Heater_Voltage','Heater_Voltage',None), ('Average_Factor','NA_Average_Factor', np.float), 
			('Minimum_Q','Minimum_Q', np.float), ('Range','Range',None), ('Added_Atten','Atten_Added_At_NA', np.float),  
			('Num_Points_Per_Scan','Num_Points_Per_Scan',np.float), ('Freq_Range', 'Freq_Range',None), 
			('Pause','Wait_Time',np.float), ('LNA', 'LNA:LNA', str), ('HEMT', 'LNA:Vg', str,'Vg'),
			('HEMT', 'LNA:Id', str,'Id'),  ('HEMT', 'LNA:Vd', str,'Vd'), ('Atten_4K', 'Atten_At_4K', np.float32),
			('Atten_NA_Output', 'Atten_NA_Output',np.float32), ('Atten_NA_Input','Atten_NA_Input',np.float32),
			('Atten_RTAmp_Input','Atten_RTAmp_Input',np.float32), ('RTAmp_In_Use', 'RTAmp_In_Use', int),
			('Elapsed_Time', 'Meaurement_Duration', np.float),('Thermometer_Configuration','Thermometer_Configuration',None),
			('Thermometer_Bias','Thermometer_Voltage_Bias', None)]

	for t in tags:
		try:
			if t[1].find(':')>-1: #The case of a dictionary
				t1 = t[1].split(':')

				#This try/except block is for the case where self.metadata.__dict__['?'] is a dictionary
				try:
					metadata.__dict__[t1[0]].update([(t1[1],self._extract_type(ScanData[t[0]], return_type = t[2],field = t[3] if len(t) > 3 else None))])
				except:
					metadata.__dict__[t1[0]] = dict([(t1[1],self._extract_type(ScanData[t[0]], return_type = t[2],field = t[3] if len(t) > 3 else None))])	
			else:
				metadata.__dict__[t[1]] = self._extract_type(ScanData[t[0]], return_type = t[2],field = t[3] if len(t) > 3 else None)
		except: 
			#the case that the field does not exist or that its in an unexpected format
			#print('Field named {0}{1} is not found. Setting value to None'.format(t[0], (':'+t[3]) if len(t) > 3 else '')) # this usesScandata nomenclature
			print('Field named {0} is not found.'.format(t[1])) # this uses self.metadata nomenclature
	try:
		metadata.Powers                = metadata.Powers.squeeze() #for case there are multiples powers
	except:
		metadata.Powers                = np.array([metadata.Powers]) # for the case there is only one power		

	
	# Remove nested array for Thermometer_Voltage_Bias data, if this data exists
	if hasattr(metadata,'Thermometer_Voltage_Bias'):
		metadata.Thermometer_Voltage_Bias  = metadata.Thermometer_Voltage_Bias.reshape((metadata.Thermometer_Voltage_Bias.shape[1],))

	if metadata.Thermometer_Configuration is not None:#if self.metadata.Thermometer_Configuration != None:
		metadata.Thermometer_Configuration = (str(metadata.Thermometer_Configuration.squeeze()[0][0]),str(metadata.Thermometer_Configuration.squeeze()[1][0]))

	# Reshape  Heater_Voltage array and  Remove final Heater voltage from self.Heater_Voltage (The final value is just the heater value at which to leave fridge )
	metadata.Heater_Voltage = metadata.Heater_Voltage.reshape((metadata.Heater_Voltage.shape[1],))
	metadata.Heater_Voltage = metadata.Heater_Voltage[0:-1]


	print('Loading Run: {0}'.format(metadata.Run))
	print('There are {0} heater voltage(s), {1} input power(s), and {2} frequecy span(s)'.format(metadata.Heater_Voltage.shape[0],metadata.Powers.shape[0], metadata.Freq_Range.shape[0]))
	heater_voltage_num = 0; power_sweep_num = 0; fsteps = 0; tpoints = 0;

	# determine fsteps = length of the freq/S21 array
	if metadata.Heater_Voltage.shape[0] == 1:
		fsteps = metadata.Freq_Range[heater_voltage_num][1]['PowerSweep'][0][0][power_sweep_num][2].squeeze()[()].size # non temp sweep, single freq_range, powersweep
		try: # Try to determine the number of temperture readings per scan. If data does not contain temp readings, pass
			tpoints = metadata.Freq_Range[heater_voltage_num][1]['PowerSweep'][0][0][power_sweep_num][3].squeeze()[()].size
		except:
			pass
	else:					
		for freq_range_num in xrange(metadata.Freq_Range.shape[0]):			
			steps = metadata.Freq_Range[freq_range_num][1]['Temp'][0][0][heater_voltage_num][1]['PowerSweep'][0][0][power_sweep_num][2].squeeze()[()].size
			fsteps = max(steps,fsteps)
			try: # Try to determine the number of temperture readings per scan. If data does not contain temp readings, pass
				points = metadata.Freq_Range[freq_range_num][1]['Temp'][0][0][heater_voltage_num][1]['PowerSweep'][0][0][power_sweep_num][3].squeeze()[()].size
				tpoints = max(points,tpoints)
			except:
				pass


	self._define_sweep_data_columns(fsteps,tpoints)

	
	metadata.Num_Powers 			= metadata.Powers.size
	metadata.Num_Heater_Voltages 	= metadata.Heater_Voltage.size
	metadata.Num_Ranges 			= metadata.Range.shape[0]
	try:
		metadata.Cable_Calibration = self._Cable_Calibration
		print('Cable Calibraiton data found and saved in Sweep_Array metadata.')
	except:
		pass

	try:
		metadata.Temperature_Calibration = self._Temperature_Calibration
		print('Temperature Calibraiton data found and saved in Sweep_Array metadata.')
	except:
		pass

	if metadata.Num_Temperatures > 0:
		print('Temperature readings found for scan(s). {0} readings per scan'.format(metadata.Num_Temperatures))
	### Examples of dealing with Freq_Range Data structure imported from Matlab .mat file			
	#    k.Freq_Range[heater_voltage_num][1]['PowerSweep']
	# 					  k.Freq_Range[0][1]['PowerSweep']
	# j.Freq_Range[0][1]['Temp'][0][0][0][1]['PowerSweep']
	# dt = np.dtype(('O', (2,3)))
	# entended = np.zeros(0,dtype = dt)
	# dt = np.dtype(('O',('O',[('Temp',('O',('O')))])))
	# dt = np.dtype([('O',[('O',[('Temp',[('O',('O'))])])])])
	# #this is the closest I can come to replecating the structure of a Temp Power Sweep
	# dt = np.dtype([('O',[('Temp','O',(1,1))],(1,2))])

	
	i=0
	self.Sweep_Array = np.zeros(metadata.Heater_Voltage.shape[0]*metadata.Powers.shape[0]*metadata.Freq_Range.shape[0], dtype = self.sweep_data_columns)
	for freq_range_num in xrange(metadata.Freq_Range.shape[0]):
			if metadata.Heater_Voltage.shape[0] == 1:
				heater_voltages = metadata.Freq_Range # non temp sweep, single freq_range, powersweep
			else:					
				heater_voltages = self._extract_type(metadata.Freq_Range[freq_range_num,1]['Temp'])
			#start here for single res powersweep
			for heater_voltage_num in xrange(heater_voltages.shape[0]):
				sweep_powers = self._extract_type(heater_voltages[heater_voltage_num,1], field = 'PowerSweep')
				for sweep in sweep_powers[:,0:sweep_powers.shape[1]]:
					self._define_sweep_array(i, Fstart = self.metadata.Range[freq_range_num,0],
												Fstop = self.metadata.Range[freq_range_num,1],
												Heater_Voltage = self.metadata.Heater_Voltage[heater_voltage_num],
												Thermometer_Voltage_Bias = self.metadata.Thermometer_Voltage_Bias[heater_voltage_num] if hasattr(metadata,'Thermometer_Voltage_Bias') else 0,#set to zero unless there is an array of temps in the ScanData
												Pinput_dB = sweep[0].squeeze()[()] - metadata.Atten_Added_At_NA if metadata.Atten_Added_At_NA != None else sweep[0].squeeze()[()], #we only want the power coming out of the source, i.e. the NA
												S21 = sweep[1].squeeze()[()],
												Frequencies = sweep[2].squeeze()[()],
												Temperature_Readings = sweep[3].squeeze()[()] if (sweep.size > 3) and (np.shape(sweep[3].squeeze()[()])[0] != 0) else np.array([0]), #set to zero unless there is an array of temps in the ScanData
												Is_Valid = True)
					i = i + 1

	if  hasattr(metadata,'Thermometer_Voltage_Bias'):
		del(metadata.__dict__['Thermometer_Voltage_Bias'])
	del(metadata.__dict__['Powers'])
	del(metadata.__dict__['Heater_Voltage'])
	del(metadata.__dict__['Range'])
	del(metadata.__dict__['Freq_Range'])

	# if self.metadata.Atten_NA_Output == None: #redundant to have both
	# 	del(self.metadata.__dict__['Atten_NA_Output'])
	# else:
	# 	del(self.metadata.__dict__['Atten_Added_At_NA'])
