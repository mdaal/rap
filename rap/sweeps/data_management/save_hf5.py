def save_hf5(self, filename = database_location, overwrite = False):
	'''Saves current self.Sweep_Array into table contained in the hdf5 file speficied by filename.
	If overwite = True, self.Sweep_Array will overwright whatever is previous table data there is.
	'''
	
	if not os.path.isfile(filename):
		print('Speficied h5 database does not exist. Creating new one.')
		pos = filename.find('/')
		if pos >= 0:
			try:
				os.makedirs(filename[0:pos+1])
			except OSError:
				print('{0} exists...'.format(filename[0:pos+1]))
		wmode = 'w'
	else:
		print('Speficied h5 database exists and will be updated.')
		wmode = 'a'
		
	db_title = 'Aggregation of Selected Data Sets'
	group_name = 'Run' + self.metadata.Run
	group_title = self.metadata.Test_Location


	try:  # for forward compatabiliity with 75uW python DAQ
		d = datetime.datetime.strptime(self.metadata.Measurement_Start_Time , '%Y%m%d%H%M')
	except:
		pass

	try:
		# case for scan data date
		d = datetime.datetime.strptime(self.metadata.Time_Created, '%B %d, %Y  %I:%M:%S.%f %p') # slightly wrong %f is microseconds. whereas we want milliseconds.
	except:
		pass
	try:
		#Case for sonnet date
		d = datetime.datetime.strptime(self.metadata.Time_Created, '%m/%d/%Y %H:%M:%S')
	except:
		pass
	sweep_data_table_name = 'T' + d.strftime('%Y%m%d%H%M')

	

	with tables.open_file(filename, mode = wmode, title = db_title ) as fileh:
		try:
			table_path = '/' + group_name + '/' + sweep_data_table_name
			sweep_data_table = fileh.get_node(table_path)

			if overwrite == True:
				print('Table {0} exists. Overwriting...'.format(table_path))
				sweep_data_table.remove()
				sweep_data_table = fileh.create_table('/'+ group_name,sweep_data_table_name,description=self.sweep_data_columns,title = 'Sweep Data Table',filters=tables.Filters(0), createparents=True)
			else:
				print('Table {0} exists. Aborting...'.format(table_path))
				return
		except:
			print('Creating table {0}'.format('/'+ group_name+'/'+sweep_data_table_name))
			sweep_data_table = fileh.create_table('/'+ group_name,sweep_data_table_name,description=self.sweep_data_columns,title = 'Sweep Data Table',filters=tables.Filters(0), createparents=True)
		
		# copy Sweep_Array to sweep_data_table
		sweep_data_table.append(self.Sweep_Array)

		# Save metadata
		for data in self.metadata.__dict__.keys():
			exec('sweep_data_table.attrs.{0} = self.metadata.{0}'.format(data))
			if self.metadata.__dict__[data] == None:
				print('table metadata {0} not defined and is set to None'.format(data))	

		sweep_data_table.flush()	

		# try:
		# 	TOC = fileh.get_node('/Contents') # is a table
		# except:
		# 	print('Creating h5 data set table of contents')
		# 	TOC = fileh.create_table('/', 'Contents', self.data_set_contents, "Table listing all tables contained in h5 file", tables.Filters(0)) #tables.Filters(0) means there is no data compression

		# TOC.append()

	# title = 'Data from Run ' + self.metadata.Run + ', Sensor: ' + self.metadata.Sensor + ', Ground Plane: ' + self.metadata.Ground_Plane

	# #determine type of  measurement...
	# if  (self.Sweep_Array.size == 1) | (np.abs(self.Sweep_Array['Fstop'] - self.Sweep_Array['Fstart']).max() >= 100e6):
	# 	groupname = 'Survey'
	# elif (np.unique(self.Sweep_Array['Heater_Voltage']).size > 1) && (np.unique(self.Sweep_Array['Pinput_dB']).size == 1):
	# 	groupname = 'T_Sweep'
	# elif (np.unique(self.Sweep_Array['Heater_Voltage']).size == 1) && (np.unique(self.Sweep_Array['Pinput_dB']).size > 1):
	# 	groupname = 'P_Sweep'
	# elif (np.unique(self.Sweep_Array['Heater_Voltage']).size > 1) && (np.unique(self.Sweep_Array['Pinput_dB']).size > 1):
	# 	groupname = 'TP_Sweep'
	# else:
	# 	groupname = 'Sweep'

	# 	groupname = 'T' + str(np.unique(self.Sweep_Array['Heater_Voltage']).size) + 'P' +  str(np.unique(self.Sweep_Array['Pinput_dB']).size)	
