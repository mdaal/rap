def load_hf5_2(self, database_filename, tablepath):
	''' This function is for loading data taken with KIDs_DAQ_75uW. It use the columns defined in that hf5 file to 
	define the columns  in  self.sweep_data_columns 
	table path is path to the database to be loaded starting from root. e.g. self.load_hf5('/Run44b/T201312102229')
	database_filename is the name of the hf5 database to be accessed for the  table informaiton'''
	
	if not os.path.isfile(database_filename):
		logging.error('Speficied h5 database does not exist. Aborting...')
		return 
	
	wmode = 'a'

	#delete previous metadata object
	del(self.metadata)
	self.metadata = metadata()
	del(self.loop)
	self.loop = loop()

	# use "with" context manage to ensure file is always closed. no need for fileh.close()
	with tables.open_file(database_filename, mode = wmode) as fileh:
		table = fileh.get_node(tablepath)	
		self.Sweep_Array = table.read()
		for key in table.attrs.keys:
			#exec('self.measurement_metadata["{0}"] = table.attrs.{0}'.format(key))
			exec('self.metadata.{0} = table.attrs.{0}'.format(key))

	#self.sweep_data_columns = self.Sweep_Array.dtype
	imported_sweep_data_columns = self.Sweep_Array.dtype
	try:
		self.metadata.Cable_Calibration = self._Cable_Calibration
		print('Cable Calibraiton data found and saved in Sweep_Array metadata.')
	except:
		pass

	fsteps = imported_sweep_data_columns['Frequencies'].shape[0]
	tpoints =  imported_sweep_data_columns['Temperature_Readings'].shape[0]
	self._define_sweep_data_columns(fsteps, tpoints)
	#self.sweep_data_columns_list

	for name in imported_sweep_data_columns.names:

		if name not in self.sweep_data_columns.names:
			self.sweep_data_columns_list.append((name,imported_sweep_data_columns[name] ))
	self.sweep_data_columns = np.dtype(self.sweep_data_columns_list)
	self.Sweep_Array = np.array(self.Sweep_Array, dtype = self.sweep_data_columns)
