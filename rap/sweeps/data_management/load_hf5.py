def load_hf5(self, tablepath, filename = database_location):
	''' table path is path to the database to be loaded starting from root. e.g. self.load_hf5('/Run44b/T201312102229')
	filename is the name of the hf5 database to be accessed for the  table informaiton'''

	if not os.path.isfile(filename):
		print('Speficied h5 database does not exist. Aborting...')
		return 
	
	wmode = 'a'
	
	#delete previous metadata object
	del(self.metadata)
	self.metadata = metadata()
	del(self.loop)
	self.loop = loop()

	# use "with" context manage to ensure file is always closed. no need for fileh.close()
	with tables.open_file(filename, mode = wmode) as fileh:
		table = fileh.get_node(tablepath)	
		self.Sweep_Array = table.read()
		for data in self.metadata.__dict__.keys():
			try:
				exec('self.metadata.{0} = table.attrs.{0}'.format(data))
			except:
				print('Table metadata is missing {0}. Setting to None'.format(data))
				exec('self.metadata.{0} = None'.format(data))
	self.sweep_data_columns = self.Sweep_Array.dtype