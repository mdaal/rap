def construct_hf5_toc(self,filename = database_location):
	''' Creates a table of contents (toc) of the hf5 database storing all the sweep_data.
	very useful for finding the name and locaiton of a table in the database'''
	if not os.path.isfile(filename):
		print('Speficied h5 database does not exist. Aborting...')
		return 
	
	wmode = 'a'

	# use "with" context manage to ensure file is always closed. no need for fileh.close()
	with tables.open_file(filename, mode = wmode) as fileh:
		table_list = [g for g in fileh.walk_nodes(classname = 'Table')]
		num_tables = len(table_list)
		TOC = np.zeros(num_tables, dtype = self.data_set_contents)
		index = 0
		for table in table_list:
			TOC['Run'][index] 						= table.get_attr('Run') 
			TOC['Time_Created'][index] 				= table.get_attr('Time_Created')
			TOC['Num_Temperature_Readings'][index]	= table.get_attr('Num_Temperatures') if table.get_attr('Num_Temperatures') !=None else 0
			#TOC['Num_Ranges'][index] 				= table.get_attr('Num_Ranges') if 'Num_Ranges' in table.attrs._v_attrnames else 1
			TOC['Num_Ranges'][index] 				= table.get_attr('Num_Ranges') if table.get_attr('Num_Ranges') !=None else 0
			TOC['Num_Powers'][index] 				= table.get_attr('Num_Powers') if table.get_attr('Num_Powers') !=None else 0
			TOC['Num_Temperatures'][index] 			= table.get_attr('Num_Heater_Voltages') if table.get_attr('Num_Heater_Voltages') !=None else 0
			TOC['Sensor'][index] 					= table.get_attr('Sensor') if table.get_attr('Sensor') !=None else ''
			TOC['Ground_Plane'][index] 				= table.get_attr('Ground_Plane') if table.get_attr('Ground_Plane') !=None  else ''
			TOC['Path'][index] 						= table._v_pathname
			index += 1 

	self.TOC = TOC
	print(TOC)
