import datetime ##--> use arrow
import tables
import os

def save_hf5(metadata, Sweep_Array, sweep_data_columns, hf5_database_path,overwrite = False):
	'''Saves current Sweep_Array into table contained in the hdf5 file speficied by hf5_database_path.
	If overwite = True, Sweep_Array will overwright whatever is previous table data there is.
	'''
	
	if not os.path.isfile(hf5_database_path):
		print('Speficied h5 database does not exist. Creating new one.')
		pos = hf5_database_path.rfind('/')
		if pos >= 0:
			try:
				os.makedirs(hf5_database_path[0:pos+1])
			except OSError:
				print('{0} exists...'.format(hf5_database_path[0:pos+1]))
		wmode = 'w'
	else:
		print('Speficied h5 database exists and will be updated.')
		wmode = 'a'
		
	db_title = 'Aggregation of Selected Data Sets'
	group_name = 'Run' + metadata.Run
	group_title = metadata.Test_Location


	try:  # for forward compatabiliity with 75uW python DAQ
		d = datetime.datetime.strptime(metadata.Measurement_Start_Time , '%Y%m%d%H%M')
	except:
		pass

	try:
		# case for scan data date
		d = datetime.datetime.strptime(metadata.Time_Created, '%B %d, %Y  %I:%M:%S.%f %p') # slightly wrong %f is microseconds. whereas we want milliseconds.
	except:
		pass
	try:
		#Case for sonnet date
		d = datetime.datetime.strptime(metadata.Time_Created, '%m/%d/%Y %H:%M:%S')
	except:
		pass
	sweep_data_table_name = 'T' + d.strftime('%Y%m%d%H%M')

	

	with tables.open_file(hf5_database_path, mode = wmode, title = db_title ) as fileh:
		try:
			table_path = '/' + group_name + '/' + sweep_data_table_name
			sweep_data_table = fileh.get_node(table_path)

			if overwrite == True:
				print('Table {0} exists. Overwriting...'.format(table_path))
				sweep_data_table.remove()
				sweep_data_table = fileh.create_table('/'+ group_name,sweep_data_table_name,description=sweep_data_columns,title = 'Sweep Data Table',filters=tables.Filters(0), createparents=True)
			else:
				print('Table {0} exists. Aborting...'.format(table_path))
				return
		except:
			print('Creating table {0}'.format('/'+ group_name+'/'+sweep_data_table_name))
			sweep_data_table = fileh.create_table('/'+ group_name,sweep_data_table_name,description=sweep_data_columns,title = 'Sweep Data Table',filters=tables.Filters(0), createparents=True)
		
		# copy Sweep_Array to sweep_data_table
		sweep_data_table.append(Sweep_Array)

		# Save metadata
		for data in metadata.__dict__.keys():
			exec('sweep_data_table.attrs.{0} = metadata.{0}'.format(data))
			if metadata.__dict__[data] == None:
				print('table metadata {0} not defined and is set to None'.format(data))	

		sweep_data_table.flush()	
