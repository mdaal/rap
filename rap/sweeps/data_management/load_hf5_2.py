from utils import  _define_sweep_data_columns

import tables
import os
import logging

def load_hf5_2(metadata, hf5_database_path, tablepath): 
	''' This function is for loading data taken with KIDs_DAQ_75uW. It use the columns defined in that hf5 file to 
	define the columns  in  sweep_data_columns 
	table path is path to the database to be loaded starting from root. e.g. load_hf5('/Run44b/T201312102229')
	hf5_database_path is the name of the hf5 database to be accessed for the  table informaiton'''
	
	if not os.path.isfile(hf5_database_path):
		logging.error('Speficied h5 database does not exist. Aborting...')
		return 
	
	wmode = 'a'



	# use "with" context manage to ensure file is always closed. no need for fileh.close()
	with tables.open_file(hf5_database_path, mode = wmode) as fileh:
		table = fileh.get_node(tablepath)	
		Sweep_Array = table.read()
		for key in table.attrs.keys:
			exec('metadata.{0} = table.attrs.{0}'.format(key))

	imported_sweep_data_columns = Sweep_Array.dtype
	fsteps = imported_sweep_data_columns['Frequencies'].shape[0]
	tpoints =  imported_sweep_data_columns['Temperature_Readings'].shape[0]
	sweep_data_columns_list, sweep_data_columns = _define_sweep_data_columns(fsteps, tpoints)
	

	for name in imported_sweep_data_columns.names:
		if name not in sweep_data_columns.names:
			sweep_data_columns_list.append((name,imported_sweep_data_columns[name] ))
	sweep_data_columns = np.dtype(sweep_data_columns_list)
	Sweep_Array = np.array(Sweep_Array, dtype = sweep_data_columns)

	return Sweep_Array, sweep_data_columns, sweep_data_columns_list
