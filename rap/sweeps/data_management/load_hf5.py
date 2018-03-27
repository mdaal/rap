import tables
import os
import logging

def load_hf5(metadata, hf5_database_path, tablepath):
    ''' table path is path to the database to be loaded starting from root. e.g. load_hf5('/Run44b/T201312102229')
    hf5_database_path is the name of the hf5 database to be accessed for the  table informaiton'''

    if not os.path.isfile(hf5_database_path):
        logging.error('Speficied h5 database does not exist. Aborting...')
        return

    wmode = 'a'


    # use "with" context manage to ensure file is always closed. no need for fileh.close()
    with tables.open_file(hf5_database_path, mode = wmode) as fileh:
        table = fileh.get_node(tablepath)
        Sweep_Array = table.read()
        for data in metadata.__dict__.keys():
            try:
                exec('metadata.{0} = table.attrs.{0}'.format(data))
            except:
                print('Table metadata is missing {0}. Setting to None'.format(data))
                exec('metadata.{0} = None'.format(data))
    sweep_data_columns = Sweep_Array.dtype

    return Sweep_Array, sweep_data_columns
