import tables
import numpy as np
import os

def construct_hf5_toc(hf5_database_path):
    ''' Creates a table of contents (toc) of the hf5 database storing all the sweep_data.
    very useful for finding the name and locaiton of a table in the database'''

    if not os.path.isfile(hf5_database_path):
        print('Speficied h5 database does not exist. Aborting...')
        return

    wmode = 'a'

    # need to get rid of data set contents... All calls of it are now in construct_hf5_toc
    data_set_contents = np.dtype([
        ("Run"                            , 'S10'),
        ("Time_Created"                    , 'S40'), # 'S40' for format December 23, 2012 12:34:65.675 PM;  'S12' for format '%Y%m%d%H%M'
        ("Num_Ranges"                    , np.uint8), # uint8 is Unsigned integer (0 to 255)
        ("Num_Powers"                    , np.uint8),
        ("Num_Temperature_Readings"        , np.uint8), # the number of temperature reading taking for each sweep
        ("Num_Temperature_Set_Points"    , np.uint8), # the number of frsdge temperatures the used in sweep
        ("Sensor"                        , 'S20'),
        ("Ground_Plane"                    , 'S30'),
        ("Path"                            , 'S100'),
        ])


    # use "with" context manage to ensure file is always closed. no need for fileh.close()
    with tables.open_file(hf5_database_path, mode = wmode) as fileh:
        table_list = [g for g in fileh.walk_nodes(classname = 'Table')]
        num_tables = len(table_list)
        TOC = np.zeros(num_tables, dtype = data_set_contents)
        index = 0
        for table in table_list:
            TOC['Run'][index]                                 = table.get_attr('Run')
            TOC['Time_Created'][index]                         = table.get_attr('Time_Created')
            TOC['Num_Temperature_Readings'][index]            = table.get_attr('Num_Temperatures') if table.get_attr('Num_Temperatures') !=None else 0
            #TOC['Num_Ranges'][index]                         = table.get_attr('Num_Ranges') if 'Num_Ranges' in table.attrs._v_attrnames else 1
            TOC['Num_Ranges'][index]                         = table.get_attr('Num_Ranges') if table.get_attr('Num_Ranges') !=None else 0
            TOC['Num_Powers'][index]                         = table.get_attr('Num_Powers') if table.get_attr('Num_Powers') !=None else 0
            TOC['Num_Temperature_Set_Points'][index]        = table.get_attr('Num_Heater_Voltages') if table.get_attr('Num_Heater_Voltages') !=None else 0
            TOC['Sensor'][index]                             = table.get_attr('Sensor') if table.get_attr('Sensor') !=None else ''
            TOC['Ground_Plane'][index]                         = table.get_attr('Ground_Plane') if table.get_attr('Ground_Plane') !=None  else ''
            TOC['Path'][index]                                 = table._v_pathname
            index += 1
    print('Column names are:')
    print(TOC.dtype.names)
    print('Rows, i.e. data, are:')
    print(TOC)
