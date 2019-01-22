from .utils import  _define_sweep_data_columns, _define_sweep_array
from ..sweep_array.pick_loop import pick_loop

import tempfile
import io
import numpy as np



def load_touchstone(metadata, filename):
    ''' The function loads S21 and Freq from  Sonnet .s2p or .s3p files into the Sweep_Array structured np array
    All Sij are extracted, but only  S21 is saved into Sweep_Array. Future editions of this code might  find need
    to load other Sij becuase S21.

    The function only loads one transmission array (S21).  pick_loop = True immediatly selectes this loop as the
    current loop.
    '''


    dt_s2p = [('Freq', np.float64), ('S11a', np.float64), ('S11b', np.float64), ('S12a', np.float64), ('S12b', np.float64),
                                    ('S21a', np.float64), ('S21b', np.float64), ('S22a', np.float64), ('S22b', np.float64)]

    dt_s3p = [('Freq', np.float64), ('S11a', np.float64), ('S11b', np.float64), ('S12a', np.float64), ('S12b', np.float64), ('S13a', np.float64), ('S13b', np.float64),
                                    ('S21a', np.float64), ('S21b', np.float64), ('S22a', np.float64), ('S22b', np.float64), ('S23a', np.float64), ('S23b', np.float64),
                                    ('S31a', np.float64), ('S31b', np.float64), ('S32a', np.float64), ('S32b', np.float64), ('S33a', np.float64), ('S33b', np.float64)]

    j = np.complex(0,1)
    data_format_converter_dict = {  b'RI': lambda Sij_r, Sij_i: Sij_r + j*Sij_i,
                                    # BUG: 'MA': lambda Sij_m, Sij_a: Sij_m * np.exp(j*np.pi/180 * Sij_a),
                                    b'DB': lambda Sij_mdb, Sij_a: np.power(10,Sij_mdb/20.0) * np.exp(j*np.pi/180 * Sij_a)
                                     }


    # if type(b'bytes') is str:
    #     kwargs = {}
    # else:
    #     kwargs = {'mode': 'w', 'encoding': 'utf-8'}
    kwargs = {}
    with tempfile.TemporaryFile(**kwargs) as tmp:
        with io.open(filename, mode='rb') as f:
            # The following while loop copies the .sNp file into a temp file, which is destroyed when closed,
            # such that the tmp file is formated in the way np.loadtxt can read the data.
            indented = False
            prev_line = b''
            m = 1. # for frequency base conversion
            data_format = b'' # The sparameter data format. Might be RI, MA, DB
            while 1:
                line  = f.readline().replace(b'\n',b'')

                pos = f.tell()
                if line == b'': # End of file is reached
                    break
                elif line.startswith(b'! Data File Written:'): # Save as Metadata
                    metadata.Time_Created = line.split(b'! Data File Written:')[1].strip().decode("utf-8", "backslashreplace") 
                    # Note: line.split(b'! Data File Written:')[1].strip() is a bytes object, e.g. b'...',=.  
                    # We want a unicode string, so we decode the bytes object using the encoding "utf-8". If a character that is not 
                    # decodable in "utf-8" is encountered, the "backslashreplace" option inserts a \xNN escape sequence. 
                    # alternatively we cound use the "ignore" to just leave the character out of the Unicode result.
                    tmp.write(line + b'\n')
                elif line.startswith(b'! From Project:') | line.startswith(b'! From Emgraph Data:'): # Save as Metadata
                    metadata.Run = line.split(b':')[1].strip().decode("utf-8", "backslashreplace")
                    #self.metadata.IS_Sonnet_Simulation = True
                    tmp.write(line + b'\n')
                elif line[0] == '#' or (type(line) is not str and line[0] == 35):
                    line  = line.replace(b'#',b'!#')
                    split_line = line.split(b' ')
                    if b'GHZ' in split_line:
                        m = 1.0e9
                    if b'S' not in split_line: # Check that file contatins S-parameters
                        raise ValueError('Not an S-parameter file.')
                    S_index = split_line.index(b'S')    # he sparameter data format is declared next in split_line list
                    if split_line[S_index+1] not in data_format_converter_dict.keys():
                        raise ValueError('The S-parameter file is "{}" format. Should be "RI" (preferred), "MA", or "DB" format.'.format(split_line[S_index+1]))
                    else:
                        data_format = split_line[S_index+1]


                    freq_convert = lambda s: s*m #Convert to Hertz
                    tmp.write(line + b'\n')

                elif line[0] == b' ': # in S matrix definition block
                    prev_line = prev_line + b' ' + line.strip() + b' '
                    next_line = f.readline()
                    # if next line is NOT indented date, then S matrix definition block is finished
                    # and we write it to tmp on a single line.
                    # for .s2p files the S matrix is fully defined on one line of f
                    # for .s3p files, the S matrix is defined in three lines. second two are indented.

                    # if not ((next_line[0] == '') | (next_line[0] == ' ')): # Changing this line to be consistent with line below...
                    if not ((next_line == b'') or (next_line[0] == b' ')):
                        tmp.write(prev_line)
                        tmp.write(b'\n')
                        prev_line = b''
                    f.seek(pos,0)

                elif line[0] == b'!':
                    tmp.write(line + b'\n')

                else:
                    tmp.write(line)
                    next_line = f.readline()
                    # add \n to line if it does not begin a S matrix definition block
                    # if not ((next_line[0] == '') | (next_line[0] == ' ')): # Changed on 7/11/17 bc Nick was have problems reading in .s2p files
                    if not ((next_line == b'') or (next_line[0] == b' ')):
                        tmp.write(b'\n')
                    f.seek(pos,0)

        tmp.seek(0,0)
        if filename.endswith('.s2p'):
            dt = dt_s2p
        elif filename.endswith('.s3p'):
            dt = dt_s3p
        Touchstone_Data = np.loadtxt(tmp, dtype=dt, comments=['!', '#'], delimiter=None, converters=None, skiprows=0, usecols=None, unpack=False, ndmin=0)

    tpoints = 0
    sweep_data_columns_list, sweep_data_columns = _define_sweep_data_columns(metadata,Touchstone_Data.size, tpoints)


    Sweep_Array = np.zeros(1, dtype = sweep_data_columns)

    _define_sweep_array(Sweep_Array, 0,
                                Fstart = freq_convert(Touchstone_Data['Freq'].min()), #Hz
                                Fstop = freq_convert(Touchstone_Data['Freq'].max()), #Hz
                                S21 = data_format_converter_dict[data_format](Touchstone_Data['S21a'],Touchstone_Data['S21b']),
                                Frequencies = freq_convert(Touchstone_Data['Freq']), #Hz
                                #Pinput_dB = 0,
                                Is_Valid = True,
                                #Mask = False, needs to be an array of lengh of S21
                                Chi_Squared = 0,
                                )


    metadata.Data_Source = filename
    #self.metadata.Min_Freq_Resolution = np.abs(Touchstone_Data['Freq'][:-1]-Touchstone_Data['Freq'][1:]).min()
    metadata.Min_Freq_Resolution = np.abs(Touchstone_Data['Freq'][0] - Touchstone_Data['Freq'][-1])/Touchstone_Data['Freq'].size #use average freq resolution
    metadata.Num_Temperatures = 0.0
    metadata.Num_Ranges = 1.0
    metadata.Num_Powers = 1.0
    metadata.Num_Heater_Voltages = 1.0
    metadata.Loop_Data_Column_Names = ("Frequencies","S21")

    return sweep_data_columns_list, sweep_data_columns, Sweep_Array
