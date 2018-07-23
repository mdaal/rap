import tempfile
import scipy as sp
import numpy as np


def _extract_type(obj, return_type=None, field=None):
    '''scandata object, obj, has a lot of single element arrays of arrays. this function
       gets the element. e.g scandata may have [[[ele]]] instead of callling ele =
       scandata[0][0][0], use this function to get ele. If ele is another structured
       numpy array with field name 'myfield', using keyword field = 'myfield' will get
       the data at field. The function will cast ele to be in the data type return_type.
       e.g. return_type = 'str' returns a string. If return_type is None, ele is returned
       as whatever type it was saved as in [[[ele]]]'''

    def cast(_obj):
        if (return_type is not None) & (_obj is not None):
            _obj = return_type(_obj)
        return _obj

    def itemloop(_obj):
        while True:
            try:
                _obj = _obj.item()
            except:
                # remove any whitespace for casting complex numbers
                if type(_obj) is str:
                    _obj = _obj.replace(" ", "").replace("i", "j")
                return cast(_obj)
        return cast(_obj)

    if field is None:
        obj = itemloop(obj)

    else:
        while obj.dtype == np.dtype('O'):
            obj = obj.item()

        if isinstance(obj.item(), str):
            obj = None
            print('Expected dictionary containing field named {0} is not found. '
                  'Returning None'.format(field))
        # if the object does not simply contain a string, e.g  [u'InP #2'], do this
        else:
            try:
                obj = obj[field]
            except:
                obj = None
                print('Field named {0} is not found. Cannot extract. '
                      'Returning None'.format(field))
        obj = itemloop(obj)
    return obj


def _define_pulse_data_columns_legacy_gui(n_points, n_triggers, n_points_noise,
                                          n_triggers_noise, n_energies):
    ''' Create the pulse_data_columns_list which is used to define the dtype of the
        Pulse_Array'''
    pulse_data_columns_list = [
        ("trace_I", np.float64, (n_triggers, n_points)),  # I traces
        ("trace_Q", np.float64, (n_triggers, n_points)),  # Q traces
        ("noise_I", np.float64, (n_triggers_noise, n_points_noise)),  # I noise
        ("noise_Q", np.float64, (n_triggers_noise, n_points_noise)),  # Q noise
        ("noise_center", np.complex128),  # noise S21 location I + jQ
        ("wavelengths", np.float64, (n_energies, )),  # wavelengths nm
        ("E", np.float64, (n_energies, )),  # energies in eV
        ("F0", np.float64),  # readout frequency
        ("sample_rate", np.float64),  # trace sample rate
        ("atten", np.float64)]  # gui attenuation

    pulse_data_columns = np.dtype(pulse_data_columns_list)

    return pulse_data_columns_list, pulse_data_columns


def _define_pulse_data_columns_JPL(n_points, n_triggers, n_points_noise,
                                   n_triggers_noise, n_energies):
    ''' Create the pulse_data_columns_list which is used to define the dtype of the
        Pulse_Array'''
    pulse_data_columns_list = [
        ("trace_I", np.float64, (n_triggers, n_points)),  # I traces
        ("trace_Q", np.float64, (n_triggers, n_points)),  # Q traces
        ("noise_I", np.float64, (n_triggers_noise, n_points_noise)),  # I noise
        ("noise_Q", np.float64, (n_triggers_noise, n_points_noise)),  # Q noise
        ("noise_center", np.complex128),  # noise S21 location I + jQ
        ("wavelengths", np.float64, (n_energies, )),  # wavelengths nm
        ("E", np.float64, (n_energies, )),  # energies in eV
        ("F0", np.float64),  # readout frequency
        ("sample_rate", np.float64),  # trace sample rate
        ("atten", np.float64)]  # programmable attenuator attenuation for pulse data

    pulse_data_columns = np.dtype(pulse_data_columns_list)

    return pulse_data_columns_list, pulse_data_columns


def _define_pulse_array(Pulse_Array, index, **field_names):
    for field_name in field_names:
        Pulse_Array[field_name][index] = field_names[field_name]

    return Pulse_Array


def _unpack_data_structure(tags, receiving_dict, struct_to_be_unpacked):
    ''' takes data specfied in tag from struct_to_be_unpacked and adds it to
        receiving_dict in the format specificed in tags the  placement of the data within
        receiving_dict is also specificed by tags, which is a list od tuples
        [(fieldname_in_struct_to_be_unpacked, destination_key_in_receiving_dict,
        format_in_receiving_dict, [optional  subfieldname_in_struct_to_be_unpacked]), ...]
    '''
    for t in tags:
        # try:
            if t[1].find(':') > -1:  # The case of a dictionary
                t1 = t[1].split(':')

                # This try/except block is for the case where metadata.__dict__['?']
                # is a dictionary
                try:
                    receiving_dict[t1[0]].update(
                        [(t1[1], _extract_type(struct_to_be_unpacked[t[0]],
                                               return_type=t[2],
                                               field=t[3] if len(t) > 3 else None))])
                except:
                    receiving_dict[t1[0]] = dict(
                        [(t1[1], _extract_type(struct_to_be_unpacked[t[0]],
                                               return_type=t[2],
                                               field=t[3] if len(t) > 3 else None))])
            else:
                receiving_dict[t[1]] = _extract_type(struct_to_be_unpacked[t[0]],
                                                     return_type=t[2],
                                                     field=t[3] if len(t) > 3 else None)
