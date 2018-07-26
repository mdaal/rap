import numpy as np
from .utils import (_define_sweep_array, _define_sweep_data_columns_JPL)
from ...calibration.iq_mixer import mixer_cal


def load_JPL_sweep_data(metadata, sweep_file, iq_mixer=None):
    # load the data from the sweep file
    npzfile = np.load(sweep_file)
    f_loop = npzfile['freqs'] * 1e9
    z_cal = npzfile['z0']
    if iq_mixer is not None:
        z_loop = mixer_cal(npzfile['z'] - z_cal, iq_mixer)
    else:
        z_loop = npzfile['z'] - z_cal
    meta = npzfile['meta'].item()

    # add the summary string to the metadata object
    metadata.summary_string = meta['summary']
    metadata.Loop_Data_Column_Names = ("frequencies", "s21")

    # create the sweep array
    sweep_data_columns_list, sweep_data_columns = _define_sweep_data_columns_JPL(meta)
    Sweep_Array = np.zeros(1, dtype=sweep_data_columns)
    _define_sweep_array(Sweep_Array, 0,
                        attenuation=meta['attenuation'],
                        power=meta['power'],
                        frequencies=f_loop,
                        s21=z_loop,
                        s21_baseline=z_cal)

    return sweep_data_columns_list, sweep_data_columns, Sweep_Array
