import os
import warnings
import numpy as np
from scipy.io import loadmat
from scipy.signal import welch, csd, resample
from datetime import datetime
from astropy.constants import h, c
from matplotlib import pyplot as plt

from .utils import (_unpack_data_structure, _define_pulse_array,
                    _define_pulse_data_columns_legacy_gui)


def load_legacy_pulse_gui_data(metadata, gui_data_path, energies, sweep,
                               use_pulse_gui_psd=True):

    # find latest pulse_config_xxxxx.mat file in gui_data_path
    file_list = os.listdir(gui_data_path)
    lastest_pulse_config_modification_time = 0
    flag = 0
    for file_ in file_list:
        if file_.startswith('pulse_config'):
            flag += 1
            parts = file_.replace('.', '_').split('_')
            if len(parts) == 3:
                config_file = os.path.join(gui_data_path, file_)
                break
            elif len(parts) == 4:
                lastest_pulse_config_modification_time = int(parts[2])
                if int(parts[2]) >= lastest_pulse_config_modification_time:
                    lastest_pulse_config_modification_time = int(parts[2])
                    config_file = os.path.join(gui_data_path, file_)
    if flag == 0:
        raise IOError('No pulse_config.mat files found in the gui_data_path')

    # load in the config file
    config = loadmat(config_file)

    # initialize data dictionaries
    data_dict = dict()
    data_dict['curr_config'] = dict()

    # define the configuration file format
    config_tags = [('__header__', 'curr_config:Time_Created', str),
                   ('curr_config', 'curr_config:adtime', float, 'adtime'),
                   ('curr_config', 'curr_config:tottime', float, 'tottime'),
                   ('curr_config', 'curr_config:atten1', float, 'atten1'),
                   ('curr_config', 'curr_config:atten2', float, 'atten2'),
                   ('curr_config', 'curr_config:samprate', float, 'samprate'),
                   ('curr_config', 'curr_config:do_plots', int, 'do_plots'),
                   ('curr_config', 'curr_config:f01', float, 'f01'),
                   ('curr_config', 'curr_config:f02', float, 'f02'),
                   ('curr_config', 'curr_config:nsig1', float, 'nsig1'),
                   ('curr_config', 'curr_config:nsig2', float, 'nsig2'),
                   ('curr_config', 'curr_config:data_path', str, 'data_path'),
                   ('curr_config', 'curr_config:numpts', int, 'numpts'),
                   ('curr_config', 'curr_config:config_file', str, 'config_file'),
                   ('curr_config', 'curr_config:sweep_data', str, 'sweep_data'),
                   ('curr_config', 'curr_config:res_nums', list, 'res_nums'),
                   ('curr_config', 'curr_config:fit_res', str, 'fit_res'),
                   ('curr_config', 'curr_config:cent1', complex, 'cent1'),
                   ('curr_config', 'curr_config:cent2', complex, 'cent2'),
                   ('curr_config', 'curr_config:pulse_path', str, 'pulse_path'),
                   ('curr_config', 'curr_config:der_flag', int, 'der_flag'),
                   ('curr_config', 'curr_config:of_flag', int, 'of_flag'),
                   ('curr_config', 'curr_config:plot_phase', int, 'plot_phase'),
                   ('curr_config', 'curr_config:take_noise_box', int, 'take_noise_box'),
                   ('curr_config', 'curr_config:noise_adtime', float, 'noise_adtime'),
                   ('curr_config', 'curr_config:noise_tottime', float, 'noise_tottime')]

    # put data from matlab config file into a dictionary with correct data types
    _unpack_data_structure(config_tags, data_dict, config)

    # extract the creation time and date in format : 'Wed Aug 09 17:15:14 2017'
    time = data_dict['curr_config']['Time_Created']
    index = time.find('Created on: ') + len('Created on: ')
    data_dict['curr_config']['Time_Created'] = time[index:]
    # remove trailing quotation mark if needed
    if data_dict['curr_config']['Time_Created'][-1] == "'":
        time = data_dict['curr_config']['Time_Created'][:-1]
        data_dict['curr_config']['Time_Created'] = time
    dt_start = datetime.strptime(data_dict['curr_config']['Time_Created'],
                                 '%a %b %d %H:%M:%S %Y')

    # find and load appropriate sweep files
    sweep.load_legacy_sweep_gui_data(gui_data_path)
    assert len(sweep.Sweep_Array) == 2, "Too much sweep data in the gui data path. " \
        "Can't determine which cooresponds to the pulse data"

    # find the pulse_data.dat file that matches the pulse_config.mat file
    if lastest_pulse_config_modification_time != 0:
        file_name = 'pulse_data_' + str(lastest_pulse_config_modification_time)
    else:
        file_name = 'pulse_data'
    if (os.path.isfile(os.path.join(gui_data_path, file_name + '.dat')) and
       os.path.isfile(os.path.join(gui_data_path, file_name + '.ns'))):
        data_file = os.path.join(gui_data_path, file_name + '.dat')
        if use_pulse_gui_psd:
            noise_file = os.path.join(gui_data_path, file_name + '.ns')
            # record the noise sample rate
            sample_rate = data_dict['curr_config']['samprate']
            # record the number of points per trace
            n_points_noise = data_dict['curr_config']['noise_adtime'] * sample_rate
            decfac = 1
        else:
            temp = str(int(sweep.Sweep_Array["Temperature"][0] / 0.001))
            group = str(sweep.Sweep_Array["Resonator_Group"][0, 0]) + 'a'
            atten = str(sweep.Sweep_Array['Noise_Chan_Input_Atn'][0])
            file_ = temp + '-' + group + '-' + atten + '.ns'
            noise_file = os.path.join(gui_data_path, file_)
            # record the noise sample rate
            sample_rate = sweep.metadata.Noise_Sample_Rate
            # record decimation factor
            decfac = sweep.metadata.Noise_Decimation_Factor
            # record the number of points per trace
            n_points_noise = (sweep.metadata.Noise_Integration_Time * sample_rate
                              / decfac)
    else:
        raise ValueError('Expecting {0} .ns and .dat to be in the gui data directory'
                         .format(file_name))
    # grab all of the I and Q noise data
    noise_traces = np.fromfile(noise_file, dtype=np.int16)
    # remove the header from the file
    noise_traces = noise_traces[4 * 12:]
    # convert the data to voltages  * 0.2 V / (2**15 - 1)
    noise_traces = noise_traces.astype(np.float64) * 0.2 / 32767.0
    # calculate the number of number of noise traces and make sure it agrees with GUI
    n_triggers_noise = len(noise_traces) / float(n_points_noise) / 4
    assert float(str(float(n_triggers_noise)).split('.')[-1]) == 0, \
        "non-integer number of noise traces found found in {0}".format(noise_file)
    assert float(str(float(n_points_noise)).split('.')[-1]) == 0, \
        "The noise adtime and sample rate do not give an integer number of " \
        "data points"
    n_triggers_noise = int(n_triggers_noise)
    n_points_noise = int(n_points_noise)
    # break I and Q noise data into their individual channels
    noise_I1 = np.zeros((n_triggers_noise, n_points_noise))
    noise_Q1 = np.zeros((n_triggers_noise, n_points_noise))
    noise_I2 = np.zeros((n_triggers_noise, n_points_noise))
    noise_Q2 = np.zeros((n_triggers_noise, n_points_noise))
    for trigger_num in range(n_triggers_noise):
        trace_num = 4 * trigger_num
        noise_I1[trigger_num, :] = noise_traces[trace_num * n_points_noise:
                                                (trace_num + 1) * n_points_noise]
        noise_Q1[trigger_num, :] = noise_traces[(trace_num + 1) * n_points_noise:
                                                (trace_num + 2) * n_points_noise]
        noise_I2[trigger_num, :] = noise_traces[(trace_num + 2) * n_points_noise:
                                                (trace_num + 3) * n_points_noise]
        noise_Q2[trigger_num, :] = noise_traces[(trace_num + 3) * n_points_noise:
                                                (trace_num + 4) * n_points_noise]
    # resample the traces if using the sweep gui data and needed
    if not use_pulse_gui_psd and (decfac != 1 or sample_rate !=
                                            data_dict['curr_config']['samprate']):
        num = int(data_dict['curr_config']['samprate'] *
                  sweep.metadata.Noise_Integration_Time)
        noise_I1 = resample(noise_I1, num, axis=1)
        noise_Q1 = resample(noise_Q1, num, axis=1)
        noise_I2 = resample(noise_I2, num, axis=1)
        noise_Q2 = resample(noise_Q2, num, axis=1)
        warnings.warn("Noise resampled to match pulse data sample rate", RuntimeWarning)
    # compute average S21 for each channel
    noise_center_1 = np.mean(noise_I1 + 1j * noise_Q1)
    noise_center_2 = np.mean(noise_I2 + 1j * noise_Q2)

    # grab all of the I and Q trace data
    triggers = np.fromfile(data_file, dtype=np.int16)
    # remove the header from the file
    triggers = triggers[4 * 14:]
    # convert the I and Q data to voltages
    triggers = triggers.astype(np.float64) * 0.2 / 32767.0  # 0.2 V / (2**15 - 1)
    # calculate the number of pulses and make sure it's an integer
    n_points = data_dict['curr_config']['numpts']
    n_triggers = len(triggers) / float(n_points) / 4
    assert float(str(float(n_triggers)).split('.')[-1]) == 0, \
        "non-integer number of pulses found in {0}".format(data_file)
    n_triggers = int(n_triggers)
    # break I and Q trace data into their individual channels
    I1 = np.zeros((n_triggers, n_points))
    Q1 = np.zeros((n_triggers, n_points))
    I2 = np.zeros((n_triggers, n_points))
    Q2 = np.zeros((n_triggers, n_points))
    for trigger_num in range(n_triggers):
        trace_num = 4 * trigger_num
        I1[trigger_num, :] = triggers[trace_num * n_points: (trace_num + 1) * n_points]
        Q1[trigger_num, :] = triggers[(trace_num + 1) * n_points:
                                      (trace_num + 2) * n_points]
        I2[trigger_num, :] = triggers[(trace_num + 2) * n_points:
                                      (trace_num + 3) * n_points]
        Q2[trigger_num, :] = triggers[(trace_num + 3) * n_points:
                                      (trace_num + 4) * n_points]
    # add some values to the metadata
    metadata.n_points_per_trigger = n_points
    metadata.n_triggers = n_triggers
    metadata.n_points_per_trigger_noise = n_points_noise
    metadata.n_triggers_noise = n_triggers_noise
    metadata.using_legacy_pulse_gui_noise = use_pulse_gui_psd
    metadata.Trace_Data_Column_Names = ("trace_I", "trace_Q", "noise_I", "noise_Q", "F0",
                                        "noise_center", "sample_rate", "E")

    # make the Pulse_Array
    wavelengths = h.to('eV s').value * c.to('nm/s').value / np.array(energies)
    pulse_data_columns_list, pulse_data_columns = _define_pulse_data_columns_legacy_gui(
        n_points, n_triggers, n_points_noise, n_triggers_noise, len(energies))
    Pulse_Array = np.zeros(2, dtype=pulse_data_columns)
    kwargs = {"trace_I": I1, "trace_Q": Q1, "wavelengths": wavelengths, "E": energies,
              "F0": data_dict['curr_config']['f01'] * 1e9,
              "atten": data_dict['curr_config']['atten1'], "noise_I": noise_I1,
              "noise_Q": noise_Q1, "noise_center": noise_center_1,
              "sample_rate": data_dict['curr_config']['samprate']}
    _define_pulse_array(Pulse_Array, 0, **kwargs)
    kwargs = {"trace_I": I2, "trace_Q": Q2, "wavelengths": wavelengths, "E": energies,
              "F0": data_dict['curr_config']['f02'] * 1e9,
              "atten": data_dict['curr_config']['atten2'], "noise_I": noise_I2,
              "noise_Q": noise_Q2, "noise_center": noise_center_2,
              "sample_rate": data_dict['curr_config']['samprate']}
    _define_pulse_array(Pulse_Array, 1, **kwargs)

    return pulse_data_columns_list, pulse_data_columns, Pulse_Array
