import numpy as np
from astropy.constants import h, c
from scipy.interpolate import interp1d
from .utils import (_define_pulse_array, _define_pulse_data_columns_JPL)
from ...calibration.iq_mixer import mixer_cal as calibration


def load_JPL_pulse_data(metadata, sweep, energies, sweep_file, noise_file, pulse_file,
                        iq_mixer=None):
    # load in the sweep
    sweep.load_JPL_sweep_data(sweep_file, iq_mixer=iq_mixer)

    # pull out the loop data
    f_loop = sweep.Sweep_Array[0]["frequencies"]
    z_loop = sweep.Sweep_Array[0]["s21"]
    z0 = sweep.Sweep_Array[0]["s21_baseline"]
    z_loop_cal = interp1d(f_loop, z_loop)
    z0_cal = interp1d(f_loop, z0)

    # define IQ mixer calibration
    if iq_mixer is not None:
        def mixer_cal(z, f):
            return calibration(z - z0_cal(f), iq_mixer)
    else:
        def mixer_cal(z, f):
            return z

    # load noise and correct for iq offset
    npzfile = np.load(noise_file)
    noise_dict = npzfile['meta'].item()
    F0 = noise_dict['frequency'] * 1e9
    z_res = z_loop_cal(F0)
    z_offset = z0_cal(F0)
    noise = mixer_cal(npzfile['I_traces'] + 1j * npzfile['Q_traces'], F0)
    noise_I = noise.real
    noise_Q = noise.imag

    # correct for voltage drift offset in each trace
    # dI = np.mean(noise_I, axis=1)[:, np.newaxis] - z_res.real
    # dQ = np.mean(noise_Q, axis=1)[:, np.newaxis] - z_res.imag
    # noise_I -= dI
    # noise_Q -= dQ
    noise_center = np.mean(noise_I + 1j * noise_Q)

    # load pulses and correct for iq offset
    npzfile = np.load(pulse_file)
    pulse_dict = npzfile['meta'].item()
    trace = mixer_cal(npzfile['I_traces'] + 1j * npzfile['Q_traces'], F0)
    trace_I = trace.real
    trace_Q = trace.imag

    # correct for voltage drift offset in each trace
    # dI = np.median(trace_I, axis=1)[:, np.newaxis] - z_res.real
    # dQ = np.median(trace_Q, axis=1)[:, np.newaxis] - z_res.imag
    # trace_I -= dI
    # trace_Q -= dQ

    # check some things
    assert pulse_dict['frequency'] == noise_dict['frequency'], \
        "noise and pulse files at different frequencies"
    assert pulse_dict['attenuation'] == noise_dict['attenuation'], \
        "noise and pulse files at different attenuations"
    assert pulse_dict['power'] == noise_dict['power'], \
        "noise and pulse files at different frequencies"
    assert pulse_dict['sample_rate'] == noise_dict['sample_rate'], \
        "noise and pulse files at different attenuations"

    # add some values to the metadata
    n_triggers = trace_I.shape[0]
    n_points = trace_I.shape[1]
    n_triggers_noise = noise_I.shape[0]
    n_points_noise = noise_I.shape[1]
    metadata.n_points_per_trigger = n_points
    metadata.n_triggers = n_triggers

    metadata.summary_string = pulse_dict['summary']
    metadata.n_points_per_trigger_noise = n_points_noise
    metadata.n_triggers_noise = n_triggers_noise
    metadata.noise_sample_rate = noise_dict['sample_rate']
    metadata.noise_summary_string = noise_dict['summary']
    metadata.Trace_Data_Column_Names = ("trace_I", "trace_Q", "noise_I", "noise_Q", "F0",
                                        "noise_center", "sample_rate", "E")
    # create the pulse array
    wavelengths = h.to('eV s').value * c.to('nm/s').value / np.array(energies)
    pulse_data_columns_list, pulse_data_columns = _define_pulse_data_columns_JPL(
        n_points, n_triggers, n_points_noise, n_triggers_noise, len(energies))
    Pulse_Array = np.zeros(1, dtype=pulse_data_columns)
    kwargs = {"trace_I": trace_I, "trace_Q": trace_Q, "wavelengths": wavelengths,
              "E": energies, "F0": F0, "atten": pulse_dict['attenuation'],
              "noise_I": noise_I, "noise_Q": noise_Q, "noise_center": noise_center,
              "sample_rate": pulse_dict['sample_rate']}
    _define_pulse_array(Pulse_Array, 0, **kwargs)

    return pulse_data_columns_list, pulse_data_columns, Pulse_Array
