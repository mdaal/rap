import numpy as np


def pick_legacy_pulse_gui_channel(metadata, trace, sweep, Pulse_Array, index):
    """
    Use this function to pick the channel to manipulate data from the Pulse_Array.
    Index is channel number to be selected.
    """
    sweep.pick_loop(index)
    trace_I, trace_Q, noise_I, noise_Q, F0, noise_center, sample_rate, energies = \
        metadata.Trace_Data_Column_Names

    trace.index = index
    trace.trace_I = Pulse_Array[index][trace_I]
    trace.trace_Q = Pulse_Array[index][trace_Q]
    trace.noise_I = Pulse_Array[index][noise_I]
    trace.noise_Q = Pulse_Array[index][noise_Q]
    trace.F0 = Pulse_Array[index][F0]
    trace.noise_center = Pulse_Array[index][noise_center]
    trace.sample_rate = Pulse_Array[index][sample_rate]
    trace.energies = PulseArray[index][energies]

    n_points = len(trace.trace_I[0, :])
    trace.time = np.linspace(0, (n_points - 1) / trace.sample_rate, n_points)
