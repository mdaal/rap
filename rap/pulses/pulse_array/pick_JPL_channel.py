import numpy as np


def pick_JPL_channel(metadata, trace, sweep, Pulse_Array):
    """
    Use this function to pick the channel to manipulate data from the Pulse_Array.
    """
    sweep.pick_JPL_loop()
    trace_I, trace_Q, noise_I, noise_Q, F0, noise_center, sample_rate, energies = \
        metadata.Trace_Data_Column_Names

    trace.index = 0
    trace.trace_I = Pulse_Array[0][trace_I]
    trace.trace_Q = Pulse_Array[0][trace_Q]
    trace.noise_I = Pulse_Array[0][noise_I]
    trace.noise_Q = Pulse_Array[0][noise_Q]
    trace.F0 = Pulse_Array[0][F0]
    trace.noise_center = Pulse_Array[0][noise_center]
    trace.sample_rate = Pulse_Array[0][sample_rate]
    trace.energies = Pulse_Array[0][energies]

    n_points = len(trace.trace_I[0, :])
    trace.time = np.linspace(0, (n_points - 1) / trace.sample_rate, n_points)