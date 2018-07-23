import numpy as np
from scipy.signal import welch, csd


def compute_noise(trace):
    # compute I/Q noise in V^2 / Hz
    n_points = trace.noise_I.shape[1]
    kwargs = {'nperseg': n_points, 'fs': trace.sample_rate, 'return_onesided': True}
    trace.noise_freqs, noise_II = welch(trace.noise_I, **kwargs)
    _, noise_QQ = welch(trace.noise_Q, **kwargs)
    _, noise_IQ = csd(trace.noise_I, trace.noise_Q, **kwargs)

    # average multiple PSDs together
    trace.noise_II = np.mean(noise_II, axis=0)
    trace.noise_QQ = np.mean(noise_QQ, axis=0)
    trace.noise_IQ = np.mean(noise_IQ, axis=0)

    # compute phase and amplitude noise in rad^2 / Hz
    if trace.trace_phase is not None:
        _, noise_PP = welch(trace.noise_phase, **kwargs)
        _, noise_AA = welch(trace.noise_amp, **kwargs)
        _, noise_PA = csd(trace.noise_phase, trace.noise_amp, **kwargs)

        # average multiple PSDs together
        trace.noise_PP = np.mean(noise_PP, axis=0)
        trace.noise_AA = np.mean(noise_AA, axis=0)
        trace.noise_PA = np.mean(noise_PA, axis=0)
