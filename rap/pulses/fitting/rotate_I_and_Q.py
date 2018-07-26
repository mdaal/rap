import numpy as np
import matplotlib.pyplot as plt


def rotate_I_and_Q(trace):
    # combine I and Q traces for the signal and noise so they can be rotated
    signal = np.array([trace.trace_I, trace.trace_Q]).swapaxes(0, 1)
    noise = np.array([trace.noise_I, trace.noise_Q]).swapaxes(0, 1)
    # rotate if a cable delay has been calculated
    if trace.cable_delay is not None:
        phi = 2 * np.pi * trace.F0 * trace.cable_delay
        R = np.array([[np.cos(phi), -np.sin(phi)],
                      [np.sin(phi), np.cos(phi)]])
        signal = R @ signal
        noise = R @ noise
        trace.noise_center *= np.exp(1j * phi)

    # center based on the loop fit
    rotated_signal = np.zeros(signal.shape)
    rotated_signal[:, 0, :] = signal[:, 0, :] - trace.center.real
    rotated_signal[:, 1, :] = signal[:, 1, :] - trace.center.imag
    rotated_noise = np.zeros(noise.shape)
    rotated_noise[:, 0, :] = noise[:, 0, :] - trace.center.real
    rotated_noise[:, 1, :] = noise[:, 1, :] - trace.center.imag
    trace.z -= trace.center
    trace.noise_center -= trace.center

    # rotate the baseline to the x-axis
    phi = np.arctan2(trace.noise_center.imag, trace.noise_center.real)
    R = np.array([[np.cos(phi), np.sin(phi)],
                  [-np.sin(phi), np.cos(phi)]])
    rotated_signal = R @ rotated_signal
    rotated_noise = R @ rotated_noise
    trace.z *= np.exp(-1j * phi)
    trace.noise_center *= np.exp(-1j * phi)

    # calculate the phase and amplitude signals
    phase = np.arctan2(rotated_signal[:, 1, :], rotated_signal[:, 0, :])
    phase_noise = np.arctan2(rotated_noise[:, 1, :], rotated_noise[:, 0, :])
    amplitude = (np.linalg.norm(rotated_signal, axis=1) - trace.radius) / trace.radius
    amp_noise = (np.linalg.norm(rotated_noise, axis=1) - trace.radius) / trace.radius

    # remove any remaining baseline in the noise by subtracting the median
    ind = int(np.floor(len(phase[0, :]) / 5))
    phase -= np.median(phase[:, :ind], axis=1)[:, np.newaxis]
    phase_noise -= np.median(phase_noise, axis=1)[:, np.newaxis]
    amplitude -= np.median(amplitude[:, :ind], axis=1)[:, np.newaxis]
    amp_noise -= np.median(amp_noise, axis=1)[:, np.newaxis]

    # save data in the trace object
    trace.trace_I = rotated_signal[:, 0, :]
    trace.trace_Q = rotated_signal[:, 1, :]
    trace.noise_I = rotated_noise[:, 0, :]
    trace.noise_Q = rotated_noise[:, 1, :]
    trace.trace_phase = phase
    trace.trace_amp = amplitude
    trace.noise_phase = phase_noise
    trace.noise_amp = amp_noise
