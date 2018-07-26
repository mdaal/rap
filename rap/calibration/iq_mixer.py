import numpy as np


def mixer_cal(z, type):
    if type == 'IR0218LC10_765402':
        # JPL large IQ mixer type_serial
        alpha = 80.4 / 79.7
        # I and Q swapped for this mixer for some reason
        gamma = 77.5 / 180 * np.pi + np.pi / 2
        cal = calibration(z, alpha, gamma)
    else:
        raise ValueError('not a valid mixer type')
    return cal


def calibration(z, alpha, gamma):
    z_corrected = (z.real / alpha +
                   1j * (-z.real * np.tan(gamma) / alpha + z.imag / np.cos(gamma)))
    return z_corrected
