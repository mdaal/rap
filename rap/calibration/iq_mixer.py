import numpy as np


def mixer_cal(z, type):
    if type == 'IR0218LC10_765402':
        # JPL large IQ mixer type_serial
        AI = 80.4
        AQ = 79.7
        gamma = 77.5 / 180 * np.pi  # positive phase for this mixer for some reason
        cal = calibration(z, AI, AQ, gamma)
    else:
        raise ValueError('not a valid mixer type')
    return cal


def calibration(z, AI, AQ, gamma):
    g = AI * z.imag / (AQ * z.real)
    # arg1 = z.real >= 0
    # arg2 = z.real < 0
    # theta = np.zeros(np.shape(arg1))
    # theta[arg1] = np.arctan((np.cos(gamma) - g[arg1]) / np.sin(gamma))
    # theta[arg2] = np.arctan((np.cos(gamma) - g[arg2]) / np.sin(gamma)) + np.pi
    theta = np.arctan((np.cos(gamma) - g) / np.sin(gamma)) + np.pi
    theta[theta >= 2 * np.pi] = theta[theta >= 2 * np.pi] - 2 * np.pi
    theta[theta < 0] = 2 * np.pi - theta[theta < 0]

    arg1 = theta >= 1e-3
    arg2 = theta < 1e-3
    r = np.zeros(np.shape(arg1))
    r[arg1] = z.real[arg1] / (AI * np.cos(theta[arg1]))
    r[arg2] = z.imag[arg2] / (AQ * np.cos(theta[arg2] - gamma))
    return r * np.exp(1j * theta)
