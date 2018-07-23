import numpy as np
from scipy.special import erfc, erfcx
import warnings


sqrt2 = np.sqrt(2)
sqrt2pi = np.sqrt(2 * np.pi)


def expgaussian(x, amplitude, center, sigma, gamma, a=0, b=1, c=0):
    """
    Return an exponentially modified Gaussian distribution.
    Use gamma = np.inf to return the usual Gaussian distribution

    Has opposite skew from:
    https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution

    a, b, c are optional polynomial coefficients which convert from amplitude to energy
    """
    try:
        result = np.zeros(x.shape)
    except AttributeError:
        x = np.array(x)
        result = np.zeros(x.shape)

    # transform with a particular calibration if provided
    x = a * x**2 + b * x + c

    z = (sigma * gamma + (x - center) / sigma) / sqrt2
    logic1 = (z < 0)
    logic2 = np.logical_and(z >= 0, z < 6.71e7)
    logic3 = (z >= 6.71e7)

    if logic1.any():
        arg1 = (sigma * gamma)**2 / 2 + (x[logic1] - center) * gamma
        arg2 = z[logic1]
        result[logic1] = amplitude * gamma / 2 * np.exp(arg1) * erfc(arg2)

    if logic2.any():
        arg1 = -0.5 * ((x[logic2] - center) / sigma)**2
        arg2 = z[logic2]
        result[logic2] = amplitude * gamma / 2 * np.exp(arg1) * erfcx(arg2)

    if logic3.any():
        arg1 = -0.5 * ((x[logic3] - center) / sigma) ** 2
        A = (amplitude / (sqrt2pi * np.abs(sigma)) /
             (1 + (x[logic3] - center) / (sigma**2 * gamma)))
        result[logic3] = A * np.exp(arg1)

    return result


def expgaussian2(x, amplitude, center, sigma, gamma, a=0, b=1, c=0):
    """
    Return an exponentially modified Gaussian distribution.
    Use gamma = np.inf to return the usual Gaussian distribution

    Has the same skew as:
    https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution

    a, b, c are optional polynomial coefficients which convert from amplitude to energy
    """
    try:
        result = np.zeros(x.shape)
    except AttributeError:
        x = np.array(x)
        result = np.zeros(x.shape)

    # transform with a particular calibration if provided
    x = a * x**2 + b * x + c

    z = (sigma * gamma - (x - center) / sigma) / sqrt2
    logic1 = (z < 0)
    logic2 = np.logical_and(z >= 0, z < 6.71e7)
    logic3 = (z >= 6.71e7)

    if logic1.any():
        arg1 = (sigma * gamma)**2 / 2 - (x[logic1] - center) * gamma
        arg2 = z[logic1]
        result[logic1] = amplitude * gamma / 2 * np.exp(arg1) * erfc(arg2)

    if logic2.any():
        arg1 = -0.5 * ((x[logic2] - center) / sigma)**2
        arg2 = z[logic2]
        result[logic2] = amplitude * gamma / 2 * np.exp(arg1) * erfcx(arg2)

    if logic3.any():
        arg1 = -0.5 * ((x[logic3] - center) / sigma) ** 2
        A = (amplitude / (sqrt2pi * np.abs(sigma)) /
             (1 + (x[logic3] - center) / (sigma**2 * gamma)))
        result[logic3] = A * np.exp(arg1)

    return result


def gaussian(x, amplitude, center, sigma, a=0, b=1, c=0):
    """
    Return a 1-dimensional Gaussian function.

    a, b, c are optional polynomial coefficients which convert from amplitude to energy
    """
    # transform with a particular calibration if provided
    x = a * x ** 2 + b * x + c
    arg = -0.5 * ((x - center) / sigma) ** 2
    return (amplitude / (sqrt2pi * np.abs(sigma))) * np.exp(arg)