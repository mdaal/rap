import numpy as np
from numpy import linalg as la
from scipy import signal as sg
from numpy import fft as fft


class Filter(object):
    """
    Class to hold the pulse filter and filter manipulation methods.
    """
    def __init__(self, template):
        if template.template_fit is None:
            raise ValueError('A template fit must be computed first')
        if template.trace_phase is None or template.trace_amp is None:
            raise ValueError('Phase and amplitude traces need to be calculated first.')

        self.template_fit = template.template_fit
        self.sample_rate = template.sample_rate

        self.noise_PP = template.noise_PP
        self.noise_PA = template.noise_PA
        self.noise_AA = template.noise_AA
        self.noise_freqs = template.noise_freqs

    def make_filter(self):
        # pull out shape parameters
        N = len(self.template_fit[0])
        shape = (2, N // 2 + 1)

        # compute template fft
        self.template_fft = np.zeros(shape, dtype=np.complex)
        self.template_fft = fft.rfft(self.template_fit)
        self.fft_freqs = fft.rfftfreq(N, d=1 / self.sample_rate)

        # assemble noise matrix
        self.S = np.array([[self.noise_PP, self.noise_PA],
                           [np.conj(self.noise_PA), self.noise_AA]])

        # compute the optimal filter: conj(template_fft) @ Sinv
        self.filter_fft = np.zeros(shape, dtype=np.complex)
        for index in range(shape[1]):
            self.filter_fft[:, index] = la.lstsq(self.S[:, :, index].T,
                                                 np.conj(self.template_fft[:, index]),
                                                 rcond=None)[0]
        # expected variance (multiplied by 4 because of single sided psd and fft)
        self.expected_variance = (np.sum((self.filter_fft @
                                          self.template_fft.T)).real /
                                  (self.sample_rate * N) * 4)**-1
        # return to time domain
        self.filter = np.zeros((2, N))
        self.filter[0, :] = fft.irfft(self.filter_fft[0, :], N)
        self.filter[1, :] = fft.irfft(self.filter_fft[1, :], N)
        # normalize the optimal filter
        norm1 = sg.convolve(self.filter[0], self.template_fit[0], mode='valid')
        norm2 = sg.convolve(self.filter[1], self.template_fit[1], mode='valid')
        norm = norm1 + norm2
        self.filter[0, :] /= norm
        self.filter[1, :] /= norm

        # compute the phase only optimal filter: conj(phase_fft) / J
        self.phase_filter_fft = (np.conj(self.template_fft[0, :]) / self.noise_PP)
        self.phase_filter = fft.irfft(self.phase_filter_fft, N)
        # normalize
        norm = sg.convolve(self.phase_filter, self.template_fit[0], mode='valid')
        self.phase_filter /= norm
        # expected variance (multiplied by 4 because of single sided psd and fft)
        self.phase_expected_variance = (np.sum(np.abs(self.template_fft[0, :])**2 /
                                               self.noise_PP) /
                                        (self.sample_rate * N) * 4)**-1

    def apply_filter(self, data, phase_only=False):
        """
        Method for convolving the two dimensional filter with the data. The data can
        either be a 2xN matrix or a 2xMxN matrix where N is the trace length and M is the
        number of traces.
        """
        shape = np.shape(data)
        if phase_only:
            if shape == np.shape(self.phase_filter):
                result = sg.convolve(self.phase_filter, data, mode='same')
            elif len(shape) == 2 and shape[1] == len(self.phase_filter):
                result = np.zeros(shape)
                for index in range(shape[0]):
                    result[index, :] = sg.convolve(self.phase_filter, data[index, :],
                                                   mode='same')
            else:
                raise ValueError("data needs to be a 1 or 2D array with one dimension "
                                 "equal in length to the filter length")
            return result

        if shape == np.shape(self.filter):
            result = np.array([sg.convolve(self.filter[0], data[0], mode='same'),
                               sg.convolve(self.filter[1], data[1], mode='same')])
        elif len(shape) == 3 and shape[0] == 2 and shape[2] == len(self.filter[0]):
            result = np.zeros(shape)
            for index in range(shape[1]):
                result[0, index, :] = sg.convolve(self.filter[0], data[0, index, :],
                                                  mode='same')
                result[1, index, :] = sg.convolve(self.filter[1], data[1, index, :],
                                                  mode='same')
        else:
            raise ValueError("data needs to be a 2 x N x M array "
                             "(last dimension optional)")
        return result
