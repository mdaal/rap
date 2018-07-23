import numpy as np
import numpy.ma as ma
import numpy.fft as fft
from scipy.integrate import simps

from matplotlib import pyplot as plt

from .template_models import (ModelND, TripleExponentialPulseModel,
                              LinearFrequencyResponseModel)
from .filter import Filter

from ..energy_resolution.compute_noise import compute_noise


class Template(object):
    '''
    Class to hold the pulse template and template manipulation methods.
    '''
    def __init__(self, trace):
        # extract needed info from trace object
        if trace.trace_phase is None or trace.trace_amp is None:
            raise ValueError('Phase and amplitude traces need to be calculated first.')
        self.trace_phase = trace.trace_phase
        self.trace_amp = trace.trace_amp
        self.sample_rate = trace.sample_rate
        self.time = trace.time
        self.data = np.array([trace.trace_phase, trace.trace_amp])

        # compute noise if it hasn't been done
        if trace.noise_PP is None:
            compute_noise(trace)
        self.noise_PP = trace.noise_PP
        self.noise_PA = trace.noise_PA
        self.noise_AA = trace.noise_AA
        self.noise_freqs = trace.noise_freqs

        self.template = None
        self.template_fit = None

    def make_template(self, Show_Plot=False):
        '''
        Make a template from phase and amplitude data. This method sets the template and
        template fit attributes, which are both 2xN arrays with the phase in the first
        row and amplitude in the second.
        '''
        # create a rough template by cutting the noise traces and averaging
        self.reset()
        self._threshold_cut()
        self._average_pulses()
        self._fit_template()
        # make a filter with the template
        filter = Filter(self)
        filter.make_filter()
        # do a better job using a filter
        self.reset()
        self.filtered_phase = filter.apply_filter(self.data[0], phase_only=True)
        self._threshold_cut(filter=True)
        self._single_pulse_cut()
        self._offset_correction()
        self._average_pulses()
        self._fit_template()

        if Show_Plot:
            self.plot_fit()

    def reset(self):
        '''
        Reset the data, template and template_fit attributes.
        '''
        self.data = np.array([self.trace_phase, self.trace_amp])
        self.template = None
        self.template_fit = None

    def plot_fit(self):
        '''
        Plot the template fit
        '''
        time = self.time * 1e6
        t0 = [self.result.params[0]['t0'], self.result.params[1]['t0']]
        # collect all the time constants from each dimension into a list
        tau_0 = []
        for key, value in self.result.params[0].items():
            if len(key) >= 3 and key[:3] == 'tau':
                tau_0.append(value.value)
        tau_1 = []
        for key, value in self.result.params[1].items():
            if len(key) >= 3 and key[:3] == 'tau':
                tau_1.append(value.value)
        tau = [np.sum(tau_0) * 1e6, np.sum(tau_1) * 1e6]

        fig = plt.figure(figsize=(12, 6))
        ax0 = plt.subplot2grid((3, 2), (0, 0), rowspan=2)
        ax0.plot(time, self.template[0], 'b-', linewidth=2, label='template data')
        ax0.plot(time, self.result.init_fit[0], 'k--', linewidth=1,
                 label='template guess')
        ax0.plot(time, self.template_fit[0], 'r-', linewidth=1, label='template fit')
        ax0.set_ylabel('phase')

        ax1 = plt.subplot2grid((3, 2), (2, 0), sharex=ax0)
        res = (self.template[0] - self.template_fit[0]) / np.max(np.abs(self.template[0]))
        ax1.plot([time[0], time[-1]], [0, 0], 'k', linewidth=1)
        ax1.plot(time, res, 'b-', linewidth=2)
        ax1.set_xlabel(r'time [$\mu$s]')
        ax1.set_ylabel('normalized residuals')
        ax1.set_xlim([np.median(time) - 2 * tau[0], np.median(time) + 5 * tau[0]])
        ax1.text(1, 0, 'rms = {:.3e}'.format(np.sqrt(np.mean(res**2))),
                 ha='right', va='bottom', transform=ax1.transAxes, fontsize=12)

        ax2 = plt.subplot2grid((3, 2), (0, 1), rowspan=2)
        ax2.plot(time, self.template[1], 'b-', linewidth=2, label='template data')
        ax2.plot(time, self.result.init_fit[1], 'k--', linewidth=1,
                 label='template guess')
        ax2.plot(time, self.template_fit[1], 'r-', linewidth=1, label='template fit')
        ax2.set_ylabel('amplitude')
        ax2.legend()

        ax3 = plt.subplot2grid((3, 2), (2, 1), sharex=ax2)
        res = (self.template[1] - self.template_fit[1]) / np.max(np.abs(self.template[1]))
        ax3.plot([time[0], time[-1]], [0, 0], 'k', linewidth=1)
        ax3.plot(time, res, 'b-', linewidth=2)
        ax3.set_xlabel(r'time [$\mu$s]')
        ax3.set_xlim([np.median(time) - 2 * tau[1], np.median(time) + 5 * tau[1]])
        ax3.text(1, 0, 'rms = {:.3e}'.format(np.sqrt(np.mean(res**2))),
                 ha='right', va='bottom', transform=ax3.transAxes, fontsize=12)

        ax0.tick_params(labelbottom=False)
        ax2.tick_params(labelbottom=False)

        plt.tight_layout()
        plt.show()

    def _threshold_cut(self, filter=False, threshold=5):
        '''
        Remove traces from the data object that don't meet the threshold condition on the
        phase trace.

        threshold is the number of standard deviations to put the threshold cut.
        '''
        if filter:
            data = np.array([-self.filtered_phase, self.data[1]])
        else:
            data = self.data
        # Compute the median average deviation use that to calculate the standard
        # deviation. This should be robust against the outliers from pulses.
        median_phase = np.median(data[0], axis=1)[:, np.newaxis]
        mad = np.median(np.abs(data[0] - median_phase))
        sigma = 1.4826 * mad

        # look for a phase < sigma * threshold around the middle of the trace
        n_points = len(data[0, 0, :])
        middle = (n_points + 1) / 2
        start = int(middle - 25)
        stop = int(middle + 25)
        phase_middle = data[0, :, start:stop]
        indices = np.where((phase_middle - median_phase < -sigma * threshold).any(axis=1))
        self.data = self.data[:, indices[0], :]
        if self.data.size == 0:
            raise RuntimeError('All data was removed by cuts')
        if filter:
            self.filtered_phase = self.filtered_phase[indices[0], :]

    def _single_pulse_cut(self):
        '''
        Remove traces with more than one pulse.
        '''
        # pull out filtered data
        data = np.array([-self.filtered_phase, self.data[1]])
        # determine pulse time constant
        tau = np.max([int(round(self._tau)), 4])
        # find the peak max index for each pulse
        t0 = np.argmin(data[0], axis=1)[:, np.newaxis]
        # find the indices in [t0-tau/2, t0 + 3 * tau] for each pulse
        indices = np.array([np.arange(0, len(data[0, 0, :]))] * len(data[0, :, 0]))
        mask = np.logical_and(indices > t0 - tau / 2, indices < t0 + 3 * tau)

        # mask the peak from the data (only use phase)
        data = ma.array(data[0], mask=mask)

        # Compute the median average deviation use that to calculate the standard
        # deviation. This should be robust against the outliers from pulses.
        median_phase = ma.median(data, axis=1)[:, np.newaxis]
        mad = ma.median(ma.abs(data - median_phase))
        sigma = 1.4826 * mad

        # look for a phase < sigma * threshold and exclude those points
        indices = np.where(np.logical_not((data - median_phase < -4 * sigma).any(axis=1)))
        self.data = self.data[:, indices[0], :]
        self.filtered_phase = self.filtered_phase[indices[0], :]
        if self.data.size == 0:
            raise RuntimeError('All data was removed by cuts')

    def _offset_correction(self, fractional_offset=True):
        '''
        Correct for trigger offsets.
        '''
        # pull out filtered data
        data = np.array([-self.filtered_phase, self.data[1]])

        # define frequency vector
        f = fft.rfftfreq(len(data[0, 0, :]))

        # find the middle index
        middle_ind = int(np.median(range(len(data[0, 0, :]))))

        # loop through triggers and shift the peaks
        for index in range(len(data[0, :, 0])):
            # pull out phase and amplitude
            phase_trace = self.data[0, index, :]
            amp_trace = self.data[1, index, :]
            # find the peak index of filtered phase data
            peak_ind = np.argmin(data[0, index, 2:-2]) + 2
            # determine the shift
            if fractional_offset:
                peak = data[0, index, peak_ind - 2: peak_ind + 3]
                poly = np.polyfit([-2, -1, 0, 1, 2], peak, 2)
                frac = -poly[1] / (2 * poly[0])
                shift = middle_ind - peak_ind - frac
            else:
                shift = middle_ind - peak_ind
            # remove the shift by applying a phase
            # same as np.roll(data, shift) for integer offsets
            self.data[0, index, :] = fft.irfft(fft.rfft(phase_trace) *
                                               np.exp(-1j * 2 * np.pi * shift * f),
                                               len(phase_trace))
            self.data[1, index, :] = fft.irfft(fft.rfft(amp_trace) *
                                               np.exp(-1j * 2 * np.pi * shift * f),
                                               len(amp_trace))

    def _average_pulses(self):
        '''
        Average the data together by summing and normalizing the phase pulse height to 1.
        '''
        # add all the data together
        self.template = np.sum(self.data, axis=1)
        # remove any small baseline error
        self.template -= np.median(self.template, axis=1)[:, np.newaxis]
        # normalize phase pulse height to 1
        self.template /= np.abs(np.min(self.template[0]))

    def _fit_template(self, model_types=[0, 0]):
        '''
        Fit a model to the template data
        '''
        # define model list and time array
        models = np.array([TripleExponentialPulseModel, LinearFrequencyResponseModel])
        time = [self.time] * 2
        # create model and parameter objects
        model = ModelND(model_list=models[model_types])
        params = model.guess(self.template, time)
        # coarse fit that has better convergence properties
        initial_result = model.fit(self.template, params, x=time, method='powell')
        # default fit method so that errors are computed
        result = model.fit(self.template, initial_result.params, x=time)
        # reset initial parameter attributes
        result.init_fit = initial_result.init_fit
        result.init_vals = initial_result.init_vals
        result.init_values = initial_result.init_values
        # save result and template fit
        self.template_fit = result.best_fit
        self.result = result
        # save a time constant for internal use
        self._tau = simps(self.template_fit[0]) / np.min(self.template_fit[0])
