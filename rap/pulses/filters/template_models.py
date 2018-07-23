import os
import numpy as np
import lmfit as lm
from numpy import fft as fft
from scipy.integrate import simps


class SingleTypeList(list):
    '''
    Class for containing a list of objects with similar methods and attributes. Methods
    and attributes which are availible for all of the contained objects are availible and
    run on all of the elements sequentially. If there is an output it is wrapped in
    another SingleTypeList.
    '''
    def __init__(self):
        super().__init__()
        self._shape = (0,)
        self._size = 0

    def __getattr__(self, name):
        if name[0] == '_':
            raise AttributeError('Can not call private functions on objects inside a '
                                 ' SingleTypeList')
        # if it's an attribute return a SingleTypeList of attributes
        if not callable(getattr(self[0], name)):
            container = SingleTypeList()
            for index, obj in enumerate(self):
                container.append(getattr(obj, name))
            return container

        # if not, return a container of the method outputs for each object
        def new_method(*args, **kwargs):
            container = SingleTypeList()
            for index, obj in enumerate(self):
                # get the method from the model in the list
                old_method = getattr(obj, name)
                # modify the args
                new_args = [arg[index] for arg in args]
                new_kwargs = {key: value[index] for key, value in kwargs.items()}
                # run the method and put the result in the list
                container.append(old_method(*new_args, **new_kwargs))
            # return None if no output from old_method
            if all(element is None for element in container):
                return None
            return container
        return new_method

    def __str__(self):

        def print_item(string, index):
            try:
                string += 'Item {}:'.format(index) + os.linesep + item.__str__() + \
                          os.linesep
            except:
                string += 'Item {}:'.format(index) + os.linesep + item.__repr__() + \
                          os.linesep
            return string

        string = ''
        for index, item in enumerate(self):
            if len(self) <= 3:
                string = print_item(string, index)
            elif index <= 1:
                string = print_item(string, index)
            elif index == 2:
                string += '...' + os.linesep
            elif index >= len(self) - 2:
                string = print_item(string, index)

        return string

    @property
    def shape(self):
        has_shape = all([hasattr(self[ind], 'shape') for ind in range(len(self))])
        if has_shape:
            same_shape = all([self[ind].shape == self[0].shape
                              for ind in range(len(self))])
            if same_shape:
                self._shape = (len(self), *self[0].shape)
        if (not has_shape) or (not same_shape):
            self._shape = (len(self),)
        return self._shape

    @property
    def size(self):
        has_size = all([hasattr(self[ind], 'size') for ind in range(len(self))])
        if has_size:
            same_size = all([self[ind].size == self[0].size
                             for ind in range(len(self))])
            if same_size:
                self._size = len(self) * self[0].size
            else:
                self._size = np.sum([self[ind].size for ind in range(len(self))])
        if not has_size:
            self._size = len(self)
        return self._size


class ModelND(SingleTypeList):
    '''
    SingleTypeList subclass for a N-D list of lmfit Model objects. Has all the methods of
    the models in model_list but all arguements to those functions need to be lists with
    the length N. This class acts just like a lmfit Model class.
    '''
    def __init__(self, independent_vars=['x'], prefix='', nan_policy='raise',
                 model_list=None, **kwargs):
        if model_list is None:
            raise ValueError('must provide a list of lmfit model classes')
        super().__init__()
        kwargs.update({'prefix': prefix, 'nan_policy': nan_policy,
                       'independent_vars': independent_vars})
        self.extend([model(**kwargs) for model in model_list])


class TripleExponentialPulseModel(lm.Model):
    '''
    lmfit style model for a pulse with an exponential rise time and a double exponential
    fall time.
    '''
    def __init__(self, independent_vars=['x'], prefix='', nan_policy='raise', **kwargs):
        kwargs.update({'prefix': prefix, 'nan_policy': nan_policy,
                       'independent_vars': independent_vars})
        super().__init__(self._pulse, **kwargs)

    def guess(self, data, x):
        dt = x[1] - x[0]
        max_tau = np.max(x) / 2
        a = np.abs(np.min(data))
        t0 = x[np.argmin(data)]
        tau0 = np.min([np.min(np.diff(x)), max_tau])
        tau = -simps(data, x=x) / a
        tau1 = np.min([tau / 5, max_tau])
        tau2 = np.min([tau, max_tau])
        alpha = 1

        initial = self._pulse(x, a, t0, tau0, tau1, tau2, alpha)
        t0 = x[np.argmin(data) - (np.argmin(initial) - np.argmin(data))]

        params = lm.Parameters()
        params.add('a', value=a, min=0, max=np.inf)
        params.add('t0', value=t0, min=0, max=np.max(x))
        params.add('tau0', value=tau0, min=0, max=max_tau)
        params.add('tau1', value=tau1, min=dt, max=max_tau)
        params.add('tau2', value=tau2, min=dt, max=max_tau)
        params.add('alpha', value=alpha, min=0, max=np.inf)

        return params

    @staticmethod
    def _pulse(x, a, t0, tau0, tau1, tau2, alpha):
        # resample x to give more precision to the amplidue calculation
        xp = np.linspace(x[0], x[-1], 10 * len(x))
        pulse_shape = np.zeros(xp.shape)
        t = xp[xp >= t0]
        pulse_shape[xp >= t0] = -((1 - np.exp(-(t - t0) / tau0)) *
                                  (alpha * np.exp(-(t - t0) / tau1) +
                                   np.exp(-(t - t0) / tau2)))
        # normalize to max pulse height of 'a'
        norm = np.abs(np.min(pulse_shape))
        if norm == 0:
            norm = 1
        pulse_shape *= a / norm
        # down sample back to original sampling rate
        pulse_shape = pulse_shape[::10]
        return pulse_shape


class LinearFrequencyResponseModel(lm.Model):
    def __init__(self, independent_vars=['x'], prefix='', nan_policy='raise', **kwargs):
        kwargs.update({'prefix': prefix, 'nan_policy': nan_policy,
                       'independent_vars': independent_vars})
        super().__init__(self._pulse, **kwargs)

    def guess(self, data, x, bandwidth=None):
        max_tau = np.max(x) / 2
        dt = x[1] - x[0]
        a = np.abs(np.min(data))
        t0 = x[np.argmin(data)]
        tau = np.min([-simps(data, x=x), max_tau]) / a
        tau0 = tau / 5
        tau1 = tau
        alpha = 1
        f = fft.rfftfreq(len(x), d=dt)
        fmax = np.max(f)
        fmin = np.min(np.abs(f))
        if bandwidth is None:
            bandwidth = 2 / 3 * fmax

        initial = self._pulse(x, a, t0, tau0, tau1, bandwidth, alpha)
        t0 = x[np.argmin(data) - (np.argmin(initial) - np.argmin(data))]

        params = lm.Parameters()
        params.add('a', value=a, min=0, max=np.inf)
        params.add('t0', value=t0, min=0, max=np.max(x))
        params.add('tau0', value=tau0, min=0, max=max_tau)
        params.add('tau1', value=tau1, min=0, max=max_tau)
        params.add('bandwidth', value=bandwidth, min=fmin, max=fmax)
        params.add('alpha', value=alpha, min=0, max=np.inf)

        return params

    @staticmethod
    def _pulse(x, a, t0, tau0, tau1, bandwidth, alpha):
        def single_pole_filter(f):
            return 1 / (1 + 1j * f / bandwidth)
        # resample x to give more precision in fft calculation
        xp = np.linspace(x[0], x[-1], 10 * len(x))
        # basic pulse shape for df/fr
        # phase is proportional to df/fr for small signals
        pulse_shape = np.zeros(xp.shape)
        t = xp[xp >= t0]
        pulse_shape[xp >= t0] = -(alpha * np.exp(-(t - t0) / tau0) +
                                  np.exp(-(t - t0) / tau1))
        # frequency transfer function
        f = fft.rfftfreq(len(xp), d=xp[1] - xp[0])
        transfer = single_pole_filter(f)
        # apply the transfer function to get the non-adiabatic response
        pulse_fft = fft.rfft(pulse_shape)
        pulse_shape = fft.irfft(pulse_fft * transfer, len(xp))
        # normalize to max pulse height of 'a'
        norm = np.abs(np.min(pulse_shape))
        if norm == 0:
            norm = 1
        pulse_shape *= a / norm
        # down sample back to original sampling rate
        pulse_shape = pulse_shape[::10]
        return pulse_shape


class LinearDissipationResponseModel(lm.Model):
    pass
