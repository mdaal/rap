import numpy as np
import lmfit as lm
from scipy.integrate import simps
from ..fitting.functions import (gaussian, expgaussian, expgaussian2)

sqrt2 = np.sqrt(2)
sqrt2pi = np.sqrt(2 * np.pi)


def calculate_amplitudes(trace, filter, phase_only=False):
    if phase_only:
        amplitude_estimator = -filter.apply_filter(trace.trace_phase, phase_only=True)
    else:
        data = np.array([trace.trace_phase, trace.trace_amp])
        data = -filter.apply_filter(data)
        # the amplitude estimator is the minimum value of the sum of the filtered traces
        amplitude_estimator = data[0, :, :] + data[1, :, :]
    # find the middle index near where the pulse is
    middle_ind = int(np.median(range(len(amplitude_estimator[0, :]))))
    amplitudes = np.min(amplitude_estimator[:, middle_ind - 10: middle_ind + 10], axis=1)

    return amplitudes


def amplitude_fit(amplitudes, energies, n_bins=100, fit_type='gaussian',
                  calibration='linear', **parameter_kwargs):
    # recast energies as a numpy array (just to be sure)
    energies = np.array(energies)

    # bin the data
    amp_counts, amp_bins = np.histogram(amplitudes, n_bins)
    bin_width = amp_bins[1] - amp_bins[0]
    amp_centers = (amp_bins + bin_width / 2)[:-1]

    # remove any bins with zero counts
    indices = (amp_counts != 0)
    amp_counts = amp_counts[indices]
    amp_centers = amp_centers[indices]

    # modified gaussian error for low count rates
    # https://doi.org/10.1016/S0168-9002(99)00269-7
    sigma_amp_counts = np.sqrt(amp_counts + 0.25) - 0.5

    # pick the fitting function
    models = {'gaussian': gaussian, 'expgaussian': expgaussian2}
    model_functions = []
    if type(fit_type) is list or type(fit_type) is np.ndarray:
        if type(fit_type) is list:
            fit_type = np.array(fit_type)
        for func in fit_type:
            model_functions.append(models[func])
    else:
        model_functions = [models[fit_type]] * energies.size
        fit_type = np.array([fit_type] * energies.size)

    # create the model
    prefixes = []
    for index, _ in enumerate(energies):
        prefixes.append('d{:d}_'.format(index))
        if index == 0:
            model = lm.Model(model_functions[index], prefix=prefixes[-1])
        else:
            model += lm.Model(model_functions[index], prefix=prefixes[-1])

    # create the model initial parameters
    params = model.make_params()
    set_gamma = False
    for index, prefix in enumerate(prefixes):
        # assume equal heights for the peaks
        amplitude = np.sum(amp_counts) * bin_width / len(prefixes)
        params[prefix + 'amplitude'].set(value=amplitude)
        # space the Gaussians between max and min bins
        max_bin = np.max(amp_centers)
        min_bin = np.min(amp_centers)
        center = min_bin + (max_bin - min_bin) / (energies.size + 2) * (index + 1)
        params[prefix + 'center'].set(value=center)
        # assume a typical R = 8 standard deviation
        sigma = energies[index] / (2.355 * 8)
        params[prefix + 'sigma'].set(value=sigma, min=1e-15)
        # fix the calibration to have no effect.
        params[prefix + 'a'].set(value=0, vary=False)
        params[prefix + 'b'].set(value=1, vary=False)
        params[prefix + 'c'].set(value=0, vary=False)
        # If gamma parameter exists start with a large gamma (gaussian assumption)
        # and only allow one to vary by default
        if prefix + 'gamma' in params.keys():
            if not set_gamma:
                params[prefix + 'gamma'].set(value=50, min=1e-15)
                set_gamma = True
                set_index = index
            else:
                params[prefix + 'gamma'].set(expr=prefixes[set_index] + 'gamma')
        # add a derived resolving power parameter
        if fit_type[index] == 'gaussian':
            scale = 2 * np.sqrt(2 * np.log(2))
        elif fit_type[index] == 'expgaussian':
            scale = 2 * np.sqrt(2 * np.log(2))
        elif fit_type[index] == 'lorentzian':
            scale = 2
        else:
            ValueError("scale not yet defined in code for fit_type {}"
                       .format(fit_type[index]))
        R_string = '-' + prefix + 'center / ({} * ' + prefix + 'sigma)'
        params.add(prefix + 'R', expr=R_string.format(scale))

    # update the parameters with any changes supplied
    for param, value in parameter_kwargs.items():
        params[param].set(**value)

    # fit the model to the data
    result = model.fit(amp_counts, params, x=amp_centers, weights=(1 / sigma_amp_counts),
                       scale_covar=False, nan_policy='omit')

    centers = np.array([result.params[prefix + 'center'] for prefix in prefixes])
    centers = np.sort(centers)
    energies = np.sort(energies)[::-1]
    if calibration == 'linear':
        calibration = np.polyfit(centers, energies, 1)
        calibration = np.append(0, calibration)
    elif calibration == 'quadratic':
        calibration = np.polyfit(centers, energies, 2)
    else:
        raise ValueError("calibration parameter must be linear or quadratic")

    return result, calibration


def binned_energy_fit(amplitudes, energies, n_bins=100, fit_type='gaussian',
                      calibration='linear', **parameter_kwargs):
    # recast energies as a numpy array (just to be sure)
    energies = np.array(energies)

    # bin the data
    amp_counts, amp_bins = np.histogram(amplitudes, n_bins)
    bin_width = amp_bins[1] - amp_bins[0]
    amp_centers = (amp_bins + bin_width / 2)[:-1]

    # remove any bins with zero counts
    indices = (amp_counts != 0)
    amp_counts = amp_counts[indices]
    amp_centers = amp_centers[indices]

    # modified gaussian error for low count rates
    # https://doi.org/10.1016/S0168-9002(99)00269-7
    sigma_amp_counts = np.sqrt(amp_counts + 0.25) - 0.5

    # pick the fitting function
    models = {'gaussian': gaussian, 'expgaussian': expgaussian}
    model_functions = []
    if type(fit_type) is list or type(fit_type) is np.ndarray:
        if type(fit_type) is list:
            fit_type = np.array(fit_type)
        for func in fit_type:
            model_functions.append(models[func])
    else:
        model_functions = [models[fit_type]] * energies.size
        fit_type = np.array([fit_type] * energies.size)

    # create the model
    prefixes = []
    for index, _ in enumerate(energies):
        prefixes.append('d{:d}_'.format(index))
        if index == 0:
            model = lm.Model(model_functions[index], prefix=prefixes[-1])
        else:
            model += lm.Model(model_functions[index], prefix=prefixes[-1])

    # create the model initial parameters
    params = model.make_params()
    set_gamma = False
    for index, prefix in enumerate(prefixes):
        # assume equal heights for the peaks
        amplitude = np.sum(amp_counts) * bin_width / len(prefixes)
        params[prefix + 'amplitude'].set(value=amplitude)
        # fix the centers at the given energies
        params[prefix + 'center'].set(value=energies[index], vary=False)
        # assume a typical R = 8 standard deviation
        sigma = energies[index] / (2.355 * 8)
        params[prefix + 'sigma'].set(value=sigma, min=1e-15)
        # set the calibration to be the same for each distribution
        # and linear with a slope of -1
        if index == 0:
            if calibration == 'linear':
                vary = False
            elif calibration == 'quadratic':
                vary = True
            else:
                raise ValueError("calibration parameter must be linear or quadratic")
            params[prefix + 'a'].set(value=0, vary=vary)
            params[prefix + 'b'].set(value=-1)
            params[prefix + 'c'].set(value=0)
        else:
            params[prefix + 'a'].set(expr=prefixes[0] + 'a')
            params[prefix + 'b'].set(expr=prefixes[0] + 'b')
            params[prefix + 'c'].set(expr=prefixes[0] + 'c')
        # If gamma parameter exists start with a large gamma (gaussian assumption)
        # and only allow one to vary by default
        if prefix + 'gamma' in params.keys():
            if not set_gamma:
                params[prefix + 'gamma'].set(value=50, min=1e-15)
                set_gamma = True
                set_index = index
            else:
                params[prefix + 'gamma'].set(expr=prefixes[set_index] + 'gamma')
        # add a derived resolving power parameter
        if fit_type[index] == 'gaussian':
            scale = 2 * np.sqrt(2 * np.log(2))
        elif fit_type[index] == 'expgaussian':
            scale = 2 * np.sqrt(2 * np.log(2))
        elif fit_type[index] == 'lorentzian':
            scale = 2
        else:
            ValueError("scale not yet defined in code for fit_type {}"
                       .format(fit_type[index]))
        R_string = prefix + 'center / ({} * ' + prefix + 'sigma)'
        params.add(prefix + 'R', expr=R_string.format(scale))

    # update the parameters with any changes supplied
    for param, value in parameter_kwargs.items():
        params[param].set(**value)

    # fit the model to the data
    result = model.fit(amp_counts, params, x=amp_centers, weights=(1 / sigma_amp_counts),
                       scale_covar=False, nan_policy='omit')

    calibration = [result.params['d0_a'].value, result.params['d0_b'].value,
                   result.params['d0_c'].value]
    return result, calibration


def energy_fit(amplitudes, energies, fit_type='gaussian', calibration='linear',
               **parameter_kwargs):
    # pick the fitting function
    models = {'gaussian': gaussian, 'expgaussian': expgaussian}
    model_functions = []
    if type(fit_type) is list or type(fit_type) is np.ndarray:
        if type(fit_type) is list:
            fit_type = np.array(fit_type)
        for func in fit_type:
            model_functions.append(models[func])
    else:
        model_functions = [models[fit_type]] * energies.size
        fit_type = np.array([fit_type] * energies.size)

    # create the model
    prefixes = []
    for index, _ in enumerate(energies):
        prefixes.append('d{:d}_'.format(index))
        if index == 0:
            model = lm.Model(model_functions[index], prefix=prefixes[-1])
        else:
            model += lm.Model(model_functions[index], prefix=prefixes[-1])

    # create the model initial parameters
    params = model.make_params()
    set_gamma = False
    for index, prefix in enumerate(prefixes):
        # assume equal heights for the peaks
        amplitude = len(amplitudes) / len(prefixes)
        params[prefix + 'amplitude'].set(value=amplitude)
        # fix the centers at the given energies
        params[prefix + 'center'].set(value=energies[index], vary=False)
        # assume a typical R = 8 standard deviation
        sigma = energies[index] / (2.355 * 8)
        params[prefix + 'sigma'].set(value=sigma, min=1e-15)
        # set the calibration to be the same for each distribution
        # and linear with a slope of -1
        if index == 0:
            if calibration == 'linear':
                vary = False
            elif calibration == 'quadratic':
                vary = True
            else:
                raise ValueError("calibration parameter must be linear or quadratic")
            params[prefix + 'a'].set(value=0, vary=vary)
            params[prefix + 'b'].set(value=-1)
            params[prefix + 'c'].set(value=0)
        else:
            params[prefix + 'a'].set(expr=prefixes[0] + 'a')
            params[prefix + 'b'].set(expr=prefixes[0] + 'b')
            params[prefix + 'c'].set(expr=prefixes[0] + 'c')
        # If gamma parameter exists start with a large gamma (gaussian assumption)
        # and only allow one to vary by default
        if prefix + 'gamma' in params.keys():
            if not set_gamma:
                params[prefix + 'gamma'].set(value=50, min=1e-15)
                set_gamma = True
                set_index = index
            else:
                params[prefix + 'gamma'].set(expr=prefixes[set_index] + 'gamma')
        # add a derived resolving power parameter
        if fit_type[index] == 'gaussian':
            scale = 2 * np.sqrt(2 * np.log(2))
        elif fit_type[index] == 'expgaussian':
            scale = 2 * np.sqrt(2 * np.log(2))
        elif fit_type[index] == 'lorentzian':
            scale = 2
        else:
            ValueError("scale not yet defined in code for fit_type {}"
                       .format(fit_type[index]))
        R_string = prefix + 'center / ({} * ' + prefix + 'sigma)'
        params.add(prefix + 'R', expr=R_string.format(scale))

    # update the parameters with any changes supplied
    for param, value in parameter_kwargs.items():
        params[param].set(**value)

    # create the chi squared function that accepts detector response amplitudes
    def chi_sq(parameters, x):
        p = parameters.valuesdict()
        dist = model.eval(parameters, x=x)
        xx = np.linspace(min(x) - 2 * np.abs(min(x)),
                         max(x) + 2 * np.abs(max(x)), 1000)
        A = simps(model.eval(parameters, x=xx), xx)
        chi2 = 2 * (-np.sum(np.log(dist)) + A)

        return chi2

    # run the minimizer
    mini = lm.Minimizer(chi_sq, params, fcn_args=[amplitudes], scale_covar=False)
    result = mini.scalar_minimize()

    calibration = [result.params['d0_a'].value, result.params['d0_b'].value,
                   result.params['d0_c'].value]

    return result, calibration, model, params
