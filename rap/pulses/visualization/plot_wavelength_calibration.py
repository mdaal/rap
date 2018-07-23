import numpy as np
from scipy.integrate import simps
from matplotlib import pyplot as plt


def plot_amplitude_fit(result, amplitudes, n_bins, plot_guess=True):
    prop_cycle = iter(plt.rcParams['axes.prop_cycle'])
    # bin the data
    amp_counts, amp_bins = np.histogram(amplitudes, n_bins)
    bin_width = amp_bins[1] - amp_bins[0]
    amp_centers = (amp_bins + bin_width / 2)[:-1]
    sigma_amp_counts = np.sqrt(amp_counts + 0.25) - 0.5

    plt.figure(figsize=[9, 6])
    plt.bar(amp_centers, amp_counts, width=1.1 * bin_width,
            yerr=sigma_amp_counts, color=next(prop_cycle)['color'])
    xx = np.linspace(min(amp_bins), max(amp_bins), 1000)
    plt.plot(xx, result.eval(x=xx), 'm-', color=next(prop_cycle)['color'], label="fit",
             zorder=3)
    components = result.eval_components(x=xx)
    prefixes = [component.prefix for component in result.components]
    try:
        for prefix in prefixes:
            plt.plot(xx, components[prefix], '--', color=next(prop_cycle)['color'],
                     label=r'$R = {:.2f} \pm {:.2f}$'
                     .format(result.params[prefix + 'R'].value,
                             result.params[prefix + 'R'].stderr))
    except AttributeError:
        for prefix in prefixes:
            plt.plot(xx, components[prefix], '--', color=next(prop_cycle)['color'],
                     label=r'$R = {:.2f}$'
                     .format(result.params[prefix + 'R'].value))
    if plot_guess:
        plt.plot(amp_centers[amp_counts != 0], result.init_fit, 'k-.', label="guess")

    try:
        d_counts = result.eval_uncertainty(sigma=1, x=xx)
        plt.fill_between(xx, result.eval(x=xx) - d_counts, result.eval(x=xx) + d_counts,
                         color="#ABABAB", label=r'1$\sigma$ confidence interval', zorder=2,
                         alpha=0.7)
    except Exception as err:
        print(err)
    plt.legend(loc=1)
    plt.ylabel('Counts')
    plt.xlabel('Detector Response')
    plt.title('Amplitude Estimation Fit')
    plt.xlim([amp_centers[0] - bin_width, amp_centers[-1] + bin_width])
    plt.ylim([0, 1.05 * max(amp_counts + sigma_amp_counts)])
    plt.show()


def plot_binned_energy_fit(result, amplitudes, n_bins, calibration, plot_guess=True):
    prop_cycle = iter(plt.rcParams['axes.prop_cycle'])
    # bin the data
    amp_counts, amp_bins = np.histogram(amplitudes, n_bins)
    bin_width = amp_bins[1] - amp_bins[0]
    amp_centers = (amp_bins + bin_width / 2)[:-1]
    energy_centers = np.polyval(calibration, amp_centers)
    diff = np.diff(energy_centers)
    bin_widths = np.diff(np.cumsum(np.append(np.append(diff[0], diff), diff[-1])))
    sigma_amp_counts = np.sqrt(amp_counts + 0.25) - 0.5

    plt.figure(figsize=[9, 6])
    plt.bar(energy_centers, amp_counts, width=1.1 * bin_widths,
            yerr=sigma_amp_counts, color=next(prop_cycle)['color'])
    xx = np.linspace(min(amp_bins), max(amp_bins), 1000)
    xxx = np.polyval(calibration, xx)
    plt.plot(xxx, result.eval(x=xx), color=next(prop_cycle)['color'] , label="fit",
             zorder=3)
    components = result.eval_components(x=xx)
    prefixes = [component.prefix for component in result.components]
    try:
        for prefix in prefixes:
            plt.plot(xxx, components[prefix], '--', color=next(prop_cycle)['color'],
                     label=r'$R = {:.2f} \pm {:.2f}$'
                     .format(result.params[prefix + 'R'].value,
                             result.params[prefix + 'R'].stderr))
    except AttributeError:
        for prefix in prefixes:
            plt.plot(xxx, components[prefix], '--', color=next(prop_cycle)['color'],
                     label=r'$R = {:.2f}$'
                     .format(result.params[prefix + 'R'].value))
    if plot_guess:
        plt.plot(energy_centers[amp_counts != 0], result.init_fit, 'k-.', label="guess")

    try:
        d_counts = result.eval_uncertainty(sigma=1, x=xx)
        plt.fill_between(xxx, result.eval(x=xx) - d_counts, result.eval(x=xx) + d_counts,
                         color="#ABABAB", label=r'1$\sigma$ confidence interval', zorder=2,
                         alpha=0.7)
    except Exception as err:
        print(err)
    plt.legend(loc=1)
    plt.ylabel('Counts')
    plt.xlabel('Energy [eV]')
    plt.title('Binned Gaussian Maximum Likelihood Fit')
    plt.xlim([energy_centers[-1] - bin_widths[-1], energy_centers[0] + bin_widths[0]])
    plt.ylim([0, 1.05 * max(amp_counts + sigma_amp_counts)])
    plt.show()


def plot_energy_fit(result, amplitudes, n_bins, model, initial_params, plot_guess=True):
    prop_cycle = iter(plt.rcParams['axes.prop_cycle'])
    # bin the data
    calibration = [result.params['d0_a'], result.params['d0_b'], result.params['d0_c']]
    energies = np.polyval(calibration, amplitudes)
    energy_counts, bin_edges = np.histogram(energies, n_bins)
    width = np.diff(bin_edges)[0]
    energy_bins = bin_edges + width / 2
    energy_bins = energy_bins[:-1]
    sigma_energy_counts = np.sqrt(energy_counts + 0.25) - 0.5

    plt.figure(figsize=[9, 6])
    plt.bar(energy_bins, energy_counts, width=1.1 * width, yerr=sigma_energy_counts,
            color=next(prop_cycle)['color'])
    xx = np.linspace(min(amplitudes), max(amplitudes), 1000)
    xxx = np.polyval(calibration, xx)
    norm = (np.abs(simps(model.eval(result.params, x=xx), xxx)) /
            simps(model.eval(result.params, x=xx), xx) / width)
    plt.plot(xxx, model.eval(result.params, x=xx) / norm, color=next(prop_cycle)['color'],
             label="fit", zorder=3)
    components = model.eval_components(params=result.params, x=xx)
    for prefix, data in components.items():
        plt.plot(xxx, data / norm, '--', color=next(prop_cycle)['color'],
                 label=r'$R = {:.2f}$'.format(result.params[prefix + 'R'].value))
    if plot_guess:
        plt.plot(xxx, model.eval(params=initial_params, x=xx) /
                 norm, 'k-.', label="guess")

    plt.legend(loc=1)
    plt.ylabel('Counts')
    plt.xlabel('Energy [eV]')
    plt.title('Binfree Poisson Maximum Likelihood Fit')
    plt.xlim([energy_bins[0] - width, energy_bins[-1] + width])
    plt.show()
