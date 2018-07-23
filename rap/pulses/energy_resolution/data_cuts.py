import numpy as np
import corner


def characterize_data(trace, filter):
    # determine the pulse arrival time
    data = np.array([trace.trace_phase, trace.trace_amp])
    filtered_data = -filter.apply_filter(data)
    # the amplitude estimator is the minimum value of the sum of the filtered traces
    amplitude_estimator = filtered_data[0] + filtered_data[1]
    trace.peak_times = np.argmin(amplitude_estimator, axis=1)

    # determine the mean of the trace prior to the pulse
    trace.prepulse_mean = np.zeros(np.shape(trace.peak_times))
    for index, peak in enumerate(trace.peak_times):
        trace.prepulse_mean[index] = np.mean(trace.trace_phase[index, :peak])

    # determine the rms value of the trace prior to the pulse
    trace.prepulse_rms = np.zeros(np.shape(trace.peak_times))
    for index, peak in enumerate(trace.peak_times):
        if peak == 0:
            trace.prepulse_rms[index] = 0
        else:
            trace.prepulse_rms[index] = np.sqrt(np.mean(
                trace.trace_phase[index, :peak]**2))

    # determine the minimum slope after the pulse peak
    trace.postpulse_min_slope = np.zeros(np.shape(trace.peak_times))
    for index, peak in enumerate(trace.peak_times):
        if peak > len(trace.trace_phase[0, :]) - 2:
            trace.postpulse_min_slope[index] = 0
        else:
            trace.postpulse_min_slope[index] = \
                np.min(np.diff(trace.trace_phase[index, peak:]))


def plot_metrics(trace):
    metrics = np.vstack([trace.peak_times, trace.prepulse_mean, trace.prepulse_rms,
                         trace.postpulse_min_slope]).T
    corner.corner(metrics, labels=['peak times', 'prepulse mean', 'prepulse rms',
                                   'postpulse min slope'],
                  plot_contours=False, plot_density=False, range=[.97, .97, .97, .97],
                  bins=[100, 20, 20, 20])


def cut_peak_time(trace, min, max):
    indices = np.where(np.logical_and(trace.peak_times >= min,
                                      trace.peak_times <= max))[0]
    trace.trace_phase = trace.trace_phase[indices, :]
    trace.trace_amp = trace.trace_amp[indices, :]
    trace.peak_times = trace.peak_times[indices]
    trace.prepulse_mean = trace.prepulse_mean[indices]
    trace.prepulse_rms = trace.prepulse_rms[indices]
    trace.postpulse_min_slope = trace.postpulse_min_slope[indices]


def cut_prepulse_mean(trace, min, max):
    indices = np.where(np.logical_and(trace.prepulse_mean >= min,
                                      trace.prepulse_mean <= max))[0]
    trace.trace_phase = trace.trace_phase[indices, :]
    trace.trace_amp = trace.trace_amp[indices, :]
    trace.peak_times = trace.peak_times[indices]
    trace.prepulse_mean = trace.prepulse_mean[indices]
    trace.prepulse_rms = trace.prepulse_rms[indices]
    trace.postpulse_min_slope = trace.postpulse_min_slope[indices]


def cut_prepulse_rms(trace, min, max):
    indices = np.where(np.logical_and(trace.prepulse_rms >= min,
                                      trace.prepulse_rms <= max))[0]
    trace.trace_phase = trace.trace_phase[indices, :]
    trace.trace_amp = trace.trace_amp[indices, :]
    trace.peak_times = trace.peak_times[indices]
    trace.prepulse_mean = trace.prepulse_mean[indices]
    trace.prepulse_rms = trace.prepulse_rms[indices]
    trace.postpulse_min_slope = trace.postpulse_min_slope[indices]


def cut_postpulse_min_slope(trace, min):
    indices = np.where(trace.postpulse_min_slope >= min)[0]
    trace.trace_phase = trace.trace_phase[indices, :]
    trace.trace_amp = trace.trace_amp[indices, :]
    trace.peak_times = trace.peak_times[indices]
    trace.prepulse_mean = trace.prepulse_mean[indices]
    trace.prepulse_rms = trace.prepulse_rms[indices]
    trace.postpulse_min_slope = trace.postpulse_min_slope[indices]
