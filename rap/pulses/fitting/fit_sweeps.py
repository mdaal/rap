from ..data_management.utils import _define_pulse_array
from .rotate_I_and_Q import rotate_I_and_Q


def fit_sweeps(pulse, Show_Plot=False, Fit_Method="circle_and_phase", Verbose=False):
    # assumes sweep and trace object have had the right loop/channel picked

    # fit loops with Fit_Method
    if Fit_Method is "circle_and_phase":
        pulse.sweep.remove_cable_delay(Show_Plot=Show_Plot, Verbose=Verbose)
        pulse.sweep.trim_loop(N=1.1, Verbose=Verbose)
        pulse.sweep.circle_fit(Show_Plot=Show_Plot)
        pulse.sweep.phase_fit(Show_Plot=Show_Plot, Verbose=Verbose)
        pulse.trace.cable_delay = pulse.sweep.metadata.Electrical_Delay
    elif Fit_Method is None:
        # allow for fitting the sweep object manually but still update the
        # trace object
        pass
    else:
        raise ValueError("Fit_Method keyword not valid")

    # update trace object with relevant values
    pulse.trace.radius = pulse.sweep.loop.r
    pulse.trace.center = pulse.sweep.loop.a + 1j * pulse.sweep.loop.b
    pulse.trace.freq = pulse.sweep.loop.freq
    pulse.trace.z = pulse.sweep.loop.z
    pulse.trace.Q = pulse.sweep.loop.Q
    pulse.trace.Qc = pulse.sweep.loop.Qc
