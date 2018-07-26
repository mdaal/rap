import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider


def plot_IQ_trace(trace):
    '''
    Plot the IQ loop and signal for a single trace object
    '''

    fig = plt.figure(figsize=(6, 8))

    ax0 = plt.subplot2grid((4, 1), (0, 0))
    phase, = ax0.plot(trace.time, trace.trace_phase[0, :] * 180 / np.pi)
    ax0.set_ylabel("phase [degrees]")
    ax0.set_title('Channel {}'.format(trace.index + 1))

    ax1 = plt.subplot2grid((4, 1), (1, 0))
    amp, = ax1.plot(trace.time, trace.trace_amp[0, :])
    ax1.set_ylabel("amplitude [loop radius]")
    ax1.set_xlabel(r"time [$\mu s$]")

    ax2 = plt.subplot2grid((4, 1), (2, 0), rowspan=2)
    phi = np.linspace(0, 2 * np.pi, 100)
    ax2.plot(trace.radius * np.cos(phi), trace.radius * np.sin(phi))
    ax2.plot(np.real(trace.z), np.imag(trace.z))
    ax2.plot(0, 0, '+')
    loop, = ax2.plot(trace.trace_I[0, :], trace.trace_Q[0, :], 'o', markersize=2)
    ax2.axis('equal')
    ax2.set_xlabel('I [Volts]')
    ax2.set_ylabel('Q [Volts]')

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)

    class Index(object):
        def __init__(self, ax_slider, ax_prev, ax_next):
            self.ind = 0
            self.num = len(trace.trace_I[:, 0])
            self.bnext = Button(ax_next, 'Next')
            self.bnext.on_clicked(self.next)
            self.bprev = Button(ax_prev, 'Previous')
            self.bprev.on_clicked(self.prev)
            self.slider = Slider(ax_slider, 'Trace Index: ', 0, self.num, valinit=0,
                                 valfmt='%d')
            self.slider.on_changed(self.update)

            position = ax_slider.get_position()
            self.slider.label.set_position((0.5, -0.5))
            self.slider.valtext.set_position((0.5, -0.5))

        def next(self, event):
            self.ind += 1
            i = self.ind % self.num
            self.slider.set_val(i)

        def prev(self, event):
            self.ind -= 1
            i = self.ind % self.num
            self.slider.set_val(i)

        def update(self, value):
            self.ind = int(value)
            i = self.ind % self.num
            loop.set_xdata(trace.trace_I[i, :])
            loop.set_ydata(trace.trace_Q[i, :])
            phase.set_ydata(trace.trace_phase[i, :] * 180 / np.pi)
            amp.set_ydata(trace.trace_amp[i, :])

            ax0.relim()
            ax0.autoscale()
            ax1.relim()
            ax1.autoscale()
            plt.draw()

    position = ax2.get_position()
    ax_slider = plt.axes([position.x0, 0.05, position.width / 2, 0.03])
    middle = position.x0 + 3 * position.width / 4
    ax_prev = plt.axes([middle - 0.18, 0.05, 0.15, 0.03])
    ax_next = plt.axes([middle + 0.02, 0.05, 0.15, 0.03])
    indexer = Index(ax_slider, ax_prev, ax_next)
    plt.show(block=True)
