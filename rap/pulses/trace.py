class trace:
    '''
    The purpose of this class is to hold data associated with pulse traces and
    fitting them
    '''
    def __init__(self):
        # channel number and associated data
        self.index = None
        self.trace_I = None
        self.trace_Q = None
        self.noise_I = None
        self.noise_Q = None
        self.noise_center = None
        self.F0 = None
        self.sample_rate = None
        self.time = None
        # fit data from that channel
        self.radius = None
        self.center = None
        self.cable_delay = None
        self.z = None
        self.freq = None
        # parameters that can be computed
        self.trace_phase = None
        self.trace_amp = None
        self.phase_noise = None
        self.amp_noise = None
        self.noise_freqs = None
        self.noise_II = None
        self.noise_IQ = None
        self.noise_QQ = None
        self.noise_PP = None  # phase noise
        self.noise_PA = None
        self.noise_AA = None  # amplitude noise
        self.calibration = None  # amplitude to energy polynomial
        self.amplitude_fit_result = None
        self.binned_fit_result = None
        self.energy_fit_result = None
        self.energy_fit_model = None
        self.energy_fit_init = None
        self.calibration_type = None

    def copy(self, trace):
        for key in self.__dict__.keys():
            self.__dict__[key] = trace.__dict__[key]

    def __del__(self):
        pass
