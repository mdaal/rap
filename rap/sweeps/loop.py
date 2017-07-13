class loop:
	'''The purpose of this class is to hold data associated with resonance loops and fitting them'''
	def __init__(self):
		self.index = None
		self.z =  None
		self.freq = None

		#output of circle fit
		self.r = None
		self.a = None #fit circle center is a+i*b
		self.b = None
		self.outer_iterations = None
		self.inner_iterations = None

		#circle fit parameters
		self.s = None
		self.Gx = None
		self.Gy = None
		self.g = None
		self.sigma = None
		self.circle_fit_exit_code = None

		#loop fit estimates quantities
		self.fr_est = None
		self.FWHM_est = None
		self.depth_est = None
		self.Q_est = None

		# intermediate fit quantities
		self.normalization = 1# used for nonlinear of  generated data probably should  get rid of
		

		# phase fit quantities
		self.Q = None 
		self.Qc = None
		self.Qi = None
		self.fr = None 
		self.FWHM = None
		self.phi = None # in  radians not degrees
		self.theta = None # in  radians not degrees
		self.R = None # outer loop radius
		self.chisquare = None
		self.pvalue = None
		self.phase_fit_method = None
		self.phase_fit_success = None
		self.phase_fit_z = None
		self.phase_fit_mask = None

		# complete fit quantities
		# self.c___ is the phase/magnitude fit
		self.cQ = None 
		self.cQc = None
		self.cQi = None
		self.cfr = None 
		self.cphi = None
		self.cR = None
		self.ctheta = None
		self.cchisquare = None
		self.cphase_fit_success = None

		# self.s___ is the self fit result where sigma^2 is determined from data
		self.sQ = None 
		self.sQc = None
		self.sQi = None
		self.sfr = None 
		self.sphi = None
		self.sR = None
		self.stheta = None
		self.schisquare = None
		self.sphase_fit_success = None


	def __del__(self):
		pass