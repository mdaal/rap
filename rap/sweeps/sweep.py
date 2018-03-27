

### super-directory imports
from ..metadata import metadata
from ..env_var import env_var

### current directory imports 
from .loop import loop
from .thermometry import thermometry

### subdirectory imports

### data_management
from .data_management.construct_hf5_toc import construct_hf5_toc
from .data_management.load_hf5 import load_hf5
from .data_management.load_hf5_2 import load_hf5_2
from .data_management.load_scandata import load_scandata
from .data_management.load_touchstone import load_touchstone
from .data_management.save_hf5 import save_hf5
from .data_management.load_legacy_sweep_gui_data import load_legacy_sweep_gui_data
# from .data_management.utils import _read_scandata_from_file
# from .data_management.utils import _download_data
# from .data_management.utils import _extract_type
# from .data_management.utils import _define_sweep_data_columns
# from .data_management.utils import _define_sweep_array

###------------------------------------------------------------
### fitting
from .fitting.circle_fit import circle_fit
from .fitting.complete_fit import complete_fit
from .fitting.decompress_gain import decompress_gain
from .fitting.downsample_loop import downsample_loop
from .fitting.normalize_loop import normalize_loop
from .fitting.phase_fit import phase_fit
from .fitting.remove_cable_delay import remove_cable_delay
from .fitting.trim_loop import trim_loop
# from .fitting.utils import _points_removed
# from .fitting.utils import _angle

from .fitting.nonlinear.nonlinear_fit import nonlinear_fit
# from .fitting.nonlinear.utils import _nonlinear_formulae

###------------------------------------------------------------
### simulation
from .simulation.generate_nonlinear_data import generate_nonlinear_data

###------------------------------------------------------------
### sweep_array
from .sweep_array.fill_sweep_array import fill_sweep_array
from .sweep_array.sweep_array_info import sweep_array_info
from .sweep_array.pick_loop import pick_loop
from .sweep_array.pick_legacy_sweep_gui_loop import pick_legacy_sweep_gui_loop
###------------------------------------------------------------
### system_calibration

from .system_calibration.fit_cable_loss import fit_cable_loss
from .system_calibration.fit_system_calibration import fit_system_calibration
# from system_calibration.utils import _construct_readout_chain

###------------------------------------------------------------
### tests


###------------------------------------------------------------
### visualization 
# from visualization._save_fig_dec import _save_fig_dec
from .visualization.plot_loop import plot_loop
from .visualization.plot_transmission import plot_transmission

###------------------------------------------------------------
### external packages
import numpy as np
import matplotlib.pyplot as plt
import os


class sweep():
	'''This class accesses sweep data and fits resonances'''
	
	
	def __init__(self, **kword):
		'''
		keywords can be:
		database_location = 'path/to/data.h5'
		'''
		self.loop = loop()
		self.metadata = metadata()	
		self.env_var = env_var(**kword)



	#def _read_scandata_from_file(self,filename_or_path): Not needed


	#def _download_data(self, URL): Not needed
		

	def plot_loop(self,  aspect='equal', show = True):
		''' Plots currently selected complex transmission in the I,Q plane. Reutrns a tuple, (fig, ax, line),
		where  fig is the figure object, ax is the axes object and line is the line object for the plotted data.

		aspect='equal' makes circles round, aspect='auto' fills the figure space.

		*Must have a loop picked in order to use this function.*
		'''

		(fig, ax, line) = plot_loop(self.metadata, self.loop,  aspect, show)

		return  (fig, ax, line)


	def plot_transmission(self, show = True):
		''' Plots currently selected complex transmission in dB as a function of frequency. Reutrns a tuple, (fig, ax, line),
		where fig is the figure object, ax is the axes object and line is the line object for the plotted data.

		*Must have a loop picked in order to use this function.*
		'''
		(fig, ax, line) = plot_transmission(self.metadata, self.loop, show)

		return  (fig, ax, line)

	# def _extract_type(self, obj, return_type = None, field = None): Not needed

	# def _define_sweep_data_columns(self, fsteps, tpoints, list_only = False): not needed

	# def _define_sweep_array(self,index,**field_names): Not needed


	def load_scandata(self, file_location, Verbose = True, **auth):
		''' file_location is the locaiton of the scandata.mat file. It can be a URL, filename or /path/filename.
		assumes that self.data is in the form of matlab ScanData Structure

		creates Sweep_Array from data contained in ScanData.mat ... The UC Berkeley data product'''

		#delete previous metadata object
		del(self.metadata)
		self.metadata = metadata()

		self.sweep_data_columns_list, self.sweep_data_columns, self.Sweep_Array = load_scandata(self.metadata, file_location, Verbose  = Verbose, **auth)

		try:
			metadata.Cable_Calibration = self._Cable_Calibration
			print('Cable Calibration data found and saved in Sweep_Array metadata.')
		except:
			pass

		try:
			metadata.Temperature_Calibration = self._Temperature_Calibration
			print('Temperature Calibration data found and saved in Sweep_Array metadata.')
		except:
			pass


	def load_touchstone(self,filename, pick_loop = True):
		''' The function loads S21 and Freq from  Sonnet .s2p or .s3p files into the Sweep_Array structured np array
		All Sij are extracted, but only  S21 is saved into Sweep_Array. Future editions of this code might  find need 
		to load other Sij becuase S21.

		The function only loads one transmission array (S21).  pick_loop = True immediatly selectes this loop as the 
		current loop.
		'''

		#delete previous metadata object
		del(self.metadata)
		self.metadata = metadata()

		self.sweep_data_columns_list, self.sweep_data_columns, self.Sweep_Array = load_touchstone(self.metadata, filename)

		if pick_loop == True: #since there is only one loop in Sweep_Array, we might as well pick it as the current loop
			self.pick_loop(0)
			#self.normalize_loop()


	def load_legacy_sweep_gui_data(self, gui_data_path):
		'''
		Import data generated by Matlab sweep_gui progam into rap data structures, e.g. Sweep_Array, metadata, sweep() objects.
		gui_data_path is the director path containing the sweep_gui program generated date.  This function finds the latest 
		sweep_config_xxxxx.mat file in director, then reads this file to dermine from which specXX-YY-ZZ.mat and spec_offresXX-YY-ZZ.mat 
		files to obtain data.

		### As of 9/13/17: still need to implement CPSD and specorth storage in Sweep_Array
		'''

		#delete previous metadata object
		del(self.metadata)
		self.metadata = metadata()

		self.sweep_data_columns_list, self.sweep_data_columns, self.Sweep_Array = load_legacy_sweep_gui_data(self.metadata, gui_data_path)


	def downsample_loop(self,N):
		''' Reduce number of loop/freq data point by every Nth point and discarding all others'''
		downsample_loop(self.loop, N)
		

	def save_hf5(self, overwrite = False):
		'''Saves current self.Sweep_Array into table contained in the hdf5 file speficied by filename.
		If overwite = True, self.Sweep_Array will overwright whatever is previous table data there is.
		'''
		save_hf5(self.metadata,self.Sweep_Array, self.sweep_data_columns,  self.env_var.database_location, overwrite)
	

	def decompress_gain(self, Compression_Calibration_Index = -1, Show_Plot = True, Verbose = True):
		''' Assumes the two lowest input powers of the power sweep are not gain compressed, thus
		cannot be used if the two lowest powers are gain compressed. '''
		decompress_gain(self.Sweep_Array, self.loop, self.metadata, Compression_Calibration_Index, Show_Plot, Verbose)


	def sweep_array_info(self):
		''' prints information about the Sweep_Array currently loaded'''
		sweep_array_info(self.Sweep_Array)


	def construct_hf5_toc(self):
		''' Creates a table of contents (toc) of the hf5 database storing all the sweep_data.
		very useful for finding the name and locaiton of a table in the database'''
		construct_hf5_toc(self.env_var.database_location)

	
	def load_hf5_2(self, database_filename):
		''' This function is for loading data taken with KIDs_DAQ_75uW. It use the columns defined in that hf5 file to 
		define the columns  in  self.sweep_data_columns 
		table path is path to the database to be loaded starting from root. e.g. self.load_hf5('/Run44b/T201312102229')
		database_filename is the name of the hf5 database to be accessed for the  table informaiton'''
		
		#delete previous metadata and loop objects
		del(self.metadata)
		self.metadata = metadata()
		del(self.loop)
		self.loop = loop()

		try:
			self.metadata.Cable_Calibration = self._Cable_Calibration
			print('Cable Calibraiton data found and saved in Sweep_Array metadata.')
		except:
			pass

		self.Sweep_Array, self.sweep_data_columns, self.sweep_data_columns_list = load_hf5_2(self.metadata, self.env_var.database_filename, tablepath)


	def load_hf5(self, tablepath):
		''' table path is path to the database to be loaded starting from root. e.g. self.load_hf5('/Run44b/T201312102229')
		filename is the name of the hf5 database to be accessed for the  table informaiton'''

		#delete previous metadata  and loop objects
		del(self.metadata)
		self.metadata = metadata()
		del(self.loop)
		self.loop = loop()

		self.Sweep_Array, self.sweep_data_columns = load_hf5(self.metadata, self.env_var.database_filename, tablepath)


	def pick_loop(self,index):
		'''Use this function to pick the current loop/transmission data from withing the Sweep_Array. 
		Index is the indes number of sweep/loop to be slected as the current loop.'''
		pick_loop(self.metadata, self.loop, self.Sweep_Array,index)
		
	def pick_legacy_sweep_gui_loop(self, index_tuple):
		'''When Sweep_Array is generated by load_legacy_sweep_gui_data, 
		this function can be used to select the current loop/transmission data from withing the Sweep_Array.
		index_tuple follows legacy gui matrlab file nameing convention... 

				i.e.  specXX-YY-ZZ.mat has data for resonator pair i = 0 or 1.
				XX is the Temp in milliKelvin
				YY is group number of the pair
				ZZ is the attenuator setting going toward the cryostat
				----> index_tuple = (XX, YY, ZZ, i) and i is 0 or 1

		issue the command obj.pick_legacy_sweep_gui_loop( (XX, YY, ZZ, i) ) to load that loop.'''
		pick_legacy_sweep_gui_loop(self.metadata, self.loop, self.Sweep_Array,index_tuple)

	def normalize_loop(self, base = 0, offset = 5):
		''' normalize loop so that mag(S21)< 1. determine normalization by averaging np.abs(S21[base:offset]).mean()
		return normalization'''
		
		normalization = normalize_loop(self.loop, base, offset)
		return normalization 


	def remove_cable_delay(self, Show_Plot = True, Verbose = True, center_freq = None, Force_Recalculate = False):
		'''
		If self.metadate.Electrical_Delay is not None, then use this value as cable delay and remove from data 

		If self.metadate.Electrical_Delay is None:
		- Determine cable delay by finding delay value, tau, which minimizes distance between adjacent S21 points. 
		Then cancel out tau in S21 data and save corrected loop in self.loop.z. Set self.metadate.Electrical_Delay = tau.
		- If S21 is large array, down sample it first before performing minimization

		If self.metadate.Electrical_Delay is None and center_freq is given:
		-If center_freq is given, this function computes the electrical delay by determining the bandwidth over which the S21 
		circle completes a full revolution starting at center_freq and ending at ~ center_freq + tau^-1. Where tau is approximated
		as the vaule deterined by minimum  distance above. 
		-center_freq should only be used when S21 is is sufficiently broadband to generously cover center_freq and ~center_freq + tau^-1.
		center_freq is in Hertz.

		Return tau in any case.

		If Force_Recalculate == True Electrical delay will be recalculated and reset in metadata
		'''

		cable_delay = remove_cable_delay(self.loop, self.metadata, Show_Plot, Verbose, center_freq, Force_Recalculate)
		return  cable_delay


	def trim_loop(self,N = 20,Verbose = True,Use_Dip = True):
		''' Find resonance frequency of transmission dip, estimate FWHM, then cut all point
		greater than N*FWHM on upper and lower side of transmission min. 

		do this to masked array.

		IF Use_Dip = True,  estimate Resonance frequency using minimum Dip 
		if Use_Dip = False, estimate Resonance frequencymax adjacent distance
		'''

		#add error if N < 1
		trim_loop(self.loop, N, Verbose, Use_Dip)

	# def _points_removed(self,initial, final): Not needed

	def circle_fit(self, Show_Plot = True):
		
		circle_fit(self.loop) #circle_fit edits the class copy of loop!
		
		if Show_Plot:
			fig, ax = self.plot_loop(show = False)[:2]		
			t = np.linspace(0, 2.0*np.pi, num=50, endpoint=True)
			j = np.complex(0,1); zc = self.loop.a + j*self.loop.b;  r = self.loop.r
			line = ax.plot(zc.real + r*np.cos(t),zc.imag + r*np.sin(t),'y-', label = 'Circle Fit')
			line = ax.plot([zc.real],[zc.imag],'yx', markersize = 10, markeredgewidth = 4, label = 'Center')
			ax.set_aspect('equal')
			plt.show()


	def phase_fit(self, Fit_Method = 'Multiple', Verbose = True, Show_Plot = True):
		'''
		Note: its best to determine angles and angle differences by starting with complex numbers 
		(interpreted as vectors) and then finding their angles with, np.angle or self._angle. It is
		not as accurate and prone to issues with domains (e.g. [-180,180]) to use arcsin or arccos.
		'''
		phase_fit(self.loop, self.env_var, Fit_Method, Verbose, Show_Plot)


	def fill_sweep_array(self, Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False, Complete_Fit = True , Remove_Gain_Compression = True, Verbose = True):
		
		fill_sweep_array(metadata, Sweep_Array, Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False, Complete_Fit = True , Remove_Gain_Compression = True, Verbose = True)

	# def _construct_readout_chain(self, F, Include_NA = True, Include_4K_to_40mK = False): Not Needed
	
	def complete_fit(self, Use_Mask = True, Verbose = False , Show_Plot = False, Save_Fig = False, Sample_Size = 100, Use_Loop_Data = False):
		'''
		Sample_Size is the number of points used to extablish \sigma^2 for the gaussian noise model

		if Use_Loop_Data = True then values of Q, Qc, fr, phi are for initial guess are taken from curret loop object. If false, values come from self.Sweep_Array
		'''

		fit, plot_dict = complete_fit(self.Sweep_Array, self.metadata, self.loop, Use_Mask, Verbose, Show_Plot, Save_Fig, Sample_Size, Use_Loop_Data)

		return fit, plot_dict


	# def _angle(self, z, deg = 0, return_offset = False): Not Needed


	def fit_system_calibration(self):
		'''compute chebyshev polynomial fits for  gain and noise values.
		save resulting polynomial coefficients list as:
		
		self.metadata.System_Calibration['device'][x + '_fit']

		where x is [gain, Tn_m ,Tn_p]... 

		use numpy.polynomial.chebyshev.chebval to evaluate polynomial
		'''
		fit_system_calibration(self.metadata)

			
	def fit_cable_loss(self, key_name, freq_range = [400e6, 1e9], Verbose = True, Show_Plot = True):
		'''produces fit to cable loss in the functional form:
		term1 + term2 + term3 = a * sqrt(f) + b * f + c
		term1 is the sum of inner and outer coaxial cable conductor losses
		term2 is due to coaxial cable dielectric loss
		term3 is a constant fudge factor
		The loss evaluates to units of dB.

		stores the  fit as dictionary
		(a,b,c,run,range_start,range_stop)= self.metadata.Cable_Calibration['One_Way_40mk']

		Two used this function load transmission for complete cable loop only (not amps or attens).
		Then call this function on that transmission data. This funciton creats the tuple (a,b,c,run,range_start,range_stop) in 
		metadata, where run is the name of the calibration run and range_start/stop is the frequency range over which the
		calibration is calculated.

		Create a function from a,b,c and it to the effect of attenuators on the input side of the cable loop.

		set freq_range = None to use full freq range	
		'''
		
		fit_cable_loss(self.metadata, self.loop, key_name, freq_range, Verbose, Show_Plot)
		
		self._Cable_Calibration = self.metadata.Cable_Calibration 
		# The _Cable_Calibration variable is used to insert cable calibration 
		# data in the metadata at the time when Sweep_Array either by load_scandata 
		# or load_hf5_2... Consider removing in the future 



	def nonlinear_fit(self, Fit_Method = 'Multiple', Verbose = True, Show_Plot = True, Save_Fig = False, Compute_Chi2 = False, Indexing = (None,None,None)):
		'''
		The indexing keyword allows for selection of the power sweep to be fit. 
		If P is the list of powers then Indexing = (Start,Stop,Step) is using only, P[Start,Stop, Step]
		'''

		fit, fig, ax = nonlinear_fit(self.metadata, self.loop, self.Sweep_Array, Fit_Method, Verbose, Show_Plot, Save_Fig, Compute_Chi2, Indexing)
		return fit, fig, ax



	# def _nonlinear_formulae(self, parameter_dict, model = 2): Not Needed


	def generate_nonlinear_data(self,  Show_Plot = True, Phase_Noise_Variance = None, Amplitude_Noise_Variance = None, Like = None, Save_Fig = False,
		curve_parameter_dict = {'f_0':700e6, 'Qtl':300e3, 'Qc':80e3, 'eta':1e-1, 'delta':1e-6, 'Zfl':30, 'Zres':50, 'phi31': np.pi/2.03, 'phiV1':np.pi/10, 'V30V30':0.01},
		sweep_parameter_dict = {'Run': 'F1', 'Pprobe_dBm_Start' :-65.0,'Pprobe_dBm_Stop': -25.0, 'Pprobe_Num_Points':10, 'numBW':40,'num': 2000, 'Up_or_Down': 'Up', 'Freq_Spacing':'Linear'}):
		'''Creates and Loads Nonlinear Data
		eta -- Q nonlinearity
		delta --  freq nonlinearity	
		V30V30 -- V^2 normalization for nonlinearity

		If another KAM.sweep object is supplied in "Like" keyword, then its metadata will copied
		'''

		#delete previous metadata object and replace with new one
		del(self.metadata)
		self.metadata = metadata()

		#delete previous loop object and replace with new one  
		del(self.loop)
		self.loop = loop()

		fig, ax, self.sweep_data_columns_list, self.sweep_data_columns, self.Sweep_Array = generate_nonlinear_data(self.metadata, Show_Plot, Phase_Noise_Variance, Amplitude_Noise_Variance, Like, Save_Fig, curve_parameter_dict, sweep_parameter_dict)
		
		# select index  default_index now that new Sweep_Array is defined
		default_index = 0
		self.pick_loop(default_index)		

		return fig, ax

	# def _save_fig_dec(self, fig, name, Use_Date = False, Make_PGF = True): Not Needed
