


### external imports
import scipy.io #for loading .mat file
import urllib2
import numpy as np
import os

def _read_scandata_from_file(filename_or_path):
		
	mat = scipy.io.loadmat(filename_or_path)

	return mat, filename_or_path


def _download_data(URL, **auth):
	''' Authenticats to URL containing data.
	Copies the .mat file licated at URL to a local file in local directory.
	.mat file is a Scan_Data matlab structure.
	returns numpy data structure contauning .mat file.
	deletes local file.'''



	if 'username' in auth.keys() and 'password' in auth.keys():
		passman = urllib2.HTTPPasswordMgrWithDefaultRealm()
		# this creates a password manager

		username = auth['username']
		password = auth['password']

		passman.add_password(None, URL, username, password)
		# because we have put None at the start it will always
		# use this username/password combination for  urls
		# for which `URL` is a super-url

		authhandler = urllib2.HTTPBasicAuthHandler(passman)
		# create the AuthHandler

		opener = urllib2.build_opener(authhandler)

		urllib2.install_opener(opener)
		# All calls to urllib2.urlopen will now use our handler
		# Make sure not to include the protocol in with the URL, or
		# HTTPPasswordMgrWithDefaultRealm will be very confused.
		# You must (of course) use it when fetching the page though.

	pagehandle = urllib2.urlopen(URL)
	# authentication is now handled automatically for us

	#import tempfile # Attempt to download data into a temp file
	#f = tempfile.NamedTemporaryFile(delete=False)
	#f.write(pagehandle.read())
	#f.close()
	#mat = scipy.io.loadmat(f.name)

	output = open('test.mat','wb')
	print('Download Initiated...')
	output.write(pagehandle.read())
	print('Download Completed...')
	output.close()
	#global mat
	mat = scipy.io.loadmat('test.mat')
	

	#this id how to tell what variables are stored in test.mat
	#print scipy.io.whosmat('test.mat')



	#html = pagehandle.read()
	#pagehandle.close()

	#soup = BeautifulSoup(html)
	#soup.contents
	os.remove('test.mat')
	return mat, URL
	
def _extract_type(obj, return_type = None, field = None):
	'''scanandata object, obj, has a lot of single element arrays of arrays. this function gets the element.
	e.g scandata may have [[[ele]]] instead of callling ele = scandata[0][0][0], use this function to get ele.
	if ele is another structured numpy array with field name 'myfield', using keyword field = 'myfield' will get
	the data at field.
	the function will cast ele to be in the data type return_typer. e.g. return_type = 'str' returns a string. 
	If return_type is None, ele is returned as whatever type it was saved as in [[[ele]]] '''
	
	def cast(_obj):
		if (return_type is not None) & (_obj is not None) : #if (return_type != None) & (_obj != None) :
			_obj = return_type(_obj)
			#pass#exec("_obj = {0}(_obj)".format(return_type))
		return _obj

	def itemloop(_obj):
		while True:
			try:
				_obj = _obj.item()
			except:
				return cast(_obj)
		return cast(_obj)


	if field == None:
		obj = itemloop(obj)

	else:
		while obj.dtype == np.dtype('O'):
			obj = obj.item()
			
		if isinstance(obj.item(), unicode): 
			obj = None
			print('Expected dictionary containing field named {0} is not found. Returning None'.format(field))				
		else: #if the object does not simply contain a string, e.g  [u'InP #2'], do this
			try:
				obj = obj[field]
			except:
				obj = None
				print('Field named {0} is not found. Cannot extract. Returning None'.format(field))	
		# try:
		# 	obj = obj[field]
		# except:
		# 	obj = None
		# 	print('Field named {0} is not found. Returning None'.format(field))
		obj = itemloop(obj)
	return obj

def _define_sweep_data_columns(metadata, fsteps, tpoints):
	''' Create the sweep_data_columns_list which is used to define the dtype of the  Sweep_Array'''
	metadata.Fsteps = fsteps
	metadata.Num_Temperatures  = tpoints

	if tpoints < 1: # we dont want a shape = (0,) array. We want at least (1,)
		tpoints = 1

	sweep_data_columns_list = [
		("Fstart"         			, np.float64), # in Hz
		("Fstop"          			, np.float64), # in Hz
		("Heater_Voltage" 			, np.float64), # in Volts
		("Pinput_dB"      			, np.float64), # in dB
		("Preadout_dB"     			, np.float64), # in dB  - The power at the input of the resonator, not inside the resonator
		("Thermometer_Voltage_Bias"	, np.float64), # in Volts
		("Temperature_Readings"    	, np.float64,(tpoints,)), # in Kelvin
		("Temperature"		    	, np.float64), # in Kelvin
		("S21"            			, np.complex128, (fsteps,)), # in complex numbers, experimental values.
		("Frequencies"    			, np.float64,(fsteps,)), # in Hz
		("Q"						, np.float64),
		("Qc"						, np.float64),
		("Fr"						, np.float64), # in Hz
		("Is_Valid"					, np.bool),
		("Chi_Squared"              , np.float64),
		("Mask"						, np.bool,(fsteps,)), # array mask selecting data used in phase fit
		("R"						, np.float64), #outer loop radius
		("r"						, np.float64), # resonance loop radius	
		("a"						, np.float64),	
		("b"						, np.float64),
		#("Normalization"			, np.float64),
		("Theta"					, np.float64),
		("Phi"						, np.float64),
		("cQ"						, np.float64),
		("cQc"						, np.float64),
		("cFr"						, np.float64), # in Hz
		("cIs_Valid"				, np.bool),
		("cChi_Squared"             , np.float64),
		("cPhi"						, np.float64),
		("cTheta"					, np.float64),
		("cR"						, np.float64),
		("sQ"						, np.float64),
		("sQc"						, np.float64),
		("sFr"						, np.float64), # in Hz
		("sIs_Valid"				, np.bool),
		("sChi_Squared"             , np.float64),
		("sPhi"						, np.float64),
		("sTheta"					, np.float64),
		("sR"						, np.float64),

		#("S21_Processed"            , np.complex128, (fsteps,)), # Processed S21 used in phase fit 
		]
	sweep_data_columns = np.dtype(sweep_data_columns_list)
	
	return sweep_data_columns_list, sweep_data_columns


def _compute_noise_spectrum_length(measurement_metadata):
	'''
	measurement_metadata is a dict.
	Computes length of noise spectrum, whcih is needed for 
	construction of the Sweep_Array data type before measurement of noise.
	'''
	frequency_segmentation = measurement_metadata['Noise_Frequency_Segmentation'] 
	frequency_resolution_per_segment = measurement_metadata['Noise_Frequency_Resolution_Per_Segment'] 
	noise_spectrum_length = 0
	for i in xrange(len(frequency_segmentation)):
		last_segmentation = frequency_segmentation[i-1] if i > 0 else 0 
		noise_spectrum_length = noise_spectrum_length + (frequency_segmentation[i] - last_segmentation)/frequency_resolution_per_segment[i]
	#try to get rid of this by correcting KIDS-DAQ-75uW Measurement_Managers code	
	noise_spectrum_length = noise_spectrum_length if measurement_metadata.has_key('Is_Legacy_Gui_Data') and measurement_metadata['Is_Legacy_Gui_Data'] else noise_spectrum_length - 1

	if divmod(noise_spectrum_length,1)[1] != 0: 
		ValueError('Noise spectrum length is not an integer is not an integer')
	return int(noise_spectrum_length)


def _define_sweep_data_columns_legacy_gui(measurement_metadata, fsteps_syn = 1):
	'''
	in the legacy gui, we use 'Synthesizer_Scan_Num_Points' = fsteps instead of the usual 'NA_Scan_Num_Points'
	'''
	#measurement_metadata['NA_Scan_Num_Points'] = fsteps
	#measurement_metadata['Num_Temperatures']  = tpoints
	#self.measurement_metadata['Synthesizer_Scan_Num_Points']  = fsteps_syn


	if measurement_metadata['Measure_On_Res_Noise']:
		noise_spectrum_length = _compute_noise_spectrum_length(measurement_metadata)
	else:
		noise_spectrum_length = 1

	# off resonance noise vector will be the same length as the on resoance noise vectors, 
	# unless no off resonance noise it to be measured. In which case, the off resoance noise 
	# are of length 1
	if measurement_metadata['Measure_Off_Res_Noise']:
		noise_spectrum_length_off_res = noise_spectrum_length
	else:
		noise_spectrum_length_off_res = 1

	# if tpoints < 1: # we dont want a shape = (0,) array. We want at least (1,)
	# 	tpoints = 1

	# if fsteps_syn < 1:
	# 	fsteps_syn = 1

	if noise_spectrum_length < 1:
		noise_spectrum_length = 1



	sweep_data_columns_list = [
		#("Fstart"         			, np.float64), # in Hz
		#("Fstop"          			, np.float64), # in Hz
		#("Heater_Voltage" 			, np.float64), # in Volts
		("Pinput_dB"      			, np.float64), # in dB
		#("Preadout_dB"     			, np.float64), # in dB  - The power at the input of the resonator, not inside the resonator
		#("Thermometer_Voltage_Bias"	, np.float64), # in Volts
		#("Temperature_Readings"    	, np.float64,(tpoints,)), # in Kelvin
		("Temperature"		    	, np.float64), # in Kelvin
		#("S21"            			, np.complex128, (fsteps,)), # in complex numbers, experimental values.
		#("Frequencies"    			, np.float64,(fsteps,)), # in Hz
		("Is_Valid"					, np.bool),
		#("Aux_Voltage"				, np.float64), # in Volts
		#("Aux_Value"				, np.float64), # value of aux signal due to Aux_Voltage, e.g. B-field in gauss if Helmholtz coil is on aux channel
		("S21_Syn"            		, np.complex128, (fsteps_syn,)), # in complex numbers
		("Frequencies_Syn"    		, np.float64,(fsteps_syn,)), # in Hz
		#("Bkgd_Freq_Syn"    		, np.float64,(bkgd_loop_num_pts,)), # in Hz
		#("Bkgd_S21_Syn"            	, np.complex128, (bkgd_loop_num_pts,)), # in complex numbers
		("Q"						, np.float64),
		("Qc"						, np.float64),
		("Fr"						, np.float64), # in Hz - Res freq estimated from na scan
		("Power_Baseline_Noise"    	, np.float64), # dBm -  off resonnce power entering digitizer for noise spectrum measurement		
		("Noise_Freq_On_Res"    	, np.float64), # Hz - freq at which on res noise was taken
		("Noise_S21_On_Res"         , np.complex128), # The I/Q point at which on res noise was taken
		#("Noise_S21_Std_On_Res"     , np.complex128), # The standard deviation in the measurement of the I/Q point at which on res noise was taken
		("Noise_II_On_Res"			, np.float64,(noise_spectrum_length,)), # PSD in V^2/Hz
		("Noise_QQ_On_Res"			, np.float64,(noise_spectrum_length,)),
		("Noise_IQ_On_Res"			, np.complex128,(noise_spectrum_length,)),
		("Noise_Freq_Off_Res"    	, np.float64), # Hz - freq at which off res noise was taken			
		("Noise_II_Off_Res"			, np.float64,(noise_spectrum_length_off_res,)),
		("Noise_QQ_Off_Res"			, np.float64,(noise_spectrum_length_off_res,)),
		("Noise_IQ_Off_Res"			, np.complex128,(noise_spectrum_length_off_res,)),
		#("Noise_Atten_Mixer_Input"  , np.float32), # actual readings of attenuation value at mixer input (box chan 2) from Attenuator box
		#("Noise_Atten_Fridge_Input" , np.float32), # actual readings of attenuation value at input to fridge (box chan 1) from Attenuator box
		("Noise_S21_Off_Res"        , np.complex128), # The I/Q point at which off res noise was taken
		#("Noise_S21_Std_Off_Res"    , np.complex128), # The standard deviation in the measurement of the I/Q point at which off res noise was taken
		("Noise_Freq_Vector"    	, np.float64,(noise_spectrum_length,)), # in Hz, the  frequencies of the noise spectrum
		("Noise_Chan_Input_Atn"		, np.uint32), # attenuator box attenuation value setting for input side of fridge (Atn Box Chan 1)
		("Noise_Chan_Output_Atn"	, np.uint32), # attenuator box attenuation value setting for output side of fridge (Atn Box Chan2)
		#phase cancellation value for carrier cuppresion - np.float64,(fsteps_syn,)
		#ampl cancellation value for carrir suppression -  np.float64,(fsteps_syn,)
		("Scan_Timestamp"			, '|S12'),
		("Resonator_Index"			, np.uint32), #unique number assigned to each unique resonator
		("Resonator_Group"			, np.uint32,(3,)), #[groupnumber, res1_index, res2_index]
		("dZdf"						, np.float64),
		]

	sweep_data_columns = np.dtype(sweep_data_columns_list)
	
	return sweep_data_columns_list, sweep_data_columns

def _define_sweep_array(Sweep_Array, index,**field_names):
	#for field_name in self.sweep_data_columns.fields.keys():
	for field_name in field_names:
		Sweep_Array[field_name][index] = field_names[field_name]
	
					