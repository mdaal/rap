
### external imports
import scipy.io #for loading .mat file

def _read_scandata_from_file(filename_or_path):
		
	mat = scipy.io.loadmat(filename_or_path)

	return mat


def _download_data(self, URL):
	''' Authenticats to URL containing data.
	Copies the .mat file licated at URL to a local file in local directory.
	.mat file is a Scan_Data matlab structure.
	returns numpy data structure contauning .mat file.
	deletes local file.'''


	passman = urllib2.HTTPPasswordMgrWithDefaultRealm()
	# this creates a password manager
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
	self.data = mat
	self.metadata.Data_Source = URL
def _extract_type(self, obj, return_type = None, field = None):
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
				print('Field named {0} is not found. Returning None'.format(field))	
		# try:
		# 	obj = obj[field]
		# except:
		# 	obj = None
		# 	print('Field named {0} is not found. Returning None'.format(field))
		obj = itemloop(obj)
	return obj
def _define_sweep_data_columns(self, fsteps, tpoints, list_only = False):
	self.metadata.Fsteps = fsteps
	self.metadata.Num_Temperatures  = tpoints

	if tpoints < 1: # we dont want a shape = (0,) array. We want at least (1,)
		tpoints = 1

	self.sweep_data_columns_list = [
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
	if list_only == False:
		self.sweep_data_columns = np.dtype(self.sweep_data_columns_list)
def _define_sweep_array(self,index,**field_names):
	#for field_name in self.sweep_data_columns.fields.keys():
	for field_name in field_names:
		self.Sweep_Array[field_name][index] = field_names[field_name]
					