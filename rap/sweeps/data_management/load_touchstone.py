def load_touchstone(metadata, filename, pick_loop = True):
	''' The function loads S21 and Freq from  Sonnet .s2p or .s3p files into the Sweep_Array structured np array
	All Sij are extracted, but only  S21 is saved into Sweep_Array. Future editions of this code might  find need 
	to load other Sij becuase S21.

	The function only loads one transmission array (S21).  pick_loop = True immediatly selectes this loop as the 
	current loop.
	'''

	import tempfile
	import io

	#delete previous metadata object
	del(metadata)
	metadata = metadata()

	dt_s2p = [('Freq', np.float64), ('S11r', np.float64), ('S11i', np.float64), ('S12r', np.float64), ('S12i', np.float64), 
									('S21r', np.float64), ('S21i', np.float64), ('S22r', np.float64), ('S22i', np.float64)]
	
	dt_s3p = [('Freq', np.float64), ('S11r', np.float64), ('S11i', np.float64), ('S12r', np.float64), ('S12i', np.float64), ('S13r', np.float64), ('S13i', np.float64),
									('S21r', np.float64), ('S21i', np.float64), ('S22r', np.float64), ('S22i', np.float64), ('S23r', np.float64), ('S23i', np.float64),
									('S31r', np.float64), ('S31i', np.float64), ('S32r', np.float64), ('S32i', np.float64), ('S33r', np.float64), ('S33i', np.float64)] 


	with tempfile.TemporaryFile() as tmp:
		with io.open(filename, mode='r') as f:
			# The following while loop copies the .sNp file into a temp file, which is destroyed when closed,
			# such that the tmp file is formated in the way np.loadtxt can read the data.
			indented = False
			prev_line = ''
			m = 1. # for frequency base conversion
			while 1: 
				line  = f.readline().replace('\n','')

				pos = f.tell()
				if line == '': # End of file is reached
					break
				elif line.startswith('! Data File Written:'): # Save as Metadata
					metadata.Time_Created = str(line.split('! Data File Written:')[1].strip())
					tmp.write(line + '\n')
				elif line.startswith('! From Project:') | line.startswith('! From Emgraph Data:'): # Save as Metadata
					metadata.Run = str(line.split(':')[1].strip())
					#self.metadata.IS_Sonnet_Simulation = True
					tmp.write(line + '\n')
				elif line[0] == '#':
					line  = line.replace('#','!#')
					if line.find('GHZ') >=-1:
						m = 1.0e9
					freq_convert = lambda s: s*m #Convert to Hertz
					tmp.write(line + '\n')	
				
				elif line[0] == ' ': # in S matrix definition block
					prev_line = prev_line + ' ' + line.strip() + ' '
					next_line = f.readline()
					# if next line is NOT indented date, then S matrix definition block is finished 
					# and we write it to tmp on a single line.
					# for .s2p files the S matrix is fully defined on one line of f
					# for .s3p files, the S matrix is defined in three lines. second two are indented.
					
					# if not ((next_line[0] == '') | (next_line[0] == ' ')): # Changing this line to be consistent with line below...
					if not ((next_line == '') or (next_line[0] == ' ')): 
						tmp.write(prev_line)
						tmp.write('\n')
						prev_line = ''
					f.seek(pos,0)
	
				elif line[0] == '!':
					tmp.write(line + '\n')

				else:
					tmp.write(line)
					next_line = f.readline()
					# add \n to line if it does not begin a S matrix definition block
					# if not ((next_line[0] == '') | (next_line[0] == ' ')): # Changed on 7/11/17 bc Nick was have problems reading in .s2p files
					if not ((next_line == '') or (next_line[0] == ' ')):
						tmp.write('\n')
					f.seek(pos,0)

		tmp.seek(0,0)
		if filename.endswith('.s2p'):
			dt = dt_s2p
		elif filename.endswith('.s3p'):
			dt = dt_s3p	
		Touchstone_Data = np.loadtxt(tmp, dtype=dt, comments='!', delimiter=None, converters=None, skiprows=0, usecols=None, unpack=False, ndmin=0)
	
	tpoints = 0
	self._define_sweep_data_columns(Touchstone_Data.size, tpoints)
	j = np.complex(0,1)

	self.Sweep_Array = np.zeros(1, dtype = self.sweep_data_columns)
	
	self._define_sweep_array(0, Fstart = freq_convert(Touchstone_Data['Freq'].min()), #Hz
								Fstop = freq_convert(Touchstone_Data['Freq'].max()), #Hz
								S21 = Touchstone_Data['S21r']+j*Touchstone_Data['S21i'],
								Frequencies = freq_convert(Touchstone_Data['Freq']), #Hz
								#Pinput_dB = 0,
								Is_Valid = True,
								#Mask = False, needs to be an array of lengh of S21
								Chi_Squared = 0,
								)


	metadata.Data_Source = filename
	#self.metadata.Min_Freq_Resolution = np.abs(Touchstone_Data['Freq'][:-1]-Touchstone_Data['Freq'][1:]).min()
	metadata.Min_Freq_Resolution = np.abs(Touchstone_Data['Freq'][0] - Touchstone_Data['Freq'][-1])/Touchstone_Data['Freq'].size #use average freq resolution

	if pick_loop == True: #since there is only one loop in Sweep_Array, we might as well pick it as the current loop
		self.pick_loop(0)
		#self.normalize_loop()
