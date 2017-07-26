def fill_sweep_array(self, Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False, Complete_Fit = True , Remove_Gain_Compression = True, Verbose = True):
	

	if Compute_Preadout == True:
		needed = ('Atten_NA_Output', 'Atten_At_4K','Cable_Calibration')

		for quantities  in needed:				
			if  self.metadata.__dict__[quantities] == None:
				if Verbose == True:
					print('{0} metadate missing. Unable to compute Preadout. Setting to 0.'.format(quantities))
				Compute_Preadout = False

			Atten_NA_Output = self.metadata.Atten_NA_Output
			Atten_At_4K = self.metadata.Atten_At_4K
			Cable_Calibration_Key = 'One_Way_40mK'
			k = self.metadata.Cable_Calibration[Cable_Calibration_Key]

		if Fit_Resonances == False:
			if Verbose == True:
				print('Resonance fit not selected. Computation of Preadout_dB requires knowledge of resonance frequency and may not work.')


		if Compute_Preadout == True:
			Preadout = lambda f: k[0]*np.sqrt(f)+k[1]*f+k[2] - Atten_NA_Output - Atten_At_4K

	if Add_Temperatures == True:
		if self.metadata.Num_Temperatures < 1:
			Temperature_Calibration = self.metadata.Temperature_Calibration
			if (self.metadata.Fridge_Base_Temp != None) & (max(self.Sweep_Array['Heater_Voltage']) == min(self.Sweep_Array['Heater_Voltage'])): #& (self.Sweep_Array.size == 1):
				#This is usually the case of a survey or power sweep: done at base temp with no Heater power
				self.Sweep_Array['Temperature'][:] = self.metadata.Fridge_Base_Temp
				print('Setting Tempreature to metadata.Fridge_Base_Temp value.')
				Add_Temperatures = False

			elif type(Temperature_Calibration) == list: 
				Temperature_Calibration = np.array(Temperature_Calibration)
				# Temperature_Calibration[:,0] is heater voltages
				# Temperature_Calibration[:,1] is temperatures voltages
				
				# becasue ScanData heater voltages are read in as numbers like 0.24999999 and 0.2500001 instread of 0.25
				# as included in the Temperature_Calibration list/array, use this 'tol' to associate closest ScanData 
				# heater voltage to voltage in Temperature_Calibration list/array.
				tol =  0.0005 
			
			else:
				if Verbose == True:
					print('Temperature_Calibration metadata is not found or not of the correct type. Unable to add temperatures.')
				Add_Temperatures = False
		else:
			tol = None
			pass


		
	num_records = self.Sweep_Array.size
	for index in xrange(num_records): 
		if Verbose == True:
			sys.stdout.write('\r {0} of {1} '.format(index+1, num_records))
			sys.stdout.flush()

		#set current loop
		self.pick_loop(index)

		if Fit_Resonances == True:
			if Remove_Gain_Compression:
				# Remove Gain Compression
				self.decompress_gain(Compression_Calibration_Index = -1, Show_Plot = False, Verbose = False)

			if self.loop.z.size > 5000:
				self.trim_loop(N = 10, Verbose = False)

			# Normalize Loop
			#self.normalize_loop() 

			# Remove Cable Delay
			self.remove_cable_delay(Show_Plot = False, Verbose = False)	# should do nothing if a delay is defined in metadata

			# Fit loop to circle
			self.circle_fit(Show_Plot = False)
			if self.loop.circle_fit_exit_code != 0:
				self._define_sweep_array(index, Is_Valid = False)

			# Fit resonance parameters
			self.phase_fit(Fit_Method = 'Multiple',Verbose = False, Show_Plot = False)
			

			self._define_sweep_array(index, Q = self.loop.Q,
											Qc = self.loop.Qc,
											Fr = self.loop.fr,
											Mask = self.loop.phase_fit_mask,
											Chi_Squared = self.loop.chisquare,
											R = self.loop.R,
											r = self.loop.r,
											a = self.loop.a,
											b = self.loop.b,
											#Normalization  = self.loop.normalization,
											Theta = self.loop.theta,
											Phi = self.loop.phi)

			if Complete_Fit:
				self.complete_fit(Use_Mask = True, Verbose = False , Show_Plot = False, Save_Fig = False, Sample_Size = 100, Use_Loop_Data = True)
				self._define_sweep_array(index, cQ = self.loop.cQ,
												cQc = self.loop.cQc,
												cFr = self.loop.cfr,
												cPhi = self.loop.cphi,
												cTheta = self.loop.ctheta,
												cR = self.loop.cR,
												cChi_Squared = self.loop.cchisquare,
												cIs_Valid = self.loop.cphase_fit_success if self.Sweep_Array['Is_Valid'][index] else self.Sweep_Array['Is_Valid'][index],

												sQ = self.loop.sQ,
												sQc = self.loop.sQc,
												sFr = self.loop.sfr,
												sPhi = self.loop.sphi,
												sTheta = self.loop.stheta,
												sR = self.loop.sR,
												sChi_Squared = self.loop.schisquare,
												sIs_Valid = self.loop.sphase_fit_success if self.Sweep_Array['Is_Valid'][index] else self.Sweep_Array['Is_Valid'][index]
												)

											


			# Only execute if phase_fit_success is False to avoid setting Is_Valid true when it was previously set fulse for a different reason, e.g bad Temp data
			if self.loop.phase_fit_success == False: 
				self._define_sweep_array(index, Is_Valid = False)

		if Compute_Preadout == True:
			if self.loop.fr != None:
				self._define_sweep_array(index, Preadout_dB = self.Sweep_Array['Pinput_dB'][index] + Preadout(self.loop.fr))
			elif np.abs(self.loop.freq[-1]-self.loop.freq[0]) > 1e9:
				if Verbose == True:
					print('Sweep bandwidth is {0} Hz. Sweep looks more like a survey. Preadout_dB is meaningless for a survey. Aborting Preadout computation... '.format(np.abs(self.loop.freq[-1]-self.loop.freq[0])))
				
			else:
				if Verbose == True:
					print('No resonance frquency (fr) on record for selected resonance. Estimating fr using sweep minimum.')
				fr = np.extract(np.abs(self.loop.z).min() == np.abs(self.loop.z),self.loop.freq)[0]
				self._define_sweep_array(index, Preadout_dB = self.Sweep_Array['Pinput_dB'][index] + fr)

		if Add_Temperatures == True:
			if self.metadata.Num_Temperatures < 1:
				condition = (self.Sweep_Array['Heater_Voltage'][index] + tol > Temperature_Calibration[:,0]) & (self.Sweep_Array['Heater_Voltage'][index] - tol < Temperature_Calibration[:,0])
				if condition.sum() >= 1:

					self.Sweep_Array['Temperature'][index] = Temperature_Calibration[condition,1][0] # <-- Needs to be updated so that duplicate voltages are handled correctly
				else:
					if Verbose == True:
						print('Unable to match unique temperature to heater voltage value for Sweep_Array[{0}]. {1} matches found.'.format(index,condition.sum() ))
			else:
				self._define_sweep_array(index, Temperature = 	self.Sweep_Array['Temperature_Readings'][index].mean()) 		
		# Clear out loop
		del(self.loop)
		self.loop = loop()
	if Verbose == True:
		print('\nSweep Array filled.')# Options selected Fit_Resonances = {0}, Compute_Preadout = {1}, Add_Temperatures = {2}'.format( Fit_Resonances,Compute_Preadout,Add_Temperatures))
