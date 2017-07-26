def complete_fit(self, Use_Mask = True, Verbose = False , Show_Plot = False, Save_Fig = False, Sample_Size = 100, Use_Loop_Data = False):
	'''
	Sample_Size is the number of points used to extablish \sigma^2 for the gaussian noise model

	if Use_Loop_Data = True then values of Q, Qc, fr, phi are for initial guess are taken from curret loop object. If false, values come from self.Sweep_Array
	'''


	if self.loop.index == None:
		print 'Loop index is not specified. please pick_loop... Aborting'
		return

	Fit_Method = 'Multiple'
	if isinstance(Fit_Method,str): #Allow for single string input for Fit_Method
	   Fit_Method={Fit_Method}




	k = constants.value('Boltzmann constant') #unit is [J/k]
	BW = self.metadata.IFBW #unit is [Hz]	 
	# SC = self.metadata.System_Calibration # contains Noise powers, gains and P1dB of readout devices
	# CC = self.metadata.Cable_Calibration # cable loss fit coefficients
	R = 50 #system impedance

	#
	#
	# Implement Gain decompression on S21!
	#
	#
	if Use_Mask:
		F = ma.array(self.Sweep_Array[self.loop.index]['Frequencies'],mask = self.Sweep_Array[self.loop.index]['Mask'])
		F = F.compressed()
		S21 = ma.array(self.Sweep_Array[self.loop.index]['S21'],mask = self.Sweep_Array[self.loop.index]['Mask'])
		S21 = S21.compressed()
	else:
		F = self.Sweep_Array[self.loop.index]['Frequencies']
		S21 = self.Sweep_Array[self.loop.index]['S21']

	P_NA_out_dB = self.Sweep_Array[self.loop.index]['Pinput_dB'] #'out' as in our of NA, change of reference point
	P_NA_out_V2 = .001 * np.power(10,P_NA_out_dB/10) *2 *R 
	P_NA_in_V2 = np.square(np.abs(S21)) * P_NA_out_V2


	#Get chain and Noise temperatures for each element of readout chan and  at each frequency
	g_s , Tn_m_s ,Tn_p_s = self._construct_readout_chain(F)


	
	sigma_squared_m = np.zeros_like(F)
	sigma_squared_p = np.zeros_like(F)
	n = len(g_s)
	
	for i in xrange(n):	
		s2_m = 4*k*Tn_m_s[i]*R*BW # This sigma for the particular stage of the readout chain
		s2_p = 4*k*Tn_p_s[i]*R*BW
		#we assume s2_p * 4 * P_NA_in_V2 = s2_m ,  s2_p measured in radian^2
		sigma_squared_m = sigma_squared_m + s2_m*np.prod(g_s[i:], axis = 0) #rememebr g is a list of np vectors
		sigma_squared_p = sigma_squared_p + s2_p*np.square(np.prod(g_s[i:], axis = 0))/(4*P_NA_in_V2) #rememeber P_NA_in_V2 is a function of S21, see above definition




	if Use_Loop_Data == False:
		#a_0,b_0  = self.Sweep_Array[self.loop.index]['a'], self.Sweep_Array[self.loop.index]['b']
		R_0        = self.Sweep_Array[self.loop.index]['R']
		theta_0    = self.Sweep_Array[self.loop.index]['Theta']
		tau_0    = self.metadata.Electrical_Delay
		Q_0      = self.Sweep_Array[self.loop.index]['Q']
		Qc_0     = self.Sweep_Array[self.loop.index]['Qc']
		fr_0     = self.Sweep_Array[self.loop.index]['Fr'] 
		phi_0    = self.Sweep_Array[self.loop.index]['Phi']#(self.Sweep_Array[self.loop.index]['Phi'] * np.pi/180) + 0*np.pi
		
	else:
		#a_0,b_0  = self.loop.a, self.loop.b
		R_0      = self.loop.R
		theta_0  = self.loop.theta
		tau_0    = self.metadata.Electrical_Delay
		Q_0      = self.loop.Q 
		Qc_0     = self.loop.Qc 
		fr_0     = self.loop.fr 
		phi_0    = self.loop.phi# (self.loop.phi * np.pi/180) + 0*np.pi
		



	#p0 is the initial guess
	#p0 = np.array([a_0,b_0,tau_0,Q_0, Qc_0, fr_0, phi_0])
	p0 = np.array([R_0, theta_0,tau_0,Q_0, Qc_0, fr_0, phi_0])

	def obj(x,s21, sigma_squared_m,sigma_squared_p ,freq):# phase / magnitude fit
		# a,b,tau,Q, Qc, fr, phi= x
		# s21_fit  = norm * np.exp(np.complex(0.,np.angle(np.complex(a,b)))) * np.exp(np.complex(0,-2*np.pi*tau)*freq) * (1 - (Q/Qc)*np.exp(np.complex(0,phi)) / (1 + np.complex(0,2*Q)*(freq-fr)/fr ) )
		R,theta,tau,Q, Qc, fr, phi= x
		s21_fit  =  R * np.exp(np.complex(0.,theta)) * np.exp(np.complex(0,-2*np.pi*tau)*freq) * (1 - (Q/Qc)*np.exp(np.complex(0,phi)) / (1 + np.complex(0,2*Q)*(freq-fr)/fr ) )

		

		# diff = s21 - s21_fit
		# frac = (diff*diff.conj()).real/sigma_squared_m
		# #frac = (np.square(diff.real)/sigma_squared_m) + (np.square(diff.imag)/sigma_squared_m)
		frac = np.square(np.abs(s21) -  np.abs(s21_fit))*P_NA_out_V2/sigma_squared_m  + np.square(np.angle(s21/s21_fit))/sigma_squared_p  #(e^ia)/(e^ib) = e^i(a-b)
		N = freq.shape[0]*1.0 - x.shape[0]
		return  frac.sum()/N

	# Dont use masked data to sample points for Gaussian variance determination.
	if Use_Mask:
		S21_Sample = self.Sweep_Array[self.loop.index]['S21']
	else:
		S21_Sample = S21

	sigma_squared = 0 
	for i in xrange(Sample_Size):
		sigma_squared = sigma_squared + np.square(np.abs(S21_Sample[i] - S21_Sample[i+1]))
	sigma_squared = sigma_squared/(2.0*Sample_Size)
	
	def obj_s(x,s21, sigma_squared ,freq): # gaussian fit
		# a,b,tau,Q, Qc, fr, phi= x
		# s21_fit  = norm * np.exp(np.complex(0.,np.angle(np.complex(a,b)))) * np.exp(np.complex(0,-2*np.pi*tau)*freq) * (1 - (Q/Qc)*np.exp(np.complex(0,phi)) / (1 + np.complex(0,2*Q)*(freq-fr)/fr ) )
		R,theta,tau,Q, Qc, fr, phi= x
		s21_fit  =  R * np.exp(np.complex(0.,theta)) * np.exp(np.complex(0,-2*np.pi*tau)*freq) * (1 - (Q/Qc)*np.exp(np.complex(0,phi)) / (1 + np.complex(0,2*Q)*(freq-fr)/fr ) )


		# diff = s21 - s21_fit
		# frac = (diff*diff.conj()).real/sigma_squared_m
		# #frac = (np.square(diff.real)/sigma_squared_m) + (np.square(diff.imag)/sigma_squared_m)
		frac = np.square(np.abs(s21_fit-s21))/sigma_squared
		N = freq.shape[0]*1.0 - x.shape[0]
		return  frac.sum()/N
	
	
	#Each fit method is saved as a lambda function in a dictionary called fit_func
	fit_func = {}
	fit_func['cPowell'] = lambda : minimize(obj, p0, args=(S21,sigma_squared_m,sigma_squared_p ,F), method='Powell', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=1e-15, callback=None, options={'disp':False})
	fit_func['sPowell'] = lambda : minimize(obj_s, p0, args=(S21,sigma_squared ,F), method='Powell', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=1e-15, callback=None, options={'disp':False})

	#fit_func['Nelder-Mead']  = lambda : minimize(obj, p0, args=(S21,sigma_squared ,F), method='Nelder-Mead', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=1e-15, callback=None, options={'disp':False, 'xtol' : 1e-6,'maxfev':1000})
	#fit_func['Newton-CG'] = lambda : minimize(obj, p0, args=(z_theta_c,f_c), method='Newton-CG', jac=jac, hess=hess, hessp=None, bounds=None, constraints=(),tol=1e-15, callback=None, options={'maxiter' : 50,'xtol': 1e-4,'disp':False})

	fit = {}
	if isinstance(Fit_Method,set):      #All string inputs for Fit_Method were changed to sets at the begining of phase_fit
	   if Fit_Method == {'Multiple'}:
	      for method in fit_func.keys():
	         fit[method] = fit_func[method]() # Execute the fit lambda function
	   else:
	      for method in Fit_Method:
	         if method not in fit_func.keys():
	            print("Unrecognized fit method. Aborting fit. \n\t Must choose one of {0} or 'Multiple'".format(fit_func.keys()))
	            return
	         else:   
	            fit[method] = fit_func[method]()
	else:
	   print("Unrecognized fit method data type. Aborting fit. \n\t Please specify using a string or a set of strings from one of {0} or 'Multiple'".format(fit_func.keys()))
	   return	         	   
	
	if Verbose:        				
		for method in fit.keys():
			print('\n{0} Minimzation Result:\n{1}\n'.format(method,fit[method]))

	
	
	# bestfit = list(fit)[0]
	# lowest = fit[bestfit].fun
	# for key in fit.keys(): 
	# 	if fit[key].fun < lowest:
	# 		lowest = fit[key].fun
	# 		bestfit = key

	cfit = 'cPowell' # phase / mag chi squared 
	#ca, cb, ctau = fit[cfit].x[0], fit[cfit].x[1], fit[cfit].x[2]
	self.loop.cR, self.loop.ctheta, ctau = cR, ctheta, ctau = fit[cfit].x[0], fit[cfit].x[1], fit[cfit].x[2]
	self.loop.cQ = cQ = fit[cfit].x[3]	
	self.loop.cQc = cQc = fit[cfit].x[4]	
	self.loop.cQi = cQi = 1.0/ ((1./self.loop.cQ ) - (1./self.loop.cQc )) 
	self.loop.cfr = cfr = fit[cfit].x[5]	
	self.loop.cphi = cphi =  fit[cfit].x[6]	
	self.loop.cchisquare = fit[cfit].fun
	self.loop.cphase_fit_success = fit[cfit].success

	sfit = 'sPowell' #gaussian chi squared
	#sa, sb, stau = fit[sfit].x[0], fit[sfit].x[1], fit[sfit].x[2]
	self.loop.sR, self.loop.stheta, stau  =sR, stheta, stau= fit[sfit].x[0], fit[sfit].x[1], fit[sfit].x[2]
	self.loop.sQ = sQ = fit[sfit].x[3]	
	self.loop.sQc = sQc = fit[sfit].x[4]	
	self.loop.sQi = sQi = 1.0/ ((1./self.loop.sQ ) - (1./self.loop.sQc )) 
	self.loop.sfr = sfr = fit[sfit].x[5]	
	self.loop.sphi = sphi =  fit[sfit].x[6]	
	self.loop.schisquare = fit[sfit].fun
	self.loop.sphase_fit_success = fit[sfit].success

	fit['sigma_squared_m'] = sigma_squared_m
	fit['sigma_squared_p'] = sigma_squared_p
	fit['sigma_squared'] = sigma_squared


	if  Show_Plot:
		ax_dict = {}
		fig = plt.figure( figsize=(6.5, 6.5), dpi=100)
		fig_dict = {fig : ax_dict}
		ax = fig.add_subplot(111,aspect='equal')
		lines = []
		s21_concurrent_c = cR * np.exp(np.complex(0.,ctheta)) * np.exp(np.complex(0,-2*np.pi*ctau)*F) * (1 - (cQ/cQc)*np.exp(np.complex(0,cphi)) / ( 1 + np.complex(1, 2*cQ)*(F-cfr)/cfr  ))
		# s21_concurrent_c = norm * np.exp(np.complex(0.,np.angle(np.complex(ca,cb)))) * np.exp(np.complex(0,-2*np.pi*ctau)*F) * (1 - (cQ/cQc)*np.exp(np.complex(0,cphi)) / ( 1 + np.complex(1, 2*cQ)*(F-cfr)/cfr  ))
		lines.append(ax.plot(s21_concurrent_c.real,s21_concurrent_c.imag, markersize  = 3, linestyle = 'None',color = 'g', marker = 'o', markerfacecolor = 'g', markeredgecolor = 'g',  label = r'Concurrent Fit -  $\sigma_{V\theta}$')[0])

		s21_concurrent_s = sR * np.exp(np.complex(0.,stheta)) * np.exp(np.complex(0,-2*np.pi*stau)*F) * (1 - (sQ/sQc)*np.exp(np.complex(0,sphi)) / ( 1 + np.complex(1, 2*sQ)*(F-sfr)/sfr  ))
		#s21_concurrent_s = norm * np.exp(np.complex(0.,np.angle(np.complex(sa,sb)))) * np.exp(np.complex(0,-2*np.pi*stau)*F) * (1 - (sQ/sQc)*np.exp(np.complex(0,sphi)) / ( 1 + np.complex(1, 2*sQ)*(F-sfr)/sfr  ))
		lines.append(ax.plot(s21_concurrent_s.real,s21_concurrent_s.imag,markersize  = 3, color = 'm',linestyle = 'None', marker = 'o', markerfacecolor = 'm', markeredgecolor = 'm',  label = r'Concurrent Fit -  $\sigma_{G}$')[0])
		lines.append(ax.plot(s21_concurrent_s[0:Sample_Size:].real,s21_concurrent_s[0:Sample_Size:].imag,'m+', label = r'_Concurrent Fit -  $\sigma_{G}$')[0])

		lines.append(ax.plot(S21.real,S21.imag,markersize  = 3,color = 'b' ,marker = 'o',  linestyle = 'None',markerfacecolor = 'b', markeredgecolor = 'b', label = r'Raw Data - $S_{21}$')[0])


		s21_stepwise  =  R_0 * np.exp(np.complex(0.,theta_0)) * np.exp(np.complex(0,-2*np.pi*tau_0)*F) * (1 - (Q_0/Qc_0)*np.exp(np.complex(0,phi_0)) /( 1 + np.complex(1, 2*Q_0)*(F-fr_0)/fr_0  ))
		#s21_stepwise  = norm * np.exp(np.complex(0.,np.angle(np.complex(a_0,b_0)))) * np.exp(np.complex(0,-2*np.pi*tau_0)*F) * (1 - (Q_0/Qc_0)*np.exp(np.complex(0,phi_0)) /( 1 + np.complex(1, 2*Q_0)*(F-fr_0)/fr_0  ))
		lines.append(ax.plot(s21_stepwise.real,s21_stepwise.imag,markersize  = 3, color = 'r', linestyle = 'None',marker = 'o', markerfacecolor = 'r', markeredgecolor = 'r', label = r'Stepwise Fit - $\hat{S}_{21}$')[0])
		ax_dict.update({ax:lines})


		ax.set_xlabel(r'$\Re[S_{21}(f)]$')
		ax.set_ylabel(r'$\Im[S_{21}(f)]$')
		ax.yaxis.labelpad = -2
		ax.legend(loc = 'upper center', fontsize=5, bbox_to_anchor=(0.5, -0.1), ncol=2,scatterpoints =1, numpoints = 1, labelspacing = .02)
		#ax.legend(loc = 'best', fontsize=9,scatterpoints =1, numpoints = 1, labelspacing = .02) 

		plot_dict = fig_dict
		plt.show()	
	else:
		plot_dict =  None
	# if  Show_Plot:	
	# 	plt.show()

	if Save_Fig == True:
		self._save_fig_dec(fig,'Concurrent_Fit_Index_{0}'.format(self.loop.index))
	


	return fit, plot_dict
