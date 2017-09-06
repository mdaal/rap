import numpy as np

def _construct_readout_chain(metadata, F, Include_NA = True, Include_4K_to_40mK = False):
	'''
	F is a frequency array.
	Constructs gain, Tn_m (T noise magnitude), and Tn_p (phase)  lists.
	Each element of list corresponds to a component, e.g. primary amp, cable 1, second amp, attenator,.
	The order of the list correspondes to the order of components in the readout chain. First element is the first component (e.g. the primary amp)
	Each element of the list is an array the same shape as F. Each element of the arrays is the gain, Tn_m (T noise magnitude), and Tn_p (phase) at that frequency.

	This method does not use self.loop. data. It only uses self.metadata

	The System_Calibration and Cable_Calibration data are input into metadate at the time of data library creating (in the file Create_Lbrary.py)

	'''
	SC = metadata.System_Calibration # contains Noise powers, gains and P1dB of readout devices
	CC = metadata.Cable_Calibration # cable loss fit coefficients

	# Chain is the string of readout cables and amplifiers/devices
	chain  = []
	
	if Include_4K_to_40mK:
		chain.append('4K_to_40mK')

	if metadata.LNA['LNA'] is not None:
		chain.append(metadata.LNA['LNA'])

	chain.append('300K_to_4K')

	if metadata.RTAmp_In_Use:
		chain.append(metadata.RTAmp) 

	
	chain.append('One_Way_300K')

	if (metadata.Atten_NA_Input is not None) and (metadata.Atten_NA_Input>0):
		chain.append('Atten_NA_Input')

	if Include_NA:
		chain.append('NA')


	passive_device_temp = {'4K_to_40mK': (4. +.04)/2, '300K_to_4K' : (290.+4.)/2, 'One_Way_300K': 290., 'Atten_NA_Input':290.}
	Tn_p_s = []
	Tn_m_s = []
	g_s = []
	for i in xrange(len(chain)):
		device = chain[i]

		if device in CC.keys():
			g = CC[device][0]*np.sqrt(F)+CC[device][1]*F+CC[device][2]
			g = np.power(10.0,g/10.0)
			g_s.append(g)
			Tn = ((1.0/g)-1)*passive_device_temp[device]
			Tn_p_s.append(Tn)
			Tn_m_s.append(Tn)
			continue

		if device in SC.keys():
			g = np.polynomial.chebyshev.chebval(F,SC[device]['g_fit'])
			g = np.power(10.0,g/10.0)
			g_s.append(g)
			Tn_p_s.append(np.polynomial.chebyshev.chebval(F,SC[device]['Tn_p_fit']))
			Tn_m_s.append(np.polynomial.chebyshev.chebval(F,SC[device]['Tn_m_fit']))
			continue

		if device is 'Atten_NA_Input':
			g =  -np.abs(metadata.Atten_NA_Input)*np.ones_like(F)
			g = np.power(10.0,g/10.0)
			g_s.append(g)
			Tn = ((1.0/g)-1)*passive_device_temp[device]
			Tn_p_s.append(Tn)
			Tn_m_s.append(Tn)
			continue

		# warn me if the component is missing from calibration data
		print('Component in readout chain is not found in calibration data!! Aborting')
		return 
	return	g_s , Tn_m_s ,Tn_p_s