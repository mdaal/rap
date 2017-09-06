import numpy as np

def _nonlinear_formulae(parameter_dict, model = 2):
	''' model 2 is paramterization based on input resonator amplitude V_3^-, e.g.: 
	parameter_dict = {'f_0':700e6, 'Qtl':300e3, 'Qc':80e3, 'eta':1e-1, 'delta':1e-6, 'Zfl':30, 'Zres':50, 'phi31': np.pi/2.03, 'phiV1':np.pi/10, 'V30V30':}
	'''
	d = parameter_dict
	k = {	'z1'     :  lambda f      : d['eta']/(d['Qtl']*d['V30V30']) + np.complex(0,1.0)*(2*d['delta']*f)/(d['V30V30']*d['f_0']),
			'z2'     :  lambda f      : (1.0/d['Qc']) + (1.0/d['Qtl']) + np.complex(0,2.0) *(f-d['f_0'])/d['f_0'],
			'z3'     :  lambda V1     : np.sqrt(np.complex(1,0)*d['Zres']/(np.pi * d['Qc'] *d['Zfl'])) * np.exp(np.complex(0,d['phi31'])) * V1 *  np.exp(np.complex(0,d['phiV1'])),
			'z1z1'   :  lambda f      : (k['z1'](f) * k['z1'](f).conjugate()).real,
			'z2z2'   :  lambda f      : (k['z2'](f) * k['z2'](f).conjugate()).real,
			'z3z3'   :  lambda V1     : (k['z3'](V1) * k['z3'](V1).conjugate()).real,
			'rez1z2c':  lambda f      : (k['z1'](f) * k['z2'](f).conjugate()).real,
			'imz1z2c':  lambda f      : (k['z1'](f) * k['z2'](f).conjugate()).imag,
			#'V3'     :  lambda S21,V1 : (S21 + (np.exp(np.complex(0,2.0*d['phi31'])) - 1.0)/2.0 )*V1*np.exp(np.complex(0,-1.0*d['phi31']))*np.sqrt(d['Zres']*d['Qc']/(d['Zfl']*np.pi)), # may have less rounding error 
			'V3'     :  lambda S21,V1 : (S21 + (np.exp(np.complex(0,2.0*d['phi31'])) - 1.0)/2.0 )*k['z3'](V1)*d['Qc']*np.exp(np.complex(0,-2.0*d['phi31'])),
			'S21'    :  lambda V3V3,f : ((1-np.exp(np.complex(0,2.0)*d['phi31']))/2 +( (1.0/d['Qc']) / (k['z2'](f) + k['z1'](f)*V3V3))*np.exp(np.complex(0,2.0)*d['phi31']))
		}
			#						   V3  = (S21 + (np.exp(np.complex(0,2.0*phi31)) - 1.0)/2.0 )*V1*np.exp(np.complex(0,-1.0*phi31))*np.sqrt(Z3*Qc/(Z1*np.pi))
			# #Now we use |V3V3|^2 = v1 to calculate the other two roots of the cubic, v2 and v3
			# v1 = V3*V3.conjugate()
			# term1 = -(z1z2c.real/z1z1) - v1/2.0
			# term2 = np.complex(0,1)*np.sqrt(4*z1z2c.imag*z1z2c.imag + 3*v1*v1*z1z1*z1z1 + 4*z1z1*z1z2c.real*v1)/(2*z1z1)
			# v2  = term1 + term2
			# v3  = term1 - term2

			# V3p = np.sqrt(v2)
			# V3m = np.sqrt(v3)

			# Note: Ztl can be removed from the calculation. In which case we use Pfl, (i.e. Vfl = sqrt(Pfl*Zfl*2)) 

	return k
