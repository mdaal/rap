from .utils import _angle, _points_removed

from scipy.stats import chisquare
import numpy as np
import numpy.ma as ma
from scipy.optimize import minimize
import matplotlib.pyplot as plt

def phase_fit(loop, env_var, Fit_Method = 'Multiple', Verbose = True, Show_Plot = True):
	'''
	Note: its best to determine angles and angle differences by starting with complex numbers 
	(interpreted as vectors) and then finding their angles with, np.angle or _angle. It is
	not as accurate and prone to issues with domains (e.g. [-180,180]) to use arcsin or arccos.
	'''

	
	if isinstance(Fit_Method,str): #Allow for single string input for Fit_Method
	   Fit_Method={Fit_Method}
	   
	
	j = np.complex(0,1)
	try:
		zc = loop.a + j*loop.b
		r = loop.r
	except:
		print('Phase fit needs loop center and radius, which are not currently defined. Aborting phase fit.')
		return



	f = f0 = loop.freq
	z = z0 = loop.z
	

	# Remove duplicate frequency elements in z and f, e.g. places where f[n] = f[n+1]
	f_adjacent_distance =  np.hstack((np.abs(f[:-1]-f[1:]), [0.0]))
	z = z1 = ma.masked_where(f_adjacent_distance==0.0, z)
	f = f1 = ma.array(f,mask = z.mask) #Syncronize mask of f to match mask of z
	


	#Estimate Resonance frequency using minimum Dip or max adjacent distance
	Use_Dip = 1 
	if Use_Dip: #better for non-linear resonances with point near loop center
		zr_mag_est = np.abs(z).min()
		zr_est_index = np.where(np.abs(z)==zr_mag_est)[0][0]
	else:
		z_adjacent_distance = np.abs(z[:-1]-z[1:])
		zr_est_index = np.argmax(z_adjacent_distance) 
		zr_mag_est = np.abs(z[zr_est_index])


	#Transmission magnitude off resonance 
	Use_Fit = 1
	if Use_Fit:
		z_max_mag = np.abs(zc)+r
	else: #suspected to be better for non-linear resonances
		z_max_mag = np.abs(z).max()

	#Depth of resonance in dB
	depth_est = 20.0*np.log10(zr_mag_est/z_max_mag)

	#Magnitude of resonance dip at half max
	res_half_max_mag = (z_max_mag+zr_mag_est)/2

	#find the indices of the closest points to this magnitude along the loop, one below zr_mag_est and one above zr_mag_est
	a = np.square(np.abs(z[:zr_est_index+1]) - res_half_max_mag)
	lower_index = np.argmin(a)#np.where(a == a.min())[0][0]
	a = np.square(np.abs(z[zr_est_index:]) - res_half_max_mag)
	upper_index = np.argmin(a) + zr_est_index

	#estimate the FWHM bandwidth of the resonance
	f_upper_FWHM = f[upper_index]
	f_lower_FWHM = f[lower_index]
	FWHM_est = np.abs(f_upper_FWHM - f_lower_FWHM)
	fr_est = f[zr_est_index]

	
	#consider refitting the circle here, or doing ellipse fit.



	#translate circle to origin, and rotate so that z[zr_est_index] has angle 0 
	z = z2 = ma.array((z.data-zc)*np.exp(-j*(_angle(zc))), mask = z.mask)

	#Compute theta_est before radious cut to prevent radius cut from removing z[f==fr_est]
	theta_est = _angle(z[zr_est_index]) #_angle(z[zr_est_index])	

	#Radius Cut: remove points that occur within r_cutoff of the origin of the centered data. 
	#(For non-linear resonances that have spurious point close to loop center)	
	r_fraction_in = 0.75
	r_fraction_out = 1.75
	r_cutoff_in  = r_fraction_in*r
	r_cutoff_out = r_fraction_out*r		
	z = z3 = ma.masked_where((np.abs(z2)<r_cutoff_in) | (np.abs(z2)>r_cutoff_out),z2, copy = True)
	# for substantially deformed loops we make sure that no more than Max_Removed_Radius_Cut points are removed from inner radious cut
	Max_Removed_Radius_Cut = 25
	while _points_removed(z2, z3)[0] > Max_Removed_Radius_Cut:
		r_fraction_in = r_fraction_in - 0.02
		r_cutoff_in  = r_fraction_in*r
		z = z3 = ma.masked_where((np.abs(z2)<r_cutoff_in) | (np.abs(z2)>r_cutoff_out),z2, copy = True)
		print 'loosening inner radius cut: r_fraction_in = {}'.format(r_fraction_in)
		if r_fraction_in <= 0:
			break
	f = f3 = ma.array(f,mask = z.mask)


	#Bandwidth Cut: cut data that is more than N * FWHM_est away from zr_mag_est
	N = 8
	z = z4 = ma.masked_where((f > fr_est + N*FWHM_est) | (fr_est - N*FWHM_est > f),z,copy = True)
	f = f4 = ma.array(f,mask = z.mask)
	z_theta,z_theta_offset =_angle(z, return_offset = True) 


	#Angle jump cut : masks points where angle jumps to next branch of angle function, 
	mask = (f > fr_est + 0.5*FWHM_est) | (f < fr_est + -0.5*FWHM_est)
	f_in_FWHM = ma.masked_where(mask,f) # or alternatively: f_in_FWHM = f; f_in_FWHM[mask] = ma.masked 
	edge1,edge2 = ma.flatnotmasked_edges(f_in_FWHM)
	angle_slope = (z_theta[edge2]-z_theta[edge1])/(f[edge2]-f[edge1]) # angle is decreasing if negative slope
	upper_cond = ((f > fr_est +  0.5*FWHM_est) & ((z_theta[edge2]<z_theta) if (angle_slope<0) else (z_theta[edge2]>z_theta))) 
	lower_cond = ((f < fr_est + -0.5*FWHM_est) & ((z_theta[edge1]>z_theta) if (angle_slope<0) else (z_theta[edge1]<z_theta))) 
	z = z5 = ma.masked_where(lower_cond|upper_cond,z, copy = True)
	f = f5 = ma.array(f,mask = z.mask)
	z_theta = z_theta5 = ma.array(z_theta,mask = z.mask)
	


	#theta_est = np.extract(f==fr_est,z_theta)[0] # The old lication of theta_est computation 
	Q_est = fr_est/FWHM_est


	#consider reducing computation by extracting only the unmasked values of z,f, and z_theta of the minimization
	#These commands return a masked array where all the masked elements are removed.
	#z = z[~z.mask]
	#f = f[~f.mask]
	#z_theta = z_theta[~z_theta.mask]

	#These commands return np array
	z_c = ma.compressed(z)
	f_c = ma.compressed(f)
	z_theta_c  = ma.compressed(z_theta)
	

	if env_var.mysys.startswith('Windows'):
		dt = np.float64
	else:	
		dt = np.float128

	def hess(x, z_theta,f): #to avoid overflow try to re write hessian so that all numbers are of order 1
		theta,fr,Q = x	
		H = np.zeros((3,3), dtype = dt)
		ff = (1-(f/fr))
		denom = (1+4.0*np.square(ff*Q))
		numer = (theta+z_theta-2.0*np.arctan(2.0*ff*Q))
		H[0,0] = (2.0*np.ones_like(z_theta)).sum()
		H[0,1] = ((-8.0*f*Q)/(np.square(fr)*denom)).sum()
		H[0,2] = ((8.0*ff)/denom).sum()
		H[1,0] = H[0,1] #((8.0*f*Q)/(np.square(fr)*denom)).sum()
		H[1,1] = ((32.0*np.square(f*Q/(np.square(fr)*denom)))  +   (64.0*np.square(f/(np.square(fr)*denom))*ff*np.power(Q,3)*numer)   +  ((16.0*f*Q/np.power(fr,3))*(numer/denom))).sum()
		H[1,2] = (((32.0*f*Q*ff)/np.square(fr*denom))  +  ((64.0*f*np.square(ff*Q)*numer)/(np.square(fr*denom)))  - ((8.0*f*numer)/(np.square(fr)*denom))).sum()
		H[2,0] = H[0,2] #((8.0*ff)/denom).sum()
		H[2,1] = H[1,2] #(((32.0*f*ff*Q)/np.square(fr*denom))  +  ((64.0*f*np.square(ff*Q)*numer)/(np.square(fr*denom)))  -  ((8.0*f*numer)/(np.square(fr)*denom))).sum()
		H[2,2] = (((32.0*np.square(ff))/np.square(denom))  +  ((64.0*np.power(ff,3)*Q*numer)/np.square(denom))).sum()				
		return H

	def jac(x,z_theta,f):
		theta,fr,Q = x
		J = np.zeros((3,),dtype = dt)    #np.zeros_like(x)
		ff = (1-(f/fr))
		denom = (1+4.0*np.square(ff*Q))
		numer = (theta+z_theta-2.0*np.arctan(2.0*ff*Q))	
		J[0] = np.sum(2.0*numer)
		J[1] = np.sum(-8.0*f*Q*numer/(np.square(fr)*denom))
		J[2] = np.sum(-8.0*ff*numer/denom)
		return J


	def obj(x,z_theta,f):
		theta,fr,Q = x
		return np.square(z_theta + theta - 2.0*np.arctan(2.0*Q*(1-f/fr))).sum()	 #<--- Need hessian of this


	def obj_ls(x,z_theta,f):
		'''object fuctinon for least squares fit'''
		theta,fr,Q = x
		residual  = z_theta + theta - 2.0*np.arctan(2.0*Q*(1-f/fr))	
		return residual

	#p0 is the initial guess
	p0 = np.array([theta_est,fr_est ,Q_est])
	
	#Each fit method is saved as a lambda function in a dictionary called fit_func
	fit_func = {}
	fit_func['Powell'] = lambda : minimize(obj, p0, args=(z_theta_c,f_c), method='Powell', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=1e-20, callback=None, options={'disp':False, 'maxiter': 70, 'maxfev': 50000, 'ftol':1e-20,'xtol':1e-20})#options={'disp':False})
	fit_func['Nelder-Mead']  = lambda : minimize(obj, p0, args=(z_theta_c,f_c), method='Nelder-Mead', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=1e-18, callback=None, options={'disp':False, 'xtol' : 1e-6,'maxfev':1000})
	fit_func['Newton-CG'] = lambda : minimize(obj, p0, args=(z_theta_c,f_c), method='Newton-CG', jac=jac, hess=hess, hessp=None, bounds=None, constraints=(),tol=1e-18, callback=None, options={'maxiter' : 50,'xtol': 1e-4,'disp':False})

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
	               				

	bestfit = list(fit)[0]
	lowest = fit[bestfit].fun
	for key in fit.keys(): 
		if fit[key].fun < lowest:
			lowest = fit[key].fun
			bestfit = key
	

	theta0 = 2*np.pi - _angle(np.exp(np.complex(0,fit[bestfit].x[0] - z_theta_offset)).conj())
	zc_m = np.abs(zc)
	R = np.sqrt(zc_m*zc_m + r*r -2.0*zc_m*r*np.cos(theta0) ) # Used in Qc

	alpha = _angle(zc)#np.angle(zc)#
	z_pivot = zc + (np.complex(-r*np.cos(theta0), r*np.sin(theta0)))*np.complex(np.cos(alpha),np.sin(alpha))# vector for origin to pizot point
	theta = _angle(z_pivot)
	phi  = np.angle(-(zc-z_pivot)*np.exp(-j*(_angle(z_pivot)))) #not that domain is [-180, +180]

	loop.R = R
	loop.phase_fit_success = fit[bestfit].success
	loop.phase_fit_z = z5.data
	loop.phase_fit_mask = z5.mask
	loop.phase_fit_method = bestfit
	loop.Q = Q = fit[bestfit].x[2]
	loop.Qc = Qc = Q*R/(2*r)
	loop.Qi = Q*Qc/(Qc-Q)
	loop.fr = fr = fit[bestfit].x[1]
	loop.FWHM = fr/Q
	loop.phi = phi # radian
	loop.theta = theta # radian
	loop.chisquare, loop.pvalue = chisquare( z_theta_c,f_exp=fit[bestfit].x[0] + 2.0*np.arctan(2.0*Q*(1-f_c/fr)))
	loop.chisquare = loop.chisquare/ f_c.shape[0]
	#estimated quantities from MAG S21 
	loop.fr_est = fr_est
	loop.FWHM_est = FWHM_est
	loop.depth_est = depth_est
	loop.Q_est = Q_est
	

	if Verbose: 
		print('Duplicates cuts:\n\t{0} duplicate frequencies removed from loop data, {1} remaining data points'.format(*_points_removed(z0,z1)))
		print('Radius cut:\n\t{2} points < r_loop*{0} or > r_loop*{1} found and removed, {3} remaining data points'.format(r_fraction_in, r_fraction_out,*_points_removed(z2,z3)))
		print('Bandwidth cut:\n\t{1} points outside of fr_est +/- {0}*FWHM_est removed, {2} remaining data points'.format(N, *_points_removed(z3,z4)))
		print('Angle jump cut:\n\t{0} points with discontinuous jumps in loop angle removed, {1} remaining data points'.format(*_points_removed(z4,z5)))
		print('Initial Guess:\n\tLoop rotation {0} deg, fr {1}, Q {2}'.format(p0[0]*180/np.pi,p0[1],p0[2] ))

		for method in fit.keys():
			print('\n{0} Minimzation Result:\n{1}\n'.format(method,fit[method]))



	if Show_Plot:
		total_removed, total_used_in_fit = _points_removed(z0,z5)
		fig1 = plt.figure( facecolor = 'w',figsize = (10,10))
		ax = fig1.add_subplot(6,1,1)
		ax.set_title('Number of points used in fit = '+str(total_used_in_fit)+', Number of points removed = ' + str(total_removed) )
		#line = ax.plot(f1[~f5.mask], np.abs(z1[~z5.mask]),'g-', label = 'Used for Fit') #fails when no points are masked
		


		if f5.mask.size <= 1:#this is the case that there are no masked points, e.g. no mask. there will allways be 1 point in the mask due to adjacent distance
			line = ax.plot(ma.compressed(f1), np.abs(ma.compressed(z1)),'g-', label = 'Used for Fit')
		else:
			line = ax.plot(f1[~f5.mask], np.abs(z1[~z5.mask]),'g-', label = 'Used for Fit')
			line = ax.plot(f1[f5.mask], np.abs(z1[z5.mask]),'r.',markersize = 2,  alpha = 0.2, label = 'Excluded Data')
		line = ax.plot([f1[zr_est_index],f1[zr_est_index]] , [np.abs(z1[zr_est_index]),np.abs(zc)+r] ,'k.', label = 'Magitude Min and Max')
		line = ax.plot([f1[lower_index], f1[upper_index], f1[upper_index]], np.abs([z1[lower_index],z1[lower_index],z1[upper_index]]),'yo-', label = 'FWHM Estimate')
		ax.set_ylabel('Magnitude')
		## Find index of closet freq point to Fr
		a = np.square(np.abs(f1 - fr))
		fr_index = np.argmin(a)
		line = ax.plot(f1[fr_index], np.abs(z1[fr_index]),'gx', markersize = 7, markeredgewidth = 4, label = 'Fr (closest)')# this is the closest point in  the cut z1 to the true fr 
		ax.legend(loc = 'best', fontsize=10,scatterpoints =1, numpoints = 1, labelspacing = .1)
		
		
		ax = fig1.add_subplot(6,1,(2,4), aspect='equal')
		t = np.linspace(0, 2.0*np.pi, num=50, endpoint=True)
		line = ax.plot([0,zc.real],[0, zc.imag],'y*-', label = 'Center Vector')	
		line = ax.plot(zc.real + r*np.cos(t),zc.imag + r*np.sin(t),'y-', label = 'Circle Fit')		
		line = ax.plot(z1.real, z1.imag,'r:', label = 'Initial Location')
		line = ax.plot(z3.real, z3.imag,'r-', label = 'Aligned w/ Origin')
		lint = ax.plot([0,z_pivot.real],[0,z_pivot.imag],'yo-', label = 'Pivot point')
		lint = ax.plot([zc.real,z_pivot.real],[zc.imag,z_pivot.imag],'yo-', label = '_zc_to_zp')#zp is zpivot
		## Find index of closet freq point to Fr
		a = np.square(np.abs(f_c - fr))
		fr_index = np.argmin(a)

		line = ax.plot(z_c[fr_index].real, z_c[fr_index].imag,'gx', markersize = 7, markeredgewidth = 4, label = 'Fr (closest)')
		line = ax.plot([0,r*np.cos(theta0)],[0,-r*np.sin(theta0)], 'b',  label = 'Fr (True)') #vector to fr
		
		line = ax.plot(z4.real, z4.imag,'g:', linewidth = 3,label = 'Bandwidth Cut')
		##pt = ax.plot([z1[0].real,z[~z.mask][0].real], [z1[0].imag,z[~z.mask][0].imag],'ko', label = 'First Point') fails when no points are masked
		pt = ax.plot([z1[0].real,ma.compressed(z5)[0].real], [z1[0].imag,ma.compressed(z5)[0].imag],'ko', label = 'First Point') #--
		pt = ax.plot(z2[zr_est_index].real, z2[zr_est_index].imag,'k*', label = 'Magnitude Min')

		#line = ax.plot(z4[z4.mask].data.real, z4[z4.mask].data.imag,'r.', alpha = 0.2, label = 'Excluded Data')
		line = ax.plot(z5[ma.getmaskarray(z5)].data.real, z5[ma.getmaskarray(z5)].data.imag,'r.', alpha = 0.2,label = 'Excluded Data')
		ax.legend(loc = 'center left', bbox_to_anchor=(1.01, 0.5), fontsize=10, scatterpoints =1, numpoints = 1, labelspacing = .1)#,numpoints)
		
		text = ('$*Resonator Properties*$\n' + '$Q =$ ' + '{0:.2f}'.format(loop.Q) +'\nf$_0$ = ' + '{0:.6f}'.format(loop.fr/1e6) 
			+  ' MHz\n$Q_c$ = ' + '{0:.2f}'.format(loop.Qc) + '\n$Q_i$ = ' + '{0:.2f}'.format(loop.Qi) + '\n|S$_{21}$|$_{min}$ = ' 
			+ '{0:.3f}'.format(loop.depth_est) + ' dB' + '\nBW$_{FWHM}$ = ' + '{0:.3f}'.format(loop.FWHM/1e3) +  ' kHz' 
			+ '\n$\chi^{2}$ = ' + '{0:.4f}'.format(loop.chisquare) + '\n$\phi$ = ' + '{0:.3f}'.format(loop.phi*180/np.pi) +' deg' + '\n' + r'$\theta$ = ' 
			+ '{0:.3f}'.format(loop.theta*180/np.pi) +' deg' +'\n$- $'+loop.phase_fit_method 
			+ ' fit $-$') 
		bbox_args = dict(boxstyle="round", fc="0.8")        
		fig1.text(0.10,0.7,text,
				ha="center", va="top", visible = True,
				bbox=bbox_args, backgroundcolor = 'w')
		# ax.text(0.01, 0.01, text,
		# 	verticalalignment='bottom', horizontalalignment='left',
		# 	transform=ax.transAxes,
		# 	color='black', fontsize=4)


		ax = fig1.add_subplot(6,1,5)
		hline = ax.axhline(y = fit[bestfit].x[0],linewidth=2, color='y', linestyle = '-.',   label = r'$\theta_{r}$')
		vline = ax.axvline(x = fit[bestfit].x[1],linewidth=2, color='y', linestyle = ':',   label = r'$f_{r}$')
		line = ax.plot(f,z_theta,'g-',linewidth = 3,label = 'Data')
		line = ax.plot(f,(-fit[bestfit].x[0] + 2.0*np.arctan(2.0*fit[bestfit].x[2]*(1-f/fit[bestfit].x[1]))),'g:', linewidth = 1, label = 'Fit ')
		#line = ax.plot(f5[~f5.mask][0],z_theta5[~z_theta5.mask][0],'ko',linewidth = 3,label = 'First Point') #Failes when  no points are masked
		line = ax.plot(ma.compressed(f5)[0],ma.compressed(z_theta5)[0],'ko',linewidth = 3,label = 'First Point')

		ax.set_ylabel('Angle [rad]')
		ax.legend(loc = 'right', fontsize=10,scatterpoints =1, numpoints = 1, labelspacing = .1)
		
		ax = fig1.add_subplot(6,1,6)
		vline = ax.axvline(x = fit[bestfit].x[1],linewidth=2, color='y', linestyle = ':',   label = r'$f_{r}$')
		style  = ['-','--',':','-.','+','x']; s = 0 #Cyclic iterable?
		for key in fit.keys():
			line = ax.plot(f,(z_theta - fit[key].x[0] - 2.0*np.arctan(2.0*fit[key].x[2]*(1-f/fit[key].x[1]))),'b'+style[s], linewidth = 3, label = 'Data - Fit ' + key)
			s += 1
		ax.set_ylabel('Angle [rad]\nresiduals')
		ax.set_xlabel('Freq [Hz]')
		ax.legend(loc = 'right', fontsize=10,scatterpoints =1, numpoints = 1, labelspacing = .1)
		plt.show()

		# fig = plt.figure( figsize=(5, 5), dpi=150)
		# ax = {}
		# ax[1] = fig.add_subplot(1,1,1)
		# #dff = (f5 - fr)/fr
		# dff = f5 
		# curve = ax[1].plot(dff,np.abs(z5))
		# ax[1].ticklabel_format(axis='x', style='sci',scilimits = (0,0), useOffset=True)	

		# for k in ax.keys():
		# 	ax[k].tick_params(axis='y', labelsize=9)
		# 	ax[k].tick_params(axis='x', labelsize=5)
		# plt.show()