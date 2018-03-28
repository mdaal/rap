from .utils import _nonlinear_formulae
from ...sweep_array.pick_loop import pick_loop
from ...loop import loop
from ...visualization.utils import _save_fig_dec
from ..remove_cable_delay import remove_cable_delay
from ..decompress_gain  import decompress_gain
from ..circle_fit import circle_fit
from ...system_calibration.utils import _construct_readout_chain


import time
from scipy.stats import chisquare
import numpy as np
from scipy import constants
from scipy.optimize import minimize
import sys

def nonlinear_fit(metadata, loop, Sweep_Array, Fit_Method = 'Multiple', Verbose = True, Show_Plot = True, Save_Fig = False, Compute_Chi2 = False, Indexing = (None,None,None)):
    '''
    The indexing keyword allows for selection of the power sweep to be fit.
    If P is the list of powers then Indexing = (Start,Stop,Step) is using only, P[Start,Stop, Step]
    '''


    R = 50 #System Impedance
    k = constants.value('Boltzmann constant') #unit is [J/k]
    BW = metadata.IFBW #unit is [Hz]

    if isinstance(Fit_Method,str): #Allow for single string input for Fit_Method
       Fit_Method={Fit_Method}

    if loop.index == None:
        print('Loop index not chosen. Setting to 0.')
        index = 0
        pick_loop(loop,Sweep_Array,index)

    Sweep_Array_Record_Index = loop.index
    V = Sweep_Array['Heater_Voltage'][Sweep_Array_Record_Index]
    Fs = Sweep_Array['Fstart'][Sweep_Array_Record_Index]

    #### NOTE:  will need to fix for the case of sweeps with  duplicate V .... will involve using np.unique
    indices = np.where( (Sweep_Array['Heater_Voltage'] == V) & (Sweep_Array['Fstart']==Fs))[0]
    P_min_index = np.where( (Sweep_Array['Heater_Voltage'] == V) & ( Sweep_Array['Fstart']==Fs) & (Sweep_Array['Pinput_dB'] == Sweep_Array['Pinput_dB'].min()))[0][0]

    ##### Q, Qc, Qtl, fr  - used for initial guess in minimization
    ##### Zfl, Zres - used in minimization, Zfl converts power to voltage
    Q   = Sweep_Array['Q'][P_min_index]
    Qc  = Sweep_Array['Qc'][P_min_index]

    Qtl = np.power( (1./Q) - (1./Qc) , -1.)
    fr = Sweep_Array['Fr'][P_min_index]
    Zfl = metadata.Feedline_Impedance
    Zres = metadata.Resonator_Impedance


    power_sweep_list = []
    invalid_power_sweep_list = []
    start, stop, step = Indexing
    for index in indices[start:stop:step]: #
        loop_i = loop()

        # Pick new loop_i
        pick_loop(loop_i,Sweep_Array,index)


        # Remove Gain Compression
        decompress_gain(Sweep_Array, loop_i, metadata, Compression_Calibration_Index = -1, Show_Plot = True, Verbose = True)

        # Normalize Loop
        Outer_Radius = Sweep_Array['R'][index]
        if (Outer_Radius <= 0) or (Outer_Radius == None):
            print('Outer loop radius non valid. Using 1')
            Outer_Radius  = 1
        loop_i.z = loop_i.z/Outer_Radius

        # Remove Cable Delay
        remove_cable_delay(loop_i, metadata, Show_Plot = True, Verbose = True, center_freq = None, Force_Recalculate = False)

        # Fit loop to circle
        circle_fit(loop_i)

        Preadout = 0.001*np.power(10, Sweep_Array['Preadout_dB'][index]/10.0) #W, Readout power at device
        V1 = np.sqrt(Preadout*2*Zfl) #V, Readout amplitude at device
        mask = Sweep_Array['Mask'][index]
        f = ma.array(loop_i.freq,mask = mask)
        z = ma.array(loop_i.z,mask = mask)
        zc = np.complex(loop_i.a,loop_i.b)
        z = z*np.exp(np.complex(0,-np.angle(zc))) #rotate to real axis, but dont translate to origin

        f_c = f.compressed()
        z_c = z.compressed()

        P_NA_out_dB = Sweep_Array[index]['Pinput_dB'] #Power out of the network analyzer, change of reference point
        P_NA_out_V2 = .001 * np.power(10,P_NA_out_dB/10) * 2 * R  #Voltage squared out of network analyzer

        if Compute_Chi2 is True: # Calculate variances for Chi2

            #z_c = z_c*Outer_Radius
            P_NA_in_V2 = np.square(np.abs(z_c)) * P_NA_out_V2
            g_s , Tn_m_s ,Tn_p_s = _construct_readout_chain(metadata, f_c) # get the gain chain

             # g_i is the total gain between the device and readout digitizer (Network Analyzer) at the frequency f_i
            sigma_squared_m = np.zeros_like(f_c)
            sigma_squared_p = np.zeros_like(f_c)
            n = len(g_s)

            for i in range(n):
                s2_m = 4*k*Tn_m_s[i]*R*BW # This sigma for the particular stage of the readout chain
                s2_p = 4*k*Tn_p_s[i]*R*BW
                #we assume s2_p * 4 * P_NA_in_V2 = s2_m ,  s2_p measured in radian^2
                sigma_squared_m = sigma_squared_m + s2_m*np.prod(g_s[i:], axis = 0) #rememebr g is a list of np vectors
                sigma_squared_p = sigma_squared_p + s2_p*np.square(np.prod(g_s[i:], axis = 0))/(4*P_NA_in_V2) #rememeber P_NA_in_V2 is a function of S21, see above definition
        else:
            sigma_squared_m = np.ones_like(f_c)
            sigma_squared_p = np.ones_like(f_c)


        if Sweep_Array['Is_Valid'][index] == True:
            power_sweep_list.append((V1,z_c,f_c,sigma_squared_m,sigma_squared_p,P_NA_out_V2,Outer_Radius))
        else:
            invalid_power_sweep_list.append((V1,z_c,f_c,sigma_squared_m,sigma_squared_p,P_NA_out_V2,Outer_Radius ))



    def progress(x):
        ''' Add a dot to stdout at the end of each iteration without removing the dot from the previous iteration or
        adding a new line.
        '''
        sys.stdout.write('.')
        sys.stdout.flush()


    V30V30 = fr #minimization will not converge if V30V30 is too small
    phiV1 = 0.0
    def obj(p):
        ''' *** Objective function to be minimized for Chi2 and other fit ***
        '''
        parameter_dict = {'f_0':p[0], 'Qtl':p[1], 'Qc':p[2], 'phi31':p[3], 'eta':p[4], 'delta':p[5], 'Zfl':Zfl, 'Zres':Zres,  'phiV1':phiV1, 'V30V30':V30V30}
        fd = _nonlinear_formulae( parameter_dict, model = 2) # get the nonlinear formulae dict, fd
        a,b,phi,tau = p[6:] # geometrical transformation parameters and tau - cable delay

        sumsq = 0
        N = 0 # total number of points in fit
        for sweep in power_sweep_list:
            V1_readout, S21_data, f,sigma_squared_m,sigma_squared_p,P_NA_out_V2 ,Outer_Radius= sweep

            V3 = fd['V3'](S21_data,V1_readout)
            v1 = V3*V3.conjugate()


            #  Compute S21 and then Impose geometrical transformations to on it
            S21_fit = (fd['S21'](v1,f) -  np.complex(a,b))/np.exp(np.complex(0,phi)+ np.complex(0,2.0*np.pi*tau)*f)

            if Compute_Chi2 is True:

                # Phase Mag approach doe not converge
                diff = np.square(( np.abs(S21_data) -  np.abs(S21_fit) ) * Outer_Radius)*P_NA_out_V2/sigma_squared_m  + np.square(np.angle(S21_data/S21_fit))/sigma_squared_p  #(e^ia)/(e^ib) = e^i(a-b)

                # Real Imaginary approach does not converge
                #diff = np.square(S21_data.real -  S21_fit.real)*P_NA_out_V2/sigma_squared_m  + np.square(S21_data.imag -  S21_fit.imag)*P_NA_out_V2/sigma_squared_m

                # Real Imaginary approach does not *without# P_NA_out_V2 does converge!
                #diff = np.square(S21_data.real -  S21_fit.real)/sigma_squared_m  + np.square(S21_data.imag -  S21_fit.imag)/sigma_squared_m

                sumsq = diff.sum()  + sumsq
                N = N + f.shape[0]*1.0
            else:
                diff = S21_data - S21_fit
                sumsq = (diff*diff.conjugate()).real.sum()  + sumsq

        if Compute_Chi2 is True:
            return sumsq/(N-p.shape[0])
        else:
            return sumsq




    phi31_est = np.pi/2
    eta_est = 0.001
    delta_est = 0.001
    a_est = 0.
    b_est = 0.
    phi_est = 0.
    tau_est = 0.0
    p0 = np.array([fr,Qtl,Qc,phi31_est,eta_est,delta_est,a_est,b_est, phi_est,tau_est ])
    #Each fit method is saved as a lambda function in a dictionary called fit_func
    fit_func = {}
    fit_func['Powell'] = lambda : minimize(obj, p0, method='Powell', jac=None, hess=None, hessp=None, bounds=None, constraints=None, tol=1e-20, callback=progress, options={'disp':False, 'maxiter': 100, 'maxfev': 50000, 'ftol':1e-14,'xtol':1e-14}) #maxfev: 11137 defaults: xtol=1e-4, ftol=1e-4,
    #fit_func['Nelder-Mead']  = lambda : minimize(obj, p0, args=(z_theta_c,f_c), method='Nelder-Mead', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=1e-15, callback=None, options={'disp':False, 'xtol' : 1e-6,'maxfev':1000})
    #fit_func['Newton-CG'] = lambda : minimize(obj, p0, args=(z_theta_c,f_c), method='Newton-CG', jac=jac, hess=hess, hessp=None, bounds=None, constraints=(),tol=1e-15, callback=None, options={'maxiter' : 50,'xtol': 1e-4,'disp':False})


    fit = {}
    start = time.time()

    for method in fit_func.keys():
        sys.stdout.write('Iterating')
        sys.stdout.flush()
        fit[method] = fit_func[method]()

    finished = time.time()
    elapsed = (finished - start )/60.0 #minutes
    print('Minimization took {:.2f} minutes'.format(elapsed))


    if fit.keys() != []: #if there is a fit object in the fit dictionary
        bestfit = list(fit)[0]
        lowest = fit[bestfit].fun # .fun is function value
        for key in fit.keys():
            if fit[key].fun < lowest:
                lowest = fit[key].fun
                bestfit = key
    else:
        bestfit = None



    if Verbose == True:
        print(fit[bestfit])

    if Show_Plot == True:
        #Determine Sweep Direction
        direction = 'up'
        if direction == 'up':
            #min |--> up sweep (like at UCB)
            extreme = np.min
        else:
            # max |--> down sweep
            extreme = np.max

        ####### Set up plot objects
        fig = plt.figure( figsize=(5, 5), dpi=150)
        ax = {}
        gs = gridspec.GridSpec(2, 2)
        ax[1] = plt.subplot(gs[0, :])
        ax[2] = plt.subplot(gs[1, 0], aspect='equal' )
        ax[3] = plt.subplot(gs[1, 1])
        note = (r'Run {run}, Resonator width {width:.0f} $\mu m$'+'\n').format(run = metadata.Run,
            width = (metadata.Resonator_Width if metadata.Resonator_Width is not None else 0)/1e-6)

        if bestfit != None:
            p = fit[bestfit].x
            parameter_dict = {'f_0':p[0], 'Qtl':p[1], 'Qc':p[2], 'phi31':p[3], 'eta':p[4], 'delta':p[5], 'Zfl':Zfl, 'Zres':Zres,  'phiV1':phiV1, 'V30V30':V30V30}
            fd = _nonlinear_formulae( parameter_dict, model = 2) # get the nonlinear formulae dict, fd
            a,b,phi,tau = p[6:]
            vline = ax[1].axvline(x = (parameter_dict['f_0']-fr)/fr,linewidth=1, color='y', linestyle = ':')#,   label = r'$f_{r}$')
            note = note + (r'$f_0$ = {f_0:3.2e} Hz, $Q_{sub1}$ = {Qtl:3.2e}, $Q_c$ = {Qc:3.2e}' +
                '\n' + r'$\phi_{sub2}$ = {ang:3.2f}$^\circ$, ${l1}$ = {et:3.2e}, ${l2}$ = {de:3.2e}').format(
                nl = '\n', et = parameter_dict['eta']/parameter_dict['V30V30'],
                de = parameter_dict['delta']/parameter_dict['V30V30'],
                l1 = r'{\eta}/{V_{3,0}^2}',
                l2  = r'{\delta}/{V_{3,0}^2}',
                ang = parameter_dict['phi31']*180/np.pi,
                sub1 = '{i}', sub2 = '{31}',**parameter_dict)


        for sweep in power_sweep_list:
            V1exp, S21exp, f ,sigma_squared_m,sigma_squared_p,P_NA_out_V2,Outer_Radius= sweep
            Pexp = 10*np.log10(V1exp*V1exp/(2 *Zfl*0.001))
            dff = (f - fr)/fr
            curve = ax[1].plot(dff,20*np.log10(np.abs(S21exp)), label = '$P_{probe}$ =' + ' {:3.2f} dBm'.format(Pexp)) # Pexp is Preadout
            curve = ax[2].plot(S21exp.real,S21exp.imag)


            if bestfit != None:
                #####Compute the experimental values of V3
                V3_exp = fd['V3'](S21exp,V1exp)

                #####Initialize arrays
                Number_of_Roots = 3
                V3V3 = np.ma.empty((f.shape[0],Number_of_Roots), dtype = np.complex128)
                V3V3_cubic = np.empty(f.shape)
                V3_cubic = np.empty(f.shape)
                S21_fit = np.empty_like(f,dtype = np.complex128)
                V3_fit = np.empty_like(f,dtype = np.complex128)

                for n in range(f.shape[0]):
                    coefs = np.array([fd['z1z1'](f[n]), 2*fd['rez1z2c'](f[n]), fd['z2z2'](f[n]), -fd['z3z3'](V1exp)])
                    V3V3[n] =np.ma.array(np.roots(coefs),mask= np.iscomplex(np.roots(coefs)),fill_value = 1)
                    V3V3_cubic[n]    = extreme(np.extract(~V3V3[n].mask,V3V3[n])).real
                    V3_cubic[n]    = np.sqrt(V3V3_cubic[n])
                    # S21_fit is adjused to take into accout fit parameters a,b,phi,tau
                    S21_fit[n]  = (fd['S21'](V3V3_cubic[n],f[n]) - np.complex(a,b))*np.exp(np.complex(0,-phi)+ np.complex(0,-tau*2.0*np.pi)*f[n])

                    # Note that V3_fit has the effect of a,b,phi,tau incorporated,
                    # So it should no be expected to equal V3_cubic
                    V3_fit[n] = fd['V3'](S21_fit[n],V1exp)

                S21_cor = np.complex(a,b)+ np.exp(np.complex(0,phi)+ np.complex(0,2.0*np.pi*tau)*f)*S21exp
                V3_cor  = fd['V3'](S21_cor,V1exp)

                curve = ax[1].plot(dff,20*np.log10(np.abs(S21_fit)), linestyle = ':', color = 'c')
                curve = ax[2].plot(S21_fit.real,S21_fit.imag, linestyle = ':', color = 'c')

                # curve = ax[3].plot(dff.real,V3_cor.real)
                # curve = ax[3].plot(dff.real,V3_cubic.real, linestyle = ':', color = 'g')


                # curve = ax[3].plot(dff,V3_exp.real)
                # curve = ax[3].plot(dff.real,V3_fit.real, linestyle = ':', color = 'c')#~np.iscomplex(V3fit)

                curve = ax[3].plot(dff,np.abs(V3_exp))
                curve = ax[3].plot(dff.real,np.abs(V3_fit), linestyle = ':', color = 'c')

        ax[1].set_title('Mag Transmission')
        ax[1].set_xlabel(r'$\delta f_0 / f_0$', color='k')
        ax[1].set_ylabel(r'$20 \cdot \log_{10}|S_{21}|$ [dB]', color='k')
        ax[1].yaxis.labelpad = 0 #-6
        ax[1].xaxis.labelpad = 3
        ax[1].ticklabel_format(axis='x', style='sci',scilimits = (0,0), useOffset=True)
        ax[1].text(0.01, 0.01, note,
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax[1].transAxes,
            color='black', fontsize=4)
        ax[1].legend(loc = 'upper center', fontsize=5, bbox_to_anchor=(.5, -1.5),  ncol=4,scatterpoints =1, numpoints = 1, labelspacing = .02)
        #bbox_to_anchor=(1.25, -0.1),bbox_transform = ax[2].transAxes,



        ax[2].set_title('Resonance Loop')
        ax[2].set_xlabel(r'$\Re$[$S_{21}$]', color='k')
        ax[2].set_ylabel(r'$\Im$[$S_{21}$]', color='k')
        ax[2].yaxis.labelpad = -4
        ax[2].ticklabel_format(axis='x', style='sci',scilimits = (0,0),useOffset=False)

        ax[3].set_title('Resonator Amplitude')
        ax[3].set_xlabel(r'$\delta f_0 / f_0$', color='k')
        ax[3].ticklabel_format(axis='x', style='sci',scilimits = (0,0),useOffset=False)

        mpl.rcParams['axes.labelsize'] = 'small' # [size in points | 'xx-small' | 'x-small' | 'small' | 'medium....

        for k in ax.keys():
            ax[k].tick_params(axis='y', labelsize=5)
            ax[k].tick_params(axis='x', labelsize=5)

        plt.subplots_adjust(left=.1, bottom=.1, right=None ,wspace=.35, hspace=.3)

        if Save_Fig == True:
            name  = 'Nonlinear_Fit_'
            if Compute_Chi2 is True:
                name = name  + 'Chi2_'
            _save_fig_dec(metadata,fig, name + 'Start_Index_'+ str(Sweep_Array_Record_Index))
        plt.subplots_adjust(top =0.90)
        plt.suptitle('Fit to Nonlinear Resonator Data', fontweight='bold')
        plt.show()


        fit.update(phiV1= phiV1, V30V30= V30V30)
    return fit, fig, ax #need to figure out a way to return all the curves too
