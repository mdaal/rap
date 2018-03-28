from ..data_management.utils import  _define_sweep_data_columns, _define_sweep_array
from ..system_calibration.utils import _construct_readout_chain
from ..fitting.nonlinear.utils import _nonlinear_formulae
from ..visualization.utils import _save_fig_dec

import numpy as np
import matplotlib.pyplot as plt
import sys


def generate_nonlinear_data(metadata, Show_Plot = True, Phase_Noise_Variance = None, Amplitude_Noise_Variance = None, Like = None, Save_Fig = False,
    curve_parameter_dict = {'f_0':700e6, 'Qtl':300e3, 'Qc':80e3, 'eta':1e-1, 'delta':1e-6, 'Zfl':30, 'Zres':50, 'phi31': np.pi/2.03, 'phiV1':np.pi/10, 'V30V30':0.01},
    sweep_parameter_dict = {'Run': 'F1', 'Pprobe_dBm_Start' :-65.0,'Pprobe_dBm_Stop': -25.0, 'Pprobe_Num_Points':10, 'numBW':40,'num': 2000, 'Up_or_Down': 'Up', 'Freq_Spacing':'Linear'}):
    '''Creates and Loads Nonlinear Data
    eta -- Q nonlinearity
    delta --  freq nonlinearity
    V30V30 -- V^2 normalization for nonlinearity

    If another KAM.sweep object is supplied in "Like" keyword, then its metadata will copied
    '''
    cd = curve_parameter_dict
    sd = sweep_parameter_dict


    # system_attenuation_before_device = -50 # dB,  Difference between Preadout and Pinput
    metadata.Electrical_Delay = 0.0
    metadata.Feedline_Impedance = cd['Zfl']
    metadata.Resonator_Impedance = cd['Zres']
    metadata.LNA = {}
    metadata.LNA['LNA'] =  'SiGe #1'
    metadata.RTAmp_In_Use = True
    metadata.Atten_At_4K = 40.
    metadata.Atten_NA_Output = 0.
    metadata.Atten_NA_Input = 0.
    Cable_Calibration_Key = 'One_Way_40mK'
    metadata.Cable_Calibration = {}
    metadata.Cable_Calibration[Cable_Calibration_Key] = (0,0,0, 'False Cable', 0, 100e9)
    metadata.IFBW = 1.0
    metadata.Num_Temperatures = 0.0
    metadata.Num_Ranges = 1.0
    metadata.Num_Powers = sd['Pprobe_Num_Points']
    metadata.Num_Heater_Voltages = 1.0
    metadata.Loop_Data_Column_Names  = ("Frequencies","S21")

    if Like is not None: #would be better to confrim that Like is an instance of KAM.sweep
            metadata.__dict__.update(Like.metadata.__dict__)

    metadata.Electrical_Delay = 0
    metadata.Time_Created =   '05/01/2015 12:00:00' # or the current datetime datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S')
    metadata.Run = sd['Run']



    Q = 1.0/ ((1.0/cd['Qtl']) + (1.0/cd['Qc']))


    ############################# Cable Calbration
    k = metadata.Cable_Calibration[Cable_Calibration_Key]
    Preadout = lambda f: k[0]*np.sqrt(f)+k[1]*f+k[2] - metadata.Atten_NA_Output - metadata.Atten_At_4K

    ############################## Make Pprobe Array - This is power at the device.
    Pprobe_dBm = np.linspace(sd['Pprobe_dBm_Start'],sd['Pprobe_dBm_Stop'], sd['Pprobe_Num_Points'])
    Pprobe = 0.001* np.power(10.0,Pprobe_dBm/10.0)
    V1V1 = Pprobe *2*cd['Zfl']
    V1 = np.sqrt(V1V1) * np.exp(np.complex(0,1)*cd['phiV1'])
    # NOTE : V1 has a phase. Its a complex number

    ################################# Create f array making sure it contains f_0
    BW = sd['numBW']*cd['f_0']/Q

    if sd['Freq_Spacing'].lower() == 'triangular': #Triangular numbers - Denser around f_0
        T = np.linspace(1, sd['num'],  num=sd['num'], endpoint=True, retstep=False, dtype=None)
        T = T*(T+1.0)/2.0
        f_plus = (T*(BW/2)/T[-1]) + cd['f_0']
        f_minus = (-T[::-1]/T[-1])*(BW/2) + cd['f_0']
        f = np.hstack((f_minus,cd['f_0'],f_plus))

    if sd['Freq_Spacing'].lower() == 'linear': #linear
        f_plus = np.linspace(cd['f_0'], cd['f_0'] + BW/2,  num=sd['num'], endpoint=True, retstep=False, dtype=None)
        f_minus = np.linspace(cd['f_0']- BW/2,cd['f_0'],   num=sd['num']-1, endpoint=False, retstep=False, dtype=None)
        f = np.hstack((f_minus,f_plus))


    if sd['Freq_Spacing'].lower() == 'log': #logerithmic - Denser around f_0, for wide band sweeps
        f_plus = np.logspace(np.log10(cd['f_0']), np.log10(cd['f_0'] + BW/2),  num=sd['num'], endpoint=True, dtype=None)
        f_minus = -f_plus[:0:-1] + 2*cd['f_0']
        f = np.hstack((f_minus,f_plus))

    #################### Initialize Arrays
    Number_of_Roots = 3
    V3V3 = np.ma.empty((f.shape[0],Number_of_Roots), dtype = np.complex128)

    V3V3_direction = np.empty(f.shape)
    S21_direction = np.empty_like(f,dtype = np.complex128)

    #################### Construct gain array
    if Like is not None: #would be better to confrim that Like is an instance of KAM.sweep
        g_s , Tn_m_s ,Tn_p_s = _construct_readout_chain(metadata, f, Include_NA = True, Include_4K_to_40mK = False)# get the gain chain
        g = np.prod(g_s, axis = 0) # results in a numpy array  that is the same length as f...  a again for each frequency
    else:
        g = np.ones_like(f)
    # g_i is the total gain between the device and readout digitizer (Network Analyzer) at the frequency f_i


    #################### Find index of f_0
    try:
        f_0_index = np.where(f == curve_parameter_dict['f_0'])[0][0]
    except:
        d2 = np.square(f - curve_parameter_dict['f_0'])
        f_0_index = np.argmin(d2)



    #################### Initialize and Configure Sweep_Array
    tpoints = 0
    sweep_data_columns_list, sweep_data_columns = _define_sweep_data_columns(metadata,f.size,tpoints)
    Sweep_Array = np.zeros(Pprobe_dBm.size, dtype = sweep_data_columns) #Sweep_Array holdes all sweep data. Its length is the number of sweeps


    fig = plt.figure( figsize=(5, 5), dpi=150)
    ax = {}
    ax[1] = fig.add_subplot(1,1,1)
    dff = (f - cd['f_0'])/cd['f_0']



    #Determine Sweep Direction
    if sd['Up_or_Down'].lower() == 'up':
        #min |--> up sweep (like at UCB)
        extreme = np.min
    else:
        # max |--> down sweep
        extreme = np.max


    print('Generating False Data...')
    for index in range(Pprobe_dBm.size):
        sys.stdout.write('\r {0} of {1} '.format(index+1, Pprobe_dBm.size))
        sys.stdout.flush()
        Phase_Noise = np.zeros_like(f) if Phase_Noise_Variance is None else np.random.normal(scale = np.sqrt(Phase_Noise_Variance), size=f.shape)
        Amplitude_Noise = np.zeros_like(f) if Amplitude_Noise_Variance is None else np.random.normal(scale = np.sqrt(Amplitude_Noise_Variance), size=f.shape)

        for n in range(f.shape[0]):
            #################### Solve for Resonator amplitude using formulae from
            fd = _nonlinear_formulae(cd, model = 2) #get the nonlinear formulae dict, fd
            coefs = np.array([fd['z1z1'](f[n]), 2*fd['rez1z2c'](f[n]), fd['z2z2'](f[n]), -fd['z3z3'](V1[index])])


            V3V3[n] =np.ma.array(np.roots(coefs),mask= np.iscomplex(np.roots(coefs)),fill_value = 1)
            V3V3_direction[n]    = extreme(np.extract(~V3V3[n].mask,V3V3[n])).real
            S21_direction[n]  = fd['S21'](V3V3_direction[n],f[n])

        S21 = S21_direction + Amplitude_Noise + np.complex(0,1)*Phase_Noise

        Pin_dB = Pprobe_dBm[index] - Preadout(cd['f_0'])

        ####################  Fill Sweep_Array
        _define_sweep_array(Sweep_Array, index,
                                    Fstart = f[0], #Hz
                                    Fstop = f[-1], #Hz
                                    S21 = S21*np.sqrt(g), # g is power gain. so sqrt(g) is voltage gain #should be np.sqrt(g)*Rprobe_V/Pin_V  <-- _V meand voltage
                                    Frequencies = f, #Hz
                                    Preadout_dB = Pprobe_dBm[index],
                                    Pinput_dB = Pin_dB,
                                    Is_Valid = True,
                                    Mask = np.zeros(f.shape, dtype=np.bool),
                                    Chi_Squared = 0,
                                    Fr = cd['f_0'], #Note! we are using the resonance freq of the lowest power S21 for all
                                    Q = Q,
                                    Qc = cd['Qc'],
                                    Heater_Voltage = 0.0,
                                    R = np.sqrt(g[f_0_index]) # remember,  V1 is the readout probe amplitude
                                    )

        curve = ax[1].plot(dff,20*np.log10(np.abs(S21)), linestyle = '-', label = '$P_{probe}$ = ' + '{0:.2f} dBm'.format(Pprobe_dBm[index]))


    ################ Configure Plot
    ax[1].set_xlabel(r'$\delta f_0 / f_0$', color='k')
    ax[1].set_ylabel(r'Mag[$S_{21}$]', color='k')
    ax[1].yaxis.labelpad = -4
    ax[1].ticklabel_format(axis='x', style='sci',scilimits = (0,0), useOffset=True)
    ax[1].legend(loc = 'right', fontsize=4,scatterpoints =1, numpoints = 1, labelspacing = .1)
    for k in ax.keys():
        ax[k].tick_params(axis='y', labelsize=9)
        ax[k].tick_params(axis='x', labelsize=5)


    if Save_Fig:
        if Like is not  None:
            like = '_Like_' + Like.metadata.Run
        else:
            like = ''
        _save_fig_dec(metadata,fig,'Mock_Data' + like)
    #plt.subplots_adjust(left=.1, bottom=.1, right=None, top=.95 ,wspace=.4, hspace=.4)
    ax[1].set_title('Mag Transmission')
    plt.suptitle('Nonlinear Resonator Plots')
    plt.show()

    return fig, ax, sweep_data_columns_list, sweep_data_columns, Sweep_Array
