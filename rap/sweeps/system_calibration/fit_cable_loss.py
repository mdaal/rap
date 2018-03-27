from ..visualization.plot_transmission import plot_transmission


import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

def fit_cable_loss(metadata, loop, key_name, freq_range = [400e6, 1e9], Verbose = True, Show_Plot = True):
    '''produces fit to cable loss in the functional form:
    term1 + term2 + term3 = a * sqrt(f) + b * f + c
    term1 is the sum of inner and outer coaxial cable conductor losses
    term2 is due to coaxial cable dielectric loss
    term3 is a constant fudge factor
    The loss evaluates to units of dB.

    stores the  fit as dictionary
    (a,b,c,run,range_start,range_stop)= self.metadata.Cable_Calibration['One_Way_40mk']

    Two used this function load transmission for complete cable loop only (not amps or attens).
    Then call this function on that transmission data. This funciton creats the tuple (a,b,c,run,range_start,range_stop) in
    metadata, where run is the name of the calibration run and range_start/stop is the frequency range over which the
    calibration is calculated.

    Create a function from a,b,c and it to the effect of attenuators on the input side of the cable loop.

    set freq_range = None to use full freq range
    '''

    f   = loop.freq
    s21 = loop.z

    if freq_range == None:
        condition = f == f
    else:
        condition = (f>freq_range[0]) & (f<freq_range[1])

    f = np.extract(condition,f)
    s21 = np.extract(condition,s21)



    def obj(x,s21,f):
        a,b,c = x
        return np.square(20*np.log10(np.abs(s21)) - a*np.sqrt(f) - b*f - c).sum() #attenuation in dB/length goes as -a*sqrt(f)-b*f-c, where c has no theoretical basis.

    p0 = np.array([-3.0e-4,-1.0e-9 ,0.5])

    res = minimize(obj, p0, args=(s21,f), method='Nelder-Mead', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=1e-15, callback=None, options={'disp':False, 'xtol' : 1e-6,'maxfev':1000})

    k = list(res.x/2.0) #devide by 2 to get one way loss
    k = k + [metadata.Run, f[0], f[-1]]

    if metadata.Cable_Calibration == None:
        cal = {}
        cal[key_name] = tuple(k)
        metadata.Cable_Calibration = cal
    else:
        metadata.Cable_Calibration[key_name] =tuple(k)

    if Verbose == True:
        print(res)

    if Show_Plot == True:
        (fig,ax,) = plot_transmission(metadata,loop,show = False)[:2]
        Cal  = lambda f: k[0]*np.sqrt(f)+k[1]*f+k[2]
        line = ax.plot(f, Cal(f)*2.0, 'r--', linewidth=3, label = 'fit - round trip')
        line = ax.plot(f, Cal(f), 'g-', linewidth=3, label = 'fit - one way')
        ax.set_xlim([freq_range[0]*0.75, freq_range[1]*1.25])
        leftvline = ax.axvline(x = freq_range[0],linewidth=2, color='k', linestyle = ':')
        rightvline = ax.axvline(x = freq_range[1],linewidth=2, color='k', linestyle = ':')
        ax.legend(loc = 'best', fontsize=10,scatterpoints =1, numpoints = 1, labelspacing = .1)
        plt.show()
