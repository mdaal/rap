import numpy as np


def fit_system_calibration(metadata):
    '''compute chebyshev polynomial fits for  gain and noise values.
    save resulting polynomial coefficients list as:

    self.metadata.System_Calibration['device'][x + '_fit']

    where x is [gain, Tn_m ,Tn_p]...

    use numpy.polynomial.chebyshev.chebval to evaluate polynomial


    '''
    max_degree = 9

    SC = metadata.System_Calibration
    #already_fit = [k + '_fit' for k in SC[key].keys()]

    # Dont fit 'freq' and 'P1dB' to 'freq'
    # dont fit specs which *are* fits already
    dont_fit  = set(['freq','P1dB'])
    for key in SC.keys():
        for spec in SC[key].keys():
            if spec.find('_fit') > -1:
                dont_fit.add(spec)

    for key in SC.keys():
        for spec in set(SC[key].keys()).difference(dont_fit): # everything in SC[key] except for dont_fit
            deg =  min(len(SC[key]['freq']) - 2,max_degree) if len(SC[key]['freq']) >2 else len(SC[key]['freq']) - 1
            coefs = np.polynomial.chebyshev.chebfit(SC[key]['freq'],  SC[key][spec], deg)
            #coefs = numpy.polynomial.polynomial.polyfit(SC[key]['freq'],  SC[key]['g'], deg)


            SC[key].update({spec + '_fit':list(coefs)})
