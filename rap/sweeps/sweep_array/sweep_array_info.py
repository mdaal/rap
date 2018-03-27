import numpy as np

def sweep_array_info(Sweep_Array):
    ''' prints information about the Sweep_Array currently loaded'''
    Input_Powers = np.unique(Sweep_Array['Pinput_dB'])
    Heater_Voltages = np.unique(Sweep_Array['Heater_Voltage'])
    Temperature_Points = np.shape(Sweep_Array['Temperature_Readings'])[1]
    Number_of_Freq_Ranges = max(np.unique(Sweep_Array['Fstart']),np.unique(Sweep_Array['Fstop']))
    print('{0:03.0f} - Total number of sweeps.\n{1:03.0f} - Number of readout powers.\n{2:03.0f} - Number of readout temperatures.\n{3:03.0f} - Number of temperatures readings.\n{4:03.0f} - Number of frequency bands.'.format(
        Sweep_Array.shape[0],
        Input_Powers.shape[0],
        Heater_Voltages.shape[0],
        Temperature_Points,
        Number_of_Freq_Ranges.shape[0]))
