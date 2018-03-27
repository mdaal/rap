from .pick_loop import pick_loop
from ..loop import loop

# from ..fitting.normalize_loop import normalize_loop
from ..fitting.remove_cable_delay import remove_cable_delay
from ..fitting.decompress_gain  import decompress_gain
from ..fitting.circle_fit import circle_fit
from ..fitting.phase_fit import phase_fit
from ..fitting.complete_fit import complete_fit
from ..fitting.trim_loop import trim_loop
from ..data_management.utils import _define_sweep_array

import numpy as np
import sys


def fill_sweep_array(metadata, Sweep_Array, Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False, Complete_Fit = True , Remove_Gain_Compression = True, Verbose = True):


    if Compute_Preadout == True:
        needed = ('Atten_NA_Output', 'Atten_At_4K','Cable_Calibration')

        for quantities  in needed:
            if  metadata.__dict__[quantities] == None:
                if Verbose == True:
                    print('{0} metadate missing. Unable to compute Preadout. Setting to 0.'.format(quantities))
                Compute_Preadout = False

            Atten_NA_Output = metadata.Atten_NA_Output
            Atten_At_4K = metadata.Atten_At_4K
            Cable_Calibration_Key = 'One_Way_40mK'
            k = metadata.Cable_Calibration[Cable_Calibration_Key]

        if Fit_Resonances == False:
            if Verbose == True:
                print('Resonance fit not selected. Computation of Preadout_dB requires knowledge of resonance frequency and may not work.')


        if Compute_Preadout == True:
            Preadout = lambda f: k[0]*np.sqrt(f)+k[1]*f+k[2] - Atten_NA_Output - Atten_At_4K

    if Add_Temperatures == True:
        if metadata.Num_Temperatures < 1:
            Temperature_Calibration = metadata.Temperature_Calibration
            if (metadata.Fridge_Base_Temp != None) & (max(Sweep_Array['Heater_Voltage']) == min(Sweep_Array['Heater_Voltage'])): #& (self.Sweep_Array.size == 1):
                #This is usually the case of a survey or power sweep: done at base temp with no Heater power
                Sweep_Array['Temperature'][:] = metadata.Fridge_Base_Temp
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


    loop_i = loop()
    num_records = Sweep_Array.size
    for index in xrange(num_records):
        if Verbose == True:
            sys.stdout.write('\r {0} of {1} '.format(index+1, num_records))
            sys.stdout.flush()

        #set current loop
        pick_loop(loop_i,Sweep_Array,index)


        if Fit_Resonances == True:
            if Remove_Gain_Compression:
                # Remove Gain Compression
                decompress_gain(Sweep_Array, loop_i, metadata, Compression_Calibration_Index = -1, Show_Plot = False, Verbose = False)

            if loop_i.z.size > 5000:
                trim_loop(loop_i, N= 20, Verbose = False)

            # Normalize Loop
            #normalize_loop()

            # Remove Cable Delay
            remove_cable_delay(loop_i, metadata, Show_Plot = False, Verbose = False) # should do nothing if a delay is defined in metadata

            # Fit loop to circle
            circle_fit(loop_i) # no 'Show_Plot = False' when using ..fitting.circle_fit
            if loop_i.circle_fit_exit_code != 0:
                _define_sweep_array(Sweep_Array,index, Is_Valid = False)


            # Fit resonance parameters
            phase_fit(loop_i, env_var, Fit_Method = 'Multiple', Verbose = False, Show_Plot = False)

            _define_sweep_array(Sweep_Array, index,
                                            Q = loop_i.Q,
                                            Qc = loop_i.Qc,
                                            Fr = loop_i.fr,
                                            Mask = loop_i.phase_fit_mask,
                                            Chi_Squared = loop_i.chisquare,
                                            R = loop_i.R,
                                            r = loop_i.r,
                                            a = loop_i.a,
                                            b = loop_i.b,
                                            #Normalization  = self.loop_i.normalization,
                                            Theta = loop_i.theta,
                                            Phi = loop_i.phi)

            if Complete_Fit:
                complete_fit(Sweep_Array, metadata, loop_i, Use_Mask = True, Verbose = False , Show_Plot = False, Save_Fig = False, Sample_Size = 100, Use_Loop_Data = True)
                _define_sweep_array(Sweep_Array, index,
                                                cQ = loop_i.cQ,
                                                cQc = loop_i.cQc,
                                                cFr = loop_i.cfr,
                                                cPhi = loop_i.cphi,
                                                cTheta = loop_i.ctheta,
                                                cR = loop_i.cR,
                                                cChi_Squared = loop_i.cchisquare,
                                                cIs_Valid = loop_i.cphase_fit_success if Sweep_Array['Is_Valid'][index] else Sweep_Array['Is_Valid'][index],

                                                sQ = loop_i.sQ,
                                                sQc = loop_i.sQc,
                                                sFr = loop_i.sfr,
                                                sPhi = loop_i.sphi,
                                                sTheta = loop_i.stheta,
                                                sR = loop_i.sR,
                                                sChi_Squared = loop_i.schisquare,
                                                sIs_Valid = loop_i.sphase_fit_success if Sweep_Array['Is_Valid'][index] else Sweep_Array['Is_Valid'][index]
                                                )




            # Only execute if phase_fit_success is False to avoid setting Is_Valid true when it was previously set fulse for a different reason, e.g bad Temp data
            if loop_i.phase_fit_success == False:
                _define_sweep_array(Sweep_Array, index, Is_Valid = False)

        if Compute_Preadout == True:
            if loop_i.fr != None:
                _define_sweep_array(Sweep_Array, index, Preadout_dB = Sweep_Array['Pinput_dB'][index] + Preadout(loop_i.fr))
            elif np.abs(loop_i.freq[-1]-loop_i.freq[0]) > 1.0e9:
                if Verbose == True:
                    print('Sweep bandwidth is {0} Hz. Sweep looks more like a survey. Preadout_dB is meaningless for a survey. Aborting Preadout computation... '.format(np.abs(loop_i.freq[-1]-loop_i.freq[0])))

            else:
                if Verbose == True:
                    print('No resonance frquency (fr) on record for selected resonance. Estimating fr using sweep minimum.')
                fr = np.extract(np.abs(loop_i.z).min() == np.abs(loop_i.z),loop_i.freq)[0]
                _define_sweep_array(Sweep_Array, index, Preadout_dB = Sweep_Array['Pinput_dB'][index] + fr)

        if Add_Temperatures == True:
            if metadata.Num_Temperatures < 1:
                condition = (Sweep_Array['Heater_Voltage'][index] + tol > Temperature_Calibration[:,0]) & (Sweep_Array['Heater_Voltage'][index] - tol < Temperature_Calibration[:,0])
                if condition.sum() >= 1:

                    Sweep_Array['Temperature'][index] = Temperature_Calibration[condition,1][0] # <-- Needs to be updated so that duplicate voltages are handled correctly
                else:
                    if Verbose == True:
                        print('Unable to match unique temperature to heater voltage value for Sweep_Array[{0}]. {1} matches found.'.format(index,condition.sum() ))
            else:
                _define_sweep_array(Sweep_Array, index, Temperature =     Sweep_Array['Temperature_Readings'][index].mean())
        # Clear out loop
        del(loop_i)

    if Verbose == True:
        print('\nSweep Array filled.')# Options selected Fit_Resonances = {0}, Compute_Preadout = {1}, Add_Temperatures = {2}'.format( Fit_Resonances,Compute_Preadout,Add_Temperatures))
