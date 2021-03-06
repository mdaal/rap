from .utils import  _extract_type, _define_sweep_array, _define_sweep_data_columns_legacy_gui

import numpy as np
import scipy.io
import datetime
import os
import warnings

# + Updated rap_demonstration.ipynb:  
# - to used importlib.reload instead of reload
# + Updated load_touchstone.py to  save  unicode strings for metadata.Time_Created, and metadata.Run instead of metadata.Run (which causes save_hf5() not to fail to parse the date created for the table name)
# +updated phase_fit to place commas in the display of the Q and f0
print("load_legacy_sweep_gui_data level init")
def load_legacy_sweep_gui_data(metadata, gui_data_path, which_copy = 'first'):
    data_dir = gui_data_path


    ### Find Latest sweep_config_xxxxx.mat and sweep_data_xxxxx.mat file in data_dir
    def get_latest_data_file(prefix, which_copy = which_copy ):
        '''
        Returns a list of files in the directory 'data_dir' which begin with 'prefix'.
        The contents of the list is determined by the  value of 'which_copy'
        which_copy= 'first'
                    'last'
                    'all'
                    'yyyymmddHHMMSS' (a specific date string)

        '''
        file_list = os.listdir(data_dir)
        modification_time = 0
        data_file = []
        found_file = False
        
        for file in file_list:
            if file.startswith(prefix):
                parts = file.replace('.', '_').split('_')

                data_file_failsafe = data_dir + os.sep + file #use this data_file is  last is actually data_file.mat
                
                if (len(parts) == 3) and (which_copy == 'first'): #data_file is data_file.mat
                    # only want this  when there is no date string, i.e. modification_time = 0
                    data_file.append(data_dir + os.sep + file)
                    found_file = True

                elif (len(parts) == 4) and (which_copy == parts[2]): #data_file is data_file_xxxxx.mat where xxxxx is the date as a string passed to which_copy
                    data_file.append(data_dir + os.sep + file)
                    found_file = True

                elif (len(parts) == 4) and (which_copy == 'last'): #data_file is data_file_xxxxx.mat where xxxxx is the latest date
                    print(which_copy, 'right here')
                    if int(parts[2]) >= modification_time:
                        modification_time = parts[2]
                        data_file.append(data_dir + os.sep + file)
                        found_file = True

                elif which_copy == 'all':
                        data_file.append(data_dir + os.sep + file)
                        found_file = True

        if found_file == False:
            data_file.append(data_file_failsafe)
            if (which_copy is not 'first') and (which_copy is not 'last'):
                warnings.warn('The specified copy of {0}, {0}_{1}.mat, was not found. Using {2}.'.format(prefix,which_copy,data_file),UserWarning)
        return  data_file

    data_dict= dict()
    data_dict['curr_config'] = dict()
    data_dict['measurement_metadata'] = dict()

    config_tags = [
            ('__header__', 'curr_config:Time_Created', str),
            ('curr_config','curr_config:savepath', str, 'savepath'),
            ('curr_config','curr_config:f0list', None, 'f0list'), #is a two row array
            ('curr_config','curr_config:starttemp', np.float, 'starttemp'), # is in mK
            ('curr_config','curr_config:stoptemp', np.float, 'stoptemp'), # is in mK
            ('curr_config','curr_config:steptemp', np.float, 'steptemp'), # is in mK
            ('curr_config','curr_config:startatten', np.float, 'startatten'),
            ('curr_config','curr_config:stopatten', np.float, 'stopatten'),
            ('curr_config','curr_config:stepatten', np.float, 'stepatten'),
            ('curr_config','curr_config:totpow', np.float, 'totpow'), #this is the total attenuation applied (input + output) in dB
            ('curr_config','curr_config:numsteps', np.float, 'numsteps'), #number of freq steps within 'span'
            ('curr_config','curr_config:samplestoavg', np.float, 'samplestoavg'), # number of samples to avg for IQ sweep
            ('curr_config','curr_config:span', np.float, 'span'), #the freq span, in MHz for fixed span scan
            ('curr_config','curr_config:adtime', np.float, 'adtime'), #individuation noise trace/sample duration
            ('curr_config','curr_config:tottime', np.float, 'tottime'), # total time to sample noise, i.e. 'adtime'*Num_Integrations
            ('curr_config','curr_config:donoise', int, 'donoise'), #is 1 if noise is to be measured, 0 otherwise
            ('curr_config','curr_config:sweeprate', np.float, 'sweeprate'), #IQ data sample rate
            ('curr_config','curr_config:noiserate', np.float, 'noiserate'), #noise data sample rate
            ('curr_config','curr_config:fixedspan', str, 'fixedspan'), #is the string u'Used fixed span' if using fixed span, as opposed to qlist
            ('curr_config','curr_config:qlist', None, 'qlist'), #is a two row array
            ('curr_config','curr_config:numfwhm', np.float, 'numfwhm'), #if using a qlist, this is the number of fwhm, determinining the freq span to scan
            ('curr_config','curr_config:quickfind', int, 'quickfind'), #is 1 to do a quick estimate of f0 and loop
            ('curr_config','curr_config:fullfit', int ,'fullfit'), #is 1 to do montecarlo for f0 and loop. NOT passed to doIQsweep
            ('curr_config','curr_config:saveraw', int ,'saveraw'), #is 1 to save time series in  binary files
            ('curr_config','curr_config:waitfortemp', int ,'waitfortemp'), # is 1 to wait for temp set point, and 1 for don't wait
            ('curr_config','curr_config:rawfreq', np.float ,'rawfreq'), #is in mK. If saving raw noise data, temp intervals at which noise traces should be saved.
            #e.g. setting to 20 saves noise data on the first step and every 20 mK after, provided raw_freq is a multiple of the step size
            ('curr_config','curr_config:spectra_check', int ,'spectra_check'), #is 1 to compute and save spectra for on and off resonance noise
            ('curr_config','curr_config:offres_check', int ,'offres_check'), #is 1 to take off res noise data
            ('curr_config','curr_config:decfac', np.float ,'decfac'), # is the decimation factor
            ('curr_config','curr_config:cpsd_box', int, 'cpsd_box'), # is 1 to calculate cross spectra between resonators
            ('curr_config','curr_config:heater_box', int, 'heater_box'), # is 1 to turn off heater
            ('curr_config','curr_config:dofilt', int, 'dofilt'), # is 1 to enable and use 100 kHz AA filter
            ('curr_config','curr_config:cdel', np.float ,'cdel'), # is the cable deley in nanoseconds
            ('curr_config','curr_config:spec_settings', None, 'spec_settings')] #is an two column array. col[0] is 'Noise_Frequency_Segmentation'; col[1] is 'Noise_Frequency_Resolution_Per_Segment'


    def _unpack_data_structure(tags, receiving_dict, struct_to_be_unpacked):
        ''' takes data specfied in tag from struct_to_be_unpacked and adds it to receiving_dict in the format specificed in tags
        the  placement of the data within receiving_dict is also specificed by tags, which is a list od tuples
        [(fieldname_in_struct_to_be_unpacked, destination_key_in_receiving_dict, format_in_receiving_dict, [optional  subfieldname_in_struct_to_be_unpacked]),...]
        '''
        for t in tags:
            # try:
                if t[1].find(':')>-1: #The case of a dictionary
                    t1 = t[1].split(':')

                    #This try/except block is for the case where metadata.__dict__['?'] is a dictionary
                    try:
                        receiving_dict[t1[0]].update([(t1[1],_extract_type(struct_to_be_unpacked[t[0]], return_type = t[2],field = t[3] if len(t) > 3 else None))])
                    except:
                        receiving_dict[t1[0]] = dict([(t1[1],_extract_type(struct_to_be_unpacked[t[0]], return_type = t[2],field = t[3] if len(t) > 3 else None))])
                else:
                    receiving_dict[t[1]] = _extract_type(struct_to_be_unpacked[t[0]], return_type = t[2],field = t[3] if len(t) > 3 else None)
            # except:
                #the case that the field does not exist or that its in an unexpected format
                #print('Field named {0}{1} is not found. Setting value to None'.format(t[0], (':'+t[3]) if len(t) > 3 else '')) # this usesScandata nomenclature
            #     print('Field named {0} is not found.'.format(t[1])) # this uses metadata nomenclature
    


    config_files = get_latest_data_file('sweep_config', which_copy) 
    sweep_data_files = get_latest_data_file('sweep_data', which_copy) #Contains datadata for Ben's fit code, i.e. zeropt and zerostd
    
    print('config_fileds and sweep_data_files', config_files,sweep_data_files)
    
    sweep_data  = scipy.io.loadmat(sweep_data_files[0])
    sweep_data_temp_list = list(np.floor(np.array(sweep_data['IQ_data'][0]['temprange'][0][0,:])*1000)/1000)
    sweep_data_atten_list = list(np.array(sweep_data['IQ_data'][0]['attenrange'][0][0,:],dtype = np.int)) #attenuator values must be ints
    print('sweep data temp list', sweep_data_temp_list)
    print('sweep data atten list', sweep_data_atten_list)

   
    config  = scipy.io.loadmat(config_files[0])

    _unpack_data_structure(config_tags, data_dict, config)


    #Now extract the creation time and Date in format : 'Wed Aug 09 17:15:14 2017'
    data_dict['curr_config']['Time_Created'] = data_dict['curr_config']['Time_Created'][data_dict['curr_config']['Time_Created'].find('Created on: ')+len('Created on: '):]
    if data_dict['curr_config']['Time_Created'][-1] == "'":
        data_dict['curr_config']['Time_Created'] = data_dict['curr_config']['Time_Created'][:-1]
    dt_start = datetime.datetime.strptime(data_dict['curr_config']['Time_Created'], '%a %b %d %H:%M:%S %Y')
    #Note: spec_settings qlist f0list are arrays

    ### get spec_setting into array format
    data_dict['measurement_metadata']['Noise_Frequency_Segmentation'] = np.array(np.squeeze(data_dict['curr_config']['spec_settings'][0,0]),dtype = np.float)
    data_dict['measurement_metadata']['Noise_Frequency_Resolution_Per_Segment'] = np.array(np.squeeze(data_dict['curr_config']['spec_settings'][0,1]),dtype = np.float)

    data_dict['measurement_metadata']['Synthesizer_Scan_Num_Points'] = int(data_dict['curr_config']['numsteps'])
    data_dict['measurement_metadata']['Measure_On_Res_Noise'] = bool(data_dict['curr_config']['donoise']) and  bool(data_dict['curr_config']['spectra_check'])
    data_dict['measurement_metadata']['Measure_Off_Res_Noise'] = bool(data_dict['curr_config']['donoise']) and bool(data_dict['curr_config']['spectra_check']) and bool(data_dict['curr_config']['offres_check'])
    data_dict['measurement_metadata']['Measure_On_Res_CPSD'] = bool(data_dict['curr_config']['donoise']) and bool(data_dict['curr_config']['spectra_check']) and bool(data_dict['curr_config']['cpsd_box'])

    ### Create a flag to let us know that this measurement_metadata was loaded from  Legacy Gui... used in _compute_noise_spectrum_length, in particular
    data_dict['measurement_metadata']['Is_Legacy_Gui_Data'] = True

    ### Construct list of output attenuator value settings
    output_atten_value_array = np.arange(data_dict['curr_config']['startatten'],data_dict['curr_config']['stopatten']+data_dict['curr_config']['stepatten'],data_dict['curr_config']['stepatten'])

    ### Construct list of temperature value settings
    print('curr_config start stop step',data_dict['curr_config']['starttemp'],data_dict['curr_config']['stoptemp']+data_dict['curr_config']['steptemp'],data_dict['curr_config']['steptemp'])
    temperature_value_array = np.arange(data_dict['curr_config']['starttemp'],data_dict['curr_config']['stoptemp']+data_dict['curr_config']['steptemp'],data_dict['curr_config']['steptemp'])

    ### Construct list of resonant frequency groups in the form [(group1res1, group1res2),(group2res1, group2res2), ...]
    resonator_group_list = [(np.float(data_dict['curr_config']['f0list'][0+c,0][0]), np.float(data_dict['curr_config']['f0list'][1+c,0][0]))  for c in range(0,data_dict['curr_config']['f0list'].shape[0],2) ]

    ### might want to sort the list of suffixes!
    spectra_filename_suffixes = ['{temp:.0f}-{resonator_group_num}-{start_atten:.0f}.mat'.format(temp = t, resonator_group_num = rg, start_atten = sa) for t in temperature_value_array for  rg in range(1,len(resonator_group_list)+1) for  sa in output_atten_value_array] # note that '{start_atten:.0f}' means don't include a decimal point
    #spectra_filename_suffixes = ['{temp}-{resonator_group_num}-{start_atten}.mat'.format(temp = str(t), resonator_group_num = rg, start_atten = int(sa)) for t in temperature_value_array for  rg in range(1,len(resonator_group_list)+1) for  sa in output_atten_value_array] # note that '{start_atten:.0f}' means don't include a decimal point
    spectra_filename_prefixes = ['spec', 'spec_offres'] if data_dict['measurement_metadata']['Measure_Off_Res_Noise'] else ['spec']
    print('curr_config temp value array', temperature_value_array)
    print('curr_config atten value array', output_atten_value_array)

    missing_spectra_filename = list()

    sweep_data_columns_list, sweep_data_columns = _define_sweep_data_columns_legacy_gui(data_dict['measurement_metadata'], fsteps_syn = data_dict['measurement_metadata']['Synthesizer_Scan_Num_Points'])

    number_of_sweeps = len(spectra_filename_suffixes) * 2 #times 2 becuase there are two resonators per group
    Sweep_Array = np.zeros(number_of_sweeps, dtype = sweep_data_columns)

    spectra_tags = [('__header__', 'Time_Created', str), ('spec1','spec1:wall', np.array, 'wall'), ('spec1','spec1:Piiall', np.array, 'Piiall'), ('spec1','spec1:Pqqall', np.array, 'Pqqall'),
            ('spec1','spec1:Piqall', np.array, 'Piqall'),('spec1','spec1:fn',np.float64, 'fn'), ('spec1','spec1:dZdf',np.float64, 'dZdf'),
            ('spec2','spec2:wall', np.array, 'wall'), ('spec2','spec2:Piiall', np.array, 'Piiall'), ('spec2','spec2:Pqqall', np.array, 'Pqqall'),
            ('spec2','spec2:Piqall', np.array, 'Piqall'),('spec2','spec2:fn',np.float64, 'fn'), ('spec2','spec2:dZdf',np.float64, 'dZdf'),
            ('specorth1','specorth1:wall', np.array, 'wall'), ('specorth1','specorth1:Piiall', np.array, 'Piiall'), ('specorth1','specorth1:Pqqall', np.array, 'Pqqall'),
            ('specorth1','specorth1:Piqall', np.array, 'Piqall'),('specorth1','specorth1:fn',np.float64, 'fn'), ('specorth1','specorth1:dZdf',np.float64, 'dZdf'),
            ('specorth2','specorth2:wall', np.array, 'wall'), ('specorth2','specorth2:Piiall', np.array, 'Piiall'), ('specorth2','specorth2:Pqqall', np.array, 'Pqqall'),
            ('specorth2','specorth2:Piqall', np.array, 'Piqall'),('specorth2','specorth2:fn',np.float64, 'fn'), ('specorth2','specorth2:dZdf',np.float64, 'dZdf'),
            ('traj1','traj1:f', np.array, 'f'),('traj1','traj1:z',np.array, 'z'), ('traj1','traj1:zr',np.complex128, 'zr'),
            ('traj2','traj2:f', np.array, 'f'),('traj2','traj2:z',np.array, 'z'), ('traj2','traj2:zr',np.complex128, 'zr'),
            ('cspec','cspec:Pi1i2all', np.array, 'Pi1i2all'), ('cspec','cspec:Pq1q2all', np.array, 'Pq1q2all'), ('cspec','cspec:Pi1q2all', np.array, 'Pi1q2all'),
            ('cspec','cspec:Pq1i2all', np.array, 'Pq1i2all')] #note  spec1:fn = specorth1:fn  and spec1:dZdf = specorth1:dZdf

    ### remove CPSD related tags if CPSD is not measured
    if data_dict['measurement_metadata']['Measure_On_Res_CPSD'] is False:
        spectra_tags = [tag for tag in spectra_tags  if tag[0] is not  'cspec']

        
    i=0
    on_res = dict()
    off_res = dict()
    dt_duration = datetime.timedelta()
    unpack_dict = {'spec': (spectra_tags, on_res), 'spec_offres': ([tag for tag in spectra_tags  if tag[0] is not  'cspec' if tag[0] is not 'traj1' if tag[0] is not 'traj2'], off_res)}
    
    ### loop through specXX-YY-ZZ.mat and spec_offresXX-YY-ZZ.mat files, pullling their data and filling Sweep_Array
    print('spectra_filename_suffixes',spectra_filename_suffixes)
    for filename_suffix in spectra_filename_suffixes:
        for filename_prefix in spectra_filename_prefixes:
            spectra_filename = data_dir + os.sep + filename_prefix + filename_suffix

            if os.path.isfile(spectra_filename):
                spectra  = scipy.io.loadmat(spectra_filename)
                _unpack_data_structure(unpack_dict[filename_prefix][0],unpack_dict[filename_prefix][1],  spectra)
            else:
                missing_spectra_filename.append(filename_prefix + filename_suffix)
                continue
        if filename_prefix + filename_suffix in missing_spectra_filename:
            continue
        temp_group_atten_list = filename_suffix.replace('.mat', '').split('-') #note: temp_group_atten_list is a list of strings

        ### get time/date the on_res file was created and then store it as a datetime object.
        on_res['Time_Created'] = on_res['Time_Created'][on_res['Time_Created'].find('Created on: ')+len('Created on: '):]
        if on_res['Time_Created'][-1] == "'":
            on_res['Time_Created'] = on_res['Time_Created'][:-1]
        dt_time_created = datetime.datetime.strptime(on_res['Time_Created'], '%a %b %d %H:%M:%S %Y')

        ### find max time interval between on_res file creation and measurement start time, to be used to compute total measurement duration
        dt_duration = dt_time_created - dt_start if dt_time_created - dt_start > dt_duration else dt_duration
        print('temp_group_atten_list',temp_group_atten_list)
        try: 
            tmp = sweep_data_temp_list.index(np.float(temp_group_atten_list[0]))
        except: 
            print('tmp not found in datafiles')

        try: 
            grp = np.int(temp_group_atten_list[1])*2 + 1
        except: 
            print('grp not found in datafiles')

        try:
            tmp = sweep_data_temp_list.index(np.float(temp_group_atten_list[0]))
            grp = np.int(temp_group_atten_list[1])*2 + 1
            atn = sweep_data_atten_list.index(np.int(temp_group_atten_list[2]))  #attenuator values must be ints
            print(tmp,atn,grp)
        except:
            print('atm not found in datafiles')
        #sweep_data['IQ_data'][0]['temps'][0][0,tmp]['attens'][0,1]['res'][0,1]['zeropt'][0,0] 
        
        _define_sweep_array(Sweep_Array, i,
                                Temperature = np.float(temp_group_atten_list[0]) * 0.001,
                                S21_Syn = on_res['traj1']['z'].flatten(),
                                Frequencies_Syn = on_res['traj1']['f'].flatten() * np.power(10.,9),
                                Noise_Freq_On_Res = on_res['spec1']['fn'] * np.power(10.,9),
                                Noise_S21_On_Res = on_res['traj1']['zr'],
                                Noise_II_On_Res = on_res['spec1']['Piiall'].flatten(),
                                Noise_QQ_On_Res = on_res['spec1']['Pqqall'].flatten(),
                                Noise_IQ_On_Res = on_res['spec1']['Piqall'].flatten(),
                                Noise_Freq_Off_Res = off_res['spec1']['fn'] * np.power(10.,9) if data_dict['measurement_metadata']['Measure_Off_Res_Noise'] else 0,
                                Noise_II_Off_Res =  off_res['spec1']['Piiall'].flatten() if data_dict['measurement_metadata']['Measure_Off_Res_Noise'] else 0,
                                Noise_QQ_Off_Res = off_res['spec1']['Pqqall'].flatten() if data_dict['measurement_metadata']['Measure_Off_Res_Noise'] else 0,
                                Noise_IQ_Off_Res = off_res['spec1']['Piqall'].flatten() if data_dict['measurement_metadata']['Measure_Off_Res_Noise'] else 0,
                                Noise_S21_Off_Res = on_res['traj1']['zr'] if data_dict['measurement_metadata']['Measure_Off_Res_Noise'] else 0,
                                Noise_Freq_Vector = on_res['spec1']['wall'].flatten(),
                                Noise_Chan_Input_Atn = np.float(temp_group_atten_list[2]), #atten 1 from matlab code, 'input' side of fridge
                                Noise_Chan_Output_Atn = np.max([data_dict['curr_config']['totpow'] - np.float(temp_group_atten_list[2]),0]), #CANNOT BE NEGATIVE!! atten 2 from matlab code
                                Scan_Timestamp = dt_time_created.strftime('%Y%m%d%H%M'),
                                Resonator_Index = 2*np.float(temp_group_atten_list[1]) - 2,
                                Resonator_Group = np.array([np.float(temp_group_atten_list[1]),i,i+1]),
                                dZdf = on_res['spec1']['dZdf']
                                )
        i = i + 1
        _define_sweep_array(Sweep_Array, i,
                                Temperature = np.float(temp_group_atten_list[0]) * 0.001,
                                S21_Syn = on_res['traj2']['z'].flatten(),
                                Frequencies_Syn = on_res['traj2']['f'].flatten() * np.power(10.,9),
                                Noise_Freq_On_Res = on_res['spec2']['fn'] * np.power(10.,9),
                                Noise_S21_On_Res = on_res['traj2']['zr'],
                                Noise_II_On_Res = on_res['spec2']['Piiall'].flatten(),
                                Noise_QQ_On_Res = on_res['spec2']['Pqqall'].flatten(),
                                Noise_IQ_On_Res = on_res['spec2']['Piqall'].flatten(),
                                Noise_Freq_Off_Res = off_res['spec2']['fn'] * np.power(10.,9) if data_dict['measurement_metadata']['Measure_Off_Res_Noise'] else 0,
                                Noise_II_Off_Res =  off_res['spec2']['Piiall'].flatten() if data_dict['measurement_metadata']['Measure_Off_Res_Noise'] else 0,
                                Noise_QQ_Off_Res = off_res['spec2']['Pqqall'].flatten() if data_dict['measurement_metadata']['Measure_Off_Res_Noise'] else 0,
                                Noise_IQ_Off_Res = off_res['spec2']['Piqall'].flatten() if data_dict['measurement_metadata']['Measure_Off_Res_Noise'] else 0,
                                Noise_S21_Off_Res = on_res['traj2']['zr'] if data_dict['measurement_metadata']['Measure_Off_Res_Noise'] else 0,
                                Noise_Freq_Vector = on_res['spec2']['wall'].flatten(),
                                Noise_Chan_Input_Atn = np.float(temp_group_atten_list[2]),
                                Noise_Chan_Output_Atn = np.max([data_dict['curr_config']['totpow'] - np.float(temp_group_atten_list[2]),0]),
                                Scan_Timestamp = dt_time_created.strftime('%Y%m%d%H%M'),
                                Resonator_Index = 2*np.float(temp_group_atten_list[1])-1,
                                Resonator_Group = np.array([np.float(temp_group_atten_list[1]),i-1,i]),
                                dZdf = on_res['spec2']['dZdf']
                                )
        i = i + 1
        on_res.clear()
        off_res.clear()




    ### New measurement_metadata keys: 'Save_Time_Series','Is_Legacy_Gui_Data','Measure_On_Res_CPSD'
    metadata.Data_Source = data_dict['curr_config']['savepath']
    if metadata.Run is None:
        metadata.Run = metadata.Data_Source.replace(':','').replace('\\','_').replace('//','_')
    metadata.Meaurement_Duration = dt_duration.total_seconds()
    metadata.Wait_Time = data_dict['curr_config']['waitfortemp']
    metadata.Num_Temp_Set_Points = temperature_value_array.size
    metadata.Time_Created = dt_start.strftime('%Y%m%d%H%M')
    metadata.Electrical_Delay = data_dict['curr_config']['cdel'] * np.power(10.,-9) if data_dict['curr_config']['cdel'] is not None else None
    metadata.Num_Heater_Voltages = len(temperature_value_array)
    metadata.Num_Powers = len(output_atten_value_array)
    metadata.Loop_Data_Column_Names  = ("Frequencies_Syn","S21_Syn")
    # metadata.Digitizer = 'NI6120'
    data_dict['measurement_metadata']['IQ_Sample_Rate'] = data_dict['curr_config']['sweeprate']
    data_dict['measurement_metadata']['Noise_Sample_Rate'] = data_dict['curr_config']['noiserate']
    data_dict['measurement_metadata']['Noise_Hrdwr_AA_Filt_In_Use'] = bool(data_dict['curr_config']['dofilt'])
    data_dict['measurement_metadata']['IQ_Num_Samples']  = data_dict['curr_config']['samplestoavg']
    data_dict['measurement_metadata']['Synthesizer_Scan_Num_BW'] = data_dict['curr_config']['numfwhm']
    data_dict['measurement_metadata']['Noise_Decimation_Factor'] = data_dict['curr_config']['decfac']
    data_dict['measurement_metadata']['Noise_Integration_Time'] = data_dict['curr_config']['adtime']
    data_dict['measurement_metadata']['Noise_Num_Integrations'] = data_dict['curr_config']['tottime']/data_dict['curr_config']['adtime']
    data_dict['measurement_metadata']['Save_Time_Series'] = bool(data_dict['curr_config']['saveraw'])

    if len(missing_spectra_filename) > 0:
        # print('The following datafiles are expected but missing from the directory:')
        # print(missing_spectra_filename)
        warnings.warn('The following datafiles are expected but missing from the directory: {}'.format(missing_spectra_filename),UserWarning)

    return sweep_data_columns_list, sweep_data_columns, Sweep_Array
    # data_dict['curr_config']['heater_box'] <-- The 'turn off heater box' option
    # data_dict['curr_config']['quickfind'] and  data_dict['curr_config']['fullfit'] <-- choice of how to fit resonator...'quick fit' versus full fit options on gui
