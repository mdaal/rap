class metadata:
	'''Every data set (e.g. as survey, power sweep, temp sweep) is stored as a pytable. 
	Each pytable has an metadata instance associated with it. 
	This specifies the contents of the metadata.
	'''
	def __init__(self):
		#metadata imported from scan data
		self.Time_Created = None
		self.Atten_Added_At_NA = None # redundant if self.Atten_NA_Input and self.Atten_NA_Output are defined; should be merged somehow
		self.NA_Average_Factor = None
		self.Fridge_Base_Temp = None
		self.Box = None
		self.Ground_Plane = None
		self.Ground_Plane_Thickness = None
		self.LNA = None
		self.IFBW = None
		self.Test_Location = None
		self.Minimum_Q = None
		self.Notes = None
		self.Num_Points_Per_Scan = None
		self.Wait_Time = None # in seconds
		self.Press = None
		self.Min_Freq_Resolution = None
		self.Run = None
		self.Sensor = None
		self.Fridge_Run_Start_Date = None
		self.Fsteps  = None #### Might want to get rid of this
		#self.IS_Sonnet_Simulation = None
		self.Data_Source = None
		self.Atten_At_4K = None
		self.Atten_NA_Output = None # positive value in dB
		self.Atten_NA_Input = None # positive value in dB
		self.Atten_RTAmp_Input = None # positive value in dB
		self.RTAmp_In_Use = None
		self.Meaurement_Duration = None # in seconds
		self.Num_Heater_Voltages = None
		self.Num_Powers = None
		self.Num_Ranges = None 	#### Might want to get rid of this	
		self.Num_Temperatures = None #number of temperature readings taken after every scan for each heater voltage/power
		self.Thermometer_Configuration = None
		self.Thermometer_Voltage_Bias = None # for UC Berkeley readout. 

		self.Num_Temp_Set_Points = None #For systems where the fridge temp is controlled by PID. number of temperature set points during the measurement
		self.Plot_Preferences =None #A list of plot preferences for the sweep

		# manual entry metadata
		self.Electrical_Delay = None # Seconds ---  computed at time of data library generation
		self.Resonator_Width = None #list if more than one
		self.Resonator_Thickness = None #list if more than one
		self.Resonator_Impedance = None
		self.Resonator_Eeff = None # Resonator Dielectric Constant
		self.Feedline_Impedance = None
		self.Cable_Calibration = None
		self.Temperature_Calibration = None # a list of tuples [(heatervoltge1, temperature), (heatervoltage2,temperature, ...)]
		self.System_Calibration = None

		self.RTAmp = None
		self.Digitizer = None


