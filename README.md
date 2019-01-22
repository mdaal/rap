# rap
Resonator Analysis Package


To Do
+ Add name of dataset on top of phase fit report
+ Add option to specify frequency interval in trim_loop and/or center freq so that multiple loops in a sweep can be ignored. 
+ Enable Python 3.6 compatability
+ Implement Pulse analysis code
+ Fully Implement Sweep_Array["Resonator_Index"]
+ Phase out Sweep_Array["Fstart"], Sweep_Array["Fstop"],
+ Phase out metadata.Fsteps, metadata.Num_Ranges (replace with num_index?)
+ replace circle_fit python code with c-code
+ Implement save to Legacy GUI format
+ Implement save/load of time series for noise
+ Implement cross power spectral density data structure
+ Implement more fit options
+ Implement fit uncertainty/error 
+ Setup logging and exception handling
+ Make Sweep_Array an object to accomodate big data sets
+ Change phase fit success from the method that gave the minimum function value to the method that gave minimum function value  AND was a successful fit.
+ Add a calibration folder on the top level (with sweeps, and pulses folder)... it might be a class...  
+ Store calibrations as functions in a dict in the calibration folder. each object references a function to apply its calibration
+ Calibrations to include alplifiers, attenuattors, cables, mixers etc. Should include a temperature/noise temperature declaration. 
+ Calibrations can also include device/resonator specific information like impedance, widths, and Eeff, etc.. 
+ calibrations will be have  sequential structure that captures the order in which signals propagate through read out  system elements (e.g. first through LNA then RT amp) 
+ analysies will be in the form of piplines. not sure where to put this code yet. 
+ piplines will also have a sequencial structure which captures the  order  that the  data is processed (eg trim loop then  circle fit then phase fit), and can incorporate 3rd party analysis packages such as 'scraps'



Resonance Fitting Methods:
- 'Diameter Correction Method' of Khalil et al. 2012  https://doi.org/10.1063/1.3692073
- Fitting error estimation of Probst et al. 2015 https://doi.org/10.1063/1.4907935
- Circle fit methods: ChernovHoussam, ChernovLesort, LevenbergMarquardtFull, LevenbergMarquardtReduced (wrap cpp code in python); some have signal to noise estimate code
- 'CPZM' method of Deng et al. https://doi.org/10.1063/1.4817512  (linear fractional transformation where coefficients also give circuit element values)
- Signal to noise est, 3dB method, Lorentzian Fit method, resonance curve area  method, inverse mapping technique, Snortland* method, (radial weighting) of Petersan et al. https://doi.org/10.1063/1.368498
