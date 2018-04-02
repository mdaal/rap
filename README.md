# rap
Resonator Analysis Package


To Do
+ Add comma to Qi, Qc , Q in phase fit plot to make it easier to read
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




Resonance Fitting Methods:
- 'Diameter Correction Method' of Khalil et al. 2012  https://doi.org/10.1063/1.3692073
- Fitting error estimation of Probst et al. 2015 https://doi.org/10.1063/1.4907935
- Circle fit methods: ChernovHoussam, ChernovLesort, LevenbergMarquardtFull, LevenbergMarquardtReduced (wrap cpp code in python); some have signal to noise estimate code
- 'CPZM' method of Deng et al. https://doi.org/10.1063/1.4817512  (linear fractional transformation where coefficients also give circuit element values)
- Signal to noise est, 3dB method, Lorentzian Fit method, resonance curve area  method, inverse mapping technique, Snortland* method, (radial weighting) of Petersan et al. https://doi.org/10.1063/1.368498
