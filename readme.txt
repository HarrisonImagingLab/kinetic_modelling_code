MATLAB code for generating kinetic modelling parameters from TACs and blood data. 

blood data should be in the form PMOD accepts, with time units in MINUTES. 

TAC data should be with time in Seconds as exported by PMOD. The first column of the TAC should be the whole brain -- it will be used to fit the blood delay (if that is chosen) and give intiial parameter guesses for all other regions. Figures of plasma fits and the first 6 region fits will be produced and saved in the current working directory. 

All files should be .csv or .txt and tab delimited. 

Results for outcome parameters in 2TCM (V_T) align very well with PMOD and fits are qualitatively the same. Weighting for 2TCM fitting is done by whole brain.

For Logan, V_T is identical to that obtained by PKIN of PMOD. Logan loads the fitted plasma activity curve from the 2TCM modelling, which is faster than integrating a function handle. 

