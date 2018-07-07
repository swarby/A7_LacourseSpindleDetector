### README ###

This summarizes each MATLAB function implemented for the A7 spindle detection algorithm.

# a7 folder main scripts #

* a7MainScript.m
	Main script to run to detect spindles 
* initA7_DEF.m
	Initiate all basic parameters required to run a7 detectors


# lib folder functions summary #

* a7SpindleDetection.m
	Compute features related to sigma power to detect spindles.
* a7subAbsPowValues.m
	Compute the absolute sigma power for different frequency band through the average sample.
* a7subDetWithIntv4featAnd.m
	Detect spindles based on the 4 features
* a7subRelSigPow.m
	Detect spindles based on the relative sigma power. 
* a7subSigmaCov.m
	Detect spindles based on the covariance between the raw data and the data filtered in the sigma band.
* a7subSigmaCorr.m
	Detect spindles based on the correlation between the raw data and the data filtered in the sigma band.
* a7subTurnOffDetSlowRatio.m
	Turn off detections when they occur in an unexpected context which can be estimated through a spectral profile.
* butterFiltZPHighPassFiltFilt.m
	High pass filter.
* butterFiltZPLowPassFiltFilt.m
	Low pass filter.
* calcPSDNormType.m
	Calculate the Power Spectral Density (PSD) of a time signal.
* cell2tab.m
    Save data into text file
* eventSmpList2DetectVect.m
	Generate a detection vector per sample from a list of events.
* event_StartsEndsDurations.m
	Identification of start and stop indexes of events (1) in a binary (logical 0/1) vector.
* findFreqIndex.m
	To find the nearest index in a frequency vector that corresponds to a frequency band.
* organiseOutputTxtFile.m
	Organise start, end and duration of events detected to be saved into a text file
* percentiles.m
	Returns percentiles of the values.
* samples2WindowsInSec.m
	Shape the data per samples into data per windows.
* windows2SamplesInSec.m
	Convert a data per window in a data per sample.


