# README #

### What is this code for? ###

This package is the A7 spindle detector algorithm that detects sleep spindles.  The A7 matlab function runs on one EEG channel at a time, but can be run on multiple EEG channels in parallel/serial for multichannel recordings.  All of the initial testing of the detector has been performed on C3-A2 in stage N2 sleep.

Three vectors loaded into matlab are required for the A7 function: EEG timeseries vector, sleep stages vector and EEG signal artifact vector.  For convenience, we have provided an example input dataset that loads the three vectors from binary .mat files and exports the spindle detections to an output file.

* Release Version 1.0
The package has been developed on Linux 14.04 with MATLAB 9.1.0.441655 (R2016b).
Note that this code may not run on previous versions of MATLAB.  For example, it does not work with R2012a due to changes to 'omitnan' flags for math functions (ie sum).

### What will I find in this code package? ###

This package is composed of a main script a7MainScript.m, a setting inita7_DEF.m function, an input, output and a library folder. The core functions of the A7 algorithm are found in the a7SpindleDetection.m file.

* Input folder
  * This is where example input vectors are found. We have provided this example data so you can inspect the format of the expected input vectors.
  * Example input : 
    * EEGvector.mat : 1 minute of broadband filtered EEG recording, sample-by-sample
    * sleepStaging.mat : sleep stages that match the EEG timeseries, sample-by-sample
                         For example this vector can be coded as: 0 or wake, 1 or N1, 2 or N2, 3 or N3 and 5 or REM etc.
			 Multiple sleep stages can be used as baseline by modifying the DEF_a7.bslSleepStaging variable. Variables in baseline should match the sleep stage vector.
    * artifactVector.mat : binary vector with/without (1/0) artifact that matches sample-by-sample the EEG timeseries.

* Output folder 
  * This is where output data should be found after running the code. For more information about output data, please look at the output_example.mat data files.
  * output_example
    * detectionVector.mat : binary vector with/without (1/0) spindle detected, in sample
    * detectionInfo.mat : matrix with:  
					1 - log10(absoluteSigmaPower), by-sample
					2 - relativeSigmaPower as a z-score, by-sample
					3 - sigmaCovariance as a z-score, by-sample
					4 - sigmaCorrelation, by-sample
					5 - artifact detection vector, by-sample
					6 - log10(slowRatioVector), by-sample
    * nremClassifier.mat : by-event binary vector that indicates if a detected spindle is In or Out of the spectral context defined by the slowRatio threshold.
    * EventDetection.txt : tab-delimited text file listing by-event detections. Each row is a spindle event.  The output has four columns: 
					1 - start of the event (sample number)
					2 - end of the event (sample number)
					3 - duration of the event in (number of samples)
					4 - spectral context classification from the NREMClassifier (binary column only if context classifier is On, no column if context classifier is Off)
					5 - sleep staging					

* Library folder
  * This folder contains all the matlab function required to run this package

** Notes **
As an example we choose a 1 min artifact free recording from N2 sleep (EEGVector.mat).


### How do I get set up? ###

To run the example data:
* Open initA7_DEF.m and modify the settings available to match your configuration.
* Run the a7MainScript.m script in matlab.
* Go to the output folder to view the results.

To run your data, you have two options:
1. Replace the input files with your data, in the same .mat format.
2. Modify the a7MainScript.m script to point to your input vectors already in matlab, instead of loading them directly from file.  In this case, you would not load external matlab files (ignore section 1), but point the a7SpindleDetection function (section 2) to the correct vectors already in matlab.

### Interpreting the Results ###
The result of the spindle detection is found in EventDetection.txt, where each row is a spindle event. Note that the output contains all spindle events found in non-artifact data, regardless of the context classification and sleep stage.

It is likely that you will want to filter the spindle detections in the results file such as:
-Remove spindle detections that are too short or too long for the definition you are using.  Note that the default duration in the A7 detector is 0.3s to 2.5s.
-Remove spindle detections that are "OUT" of the normal spindle context (contextClassifier=0).  Alternatively, you may choose to keep all spindle detections ("IN"=1 and "OUT"=0) or analyze "OUT" of context spindles separately, depending on your analysis plan.  See the Lacourse et al. manuscript for details about the context classifier.
-Filter for spindles that are in stage 2 or NREM (etc) depending on your analysis plan.


### Contributors ###

* The code was written by Karine Lacourse (M.Ing.)
* The code was reviewed by Jacques Delfrate (M.Sc.A.)
* The project was under the supervision of Dr. Simon Warby (Ph.D.)

### Who do I talk to if I need help? ###

* Contact Simon Warby : simon.c.warby@umontreal.ca

### REMARKS : ###
    Free use and modification of this code is permitted, provided that
    any modifications are also freely distributed.

    When using this code or modifications of this code, please cite:
        Lacourse K, Defrate J, Beaudry J, Peppard P, Warby SC. 2018. A sleep spindle 
        detection algorithm that emulates human expert spindle scoring.  [Full citation 
        to be determined].
