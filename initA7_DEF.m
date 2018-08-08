% initA7_DEF initiates all basic parameters required to run a7 detectors
%  Purpose :
%       paths, input and output definition
%       Turn On/Off context classifier
%       Select appropriate sampling rate
%       Select sleep staging baseline 
%       (can be one, several or no sleep stages)
%  Input :
%       DEF_a7: structure that contain important a7 features
%  Output :
%       DEF_a7: structure with all the parameters set as defined below
%
%
% Author: Jacques Delfrate
% Date : 2018-02-13

clear variables;
clc;

 % Path
    % Define inputs directory
    DEF_a7.inputPath            = './input/';                % input path directory
    % Define inputs files in the input directory
    DEF_a7.EEGvector            = 'EEGvector.mat';           % input eeg vector name
    DEF_a7.sleepStaging         = 'sleepStaging.mat';        % input sleep stages vector name
    DEF_a7.artifactVector       = 'artifactVector.mat';      % input artifact vector name
    % Define outputs
    DEF_a7.outputGrandParentDir = './output/';               % output path
    DEF_a7.outputTxtFile        = 'EventDetection.txt';      % Summury of detection
    DEF_a7.outputTS             = 'detectionInfo.mat';       % detection info output file name
    DEF_a7.outputDetectInfo     = 'detectionVector.mat';     % TS output file name
    DEF_a7.outputNREMClass      = 'nremClassifier.mat';      % context classifier info output file name
    DEF_a7.detectorDef          = 'DEF_a7.mat';              % detector info output file name
    
 % A7 default parameters
    % Sampling rate of the eeg vector.
    DEF_a7.standard_sampleRate = 100;
    % Sleep stage considered in the baseline.  
    % Set to ['N2', 'N3'] to consider stage 'N2' and 'N3'  as baseline, 
    % ['2'] for just stage '2' etc.
    % Note that the stage units must match what is found in the
    % sleepStaging input file (sleepStaging.mat) 
    % (ie use '2', or 'N2' if '2' or 'N2' is used etc.)
    DEF_a7.bslSleepStaging = ['N2'];

