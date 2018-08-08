function [detInfoTS, outputFile] = organiseOutputTxtFile(detVect, NREMClass, ...
    sleepStageVect, DEF_a7, slowRInfoTS, detInfoTS)

% Purpose
%   Organise start, end and duration of events detected to be saved into a
%   text file
% Input
%   detVect        : binary vector of detected events in sample
%   NREMClass      : binary vector of context classifier
%   sleepStageVect : sleep stage vector
%   DEF_a7         : definition parameter of the current detector
%   slowRInfoTS    : Time series with slow ratio (context classifier)
%                    information
% Output
%   detInfoTS      : matrix of detection info converted into a time series 
%                    (with the same size than timeSeries)
%                    first 4 rows are the features used to detect spindles
%                    1 : PSDSigmaLog (sigma power log10 transformed)
%                    2 : relSigPow (z-score of relative sigma power from 30 sec clean around the PSA window)
%                    3 : sigmaCov (z-score of the covariance of sigma from 30 sec clean closed to the PSA window)
%                    4 : sigmaCorr (correlation of sigma)
%                    5 : alphaLog (alpha power log10 transformed)
%                    6 : PSDLowFreq (delta + theta power for slow ratio)
%                    7 : PSDHighFreq (beta power for slow ratio)
%                    8 : artifact vector (logical vector 0: no artifact, 1: artifact) 
%   outputTxtFile  : outputMatrix organized to be saved
% 
% Requirement
%   event_StartsEndsDuration.m
%
% Author : Jacques Delfrate 2018-03-07
% Change Log:
%        - 

% Computes start, end and duration of detected events
[starts, ends, durations] = event_StartsEndsDurations(detVect);

% Organise output matrix
startSec    = round(starts/DEF_a7.standard_sampleRate, 1);
endSec      = ends/DEF_a7.standard_sampleRate;
durationSec = durations/DEF_a7.standard_sampleRate;

% get sleep staging for each spindle
sleepStage = sleepStageVect(starts, 1);

if ~isempty(slowRInfoTS)
    % slow ratio output    
    detInfoTS(6,:) = slowRInfoTS;
    detInfoTS      = detInfoTS';
    NREMClass      = cell2mat(NREMClass);

    % Column label
    titleLabel = {'start_sample', 'end_sample', 'duration_sample', ...
        'start_sec', 'end_sec', 'duration_sec', ...
        'contextClassifier', 'sleepStage'};

    outputFile = num2cell([starts, ends, durations, startSec, endSec, ...
        durationSec, NREMClass]);
    sleepStage = cellstr(sleepStage);
    outputFile = [outputFile, sleepStage];
    outputFile = [titleLabel; outputFile];   
    
else 
    % get sleep staging for each spindle
    sleepStage = sleepStageVect(starts, 1);

    % Column label
    titleLabel = {'start_sample', 'end_sample', 'duration_sample', ...
        'start_sec', 'end_sec', 'duration_sec', 'sleepStage'};

    outputFile = num2cell([starts, ends, durations, startSec, endSec, ...
        durationSec]);
    sleepStage = cellstr(sleepStage);
    outputFile = [outputFile, sleepStage];
    outputFile = [titleLabel; outputFile];    
    
end
end

