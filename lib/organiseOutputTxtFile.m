function outputFile = organiseOutputTxtFile(detVect, NREMClass, ...
    sleepStageVect, DEF_a7)

% Purpose
%   Organise start, end and duration of events detected to be saved into a
%   text file
% Input
%   detVect   : binary vector of detected events in sample
%   NREMClass : binary vector of context classifier
% Output
%   outputTxtFile : outputMatrix organized to be saved
% 
% Requirement
%   event_StartsEndsDuration.m
%
% Author : Jacques Delfrate 2018-03-07

% Computes start, end and duration of detected events
[starts, ends, durations] = event_StartsEndsDurations(detVect);

% get sleep staging for each spindle
sleepStage = sleepStageVect(starts, 1);

% Column label
titleLabel = {'start_sample', 'end_sample', 'duration_sample', ...
    'start_sec', 'end_sec', 'duration_sec', ...
    'contextClassifier', 'sleepStage'};

% Organise output matrix
startSec    = round(starts/DEF_a7.standard_sampleRate, 1);
endSec      = ends/DEF_a7.standard_sampleRate;
durationSec = durations/DEF_a7.standard_sampleRate;

outputFile = [starts, ends, durations, startSec, endSec, ...
    durationSec, NREMClass, sleepStage];
outputFile = [titleLabel; num2cell(outputFile)];
end

