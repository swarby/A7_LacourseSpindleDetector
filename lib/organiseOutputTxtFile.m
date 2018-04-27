function outputFile = organiseOutputTxtFile(detVect, NREMClass, ...
    sleepStageVect)

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
titleLabel = {'start', 'end', 'duration', 'contextClassifier', ...
    'sleepStage'};

% Organise output matrix
outputFile = [starts, ends, durations, NREMClass, sleepStage];
outputFile = [titleLabel; num2cell(outputFile)];
end
