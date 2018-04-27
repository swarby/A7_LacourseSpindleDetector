function [ detectionVect ] = eventSmpList2DetectVect( NREMEventsLst, nSamples )
% Purpose: Generate a detection vector per sample from a list of events.
% Input 
%   NREMEventsLst : a matrix of event list [nEvents X 2]
%                   [start sample event 1, end sample event 1]
%                   [start sample event 2, end sample event 2]
%                   [                    ...                 ]
%   nSamples : number of samples in the detectionVect 
%               (usually from the number of samples in time series)
%
% Author: Karine Lacourse 2015-09-22
%
    detectionVect = zeros(nSamples,1);
    for iEvent = 1: size(NREMEventsLst,1)
        detectionVect(NREMEventsLst(iEvent,1):NREMEventsLst(iEvent,2))=1;
    end

end

