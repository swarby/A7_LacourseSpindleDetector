% event_StartsEndsDurations.m
% 
% Purpose:
% Identification of start and stop indexes of events (1) in a binary (logical 0/1) vector.
% Also outputs the duration of the events, ie the number of consecutive events.
%
% Usage:
% [STARTS, ENDS, DURATIONS] = event_StartsEndsDurations(dataVector)
%
% dataVector - input timeseries, needs to be binary or logical (0/1).
%
% STARTS - index of the start positions of events
% ENDS - index of the stop position
% DURATIONS - number of consecutive events
%
% Requirements:
%
%
% Notes:
% The unit of duration is the 'number of events', ie usually samples.
%
% Examples:
% [mx(:,1), mx(:,2), mx(:,3)] = event_StartsEndsDurations(mockdata)
%
% 
% Authors:
% Julie Christensen 2014-06-15
% modified by Simon Warby 2014-05-05
% 
% Changelog:
% Changed function name; added duration to output SCW 2014-06-16


function [STARTS, ENDS, DURATIONS] = event_StartsEndsDurations(dataVector)

%% Input checks
if nargin==1
    if length(find(dataVector==0|dataVector==1))<length(dataVector)
        error(['Only binary vectors allowed for single input to '...
            '.']);
    end  
    if ~isvector(dataVector)
            error('Only vectors allowed as input.');
    end
else
    if any(~isvector(dataVector))
        error('Only vectors allowed as input.');
    end
end


%% Converting to row vectors to ensure similar inputs
if ~isrow(dataVector)
    dataVector = dataVector';
end


%% Finding edges
cross_edges = dataVector(1:end-1)-dataVector(2:end);
STARTS = find(cross_edges==-1)+1;
ENDS = find(cross_edges==1);

if dataVector(1)==1
    STARTS = [1 STARTS];
end
if dataVector(end)==1
    ENDS = [ENDS length(dataVector)];
end

DURATIONS = (ENDS+1) - STARTS ; 

STARTS = STARTS' ;
ENDS = ENDS' ;
DURATIONS = DURATIONS' ; 

% Bind into a matrix
% x = [STARTS, ENDS, DURATIONS]
% Bind directly into a 3 colum matrix
% [mx(:,1), mx(:,2), mx(:,3)] = event_StartsEndsDurations(mockdata)

end   


%% Testing
%   mockdata = [0 0 1 0 0 0 1 0 0 0 1 1 1 1 1 0 0 0 0 ]
%   mockdata = [1 0 1 0 0 0 1 0 0 0 1 1 1 1 1 0 1 1 0 1]
%   [start_i, end_i, duration_i] = event_StartsEndsDurations(mockdata)
%   x = [start_i, end_i, duration_i]

