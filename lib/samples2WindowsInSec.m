function windowsPerSamples = samples2WindowsInSec(...
    dataPerSamples, nWindows, winLengthInSec, winStepInSec, samplingRate)
% Purpose:
% Shape the data per samples into data per windows. The input could be a
% tall vector as a time series and the output would be a matrix of the 
% samples included in each windows [nWindows x winLengthInSample]. 
% The output information is redundant if there is overlap between windows
% because the same sample would be included in more than one window.
%
% Inputs:
%   dataPerSamples:     Data vector per samples ex) time series 
%                       (can be wide or tall)
%   nWindows:           (int) The number of windows to generate
%   winLengthInSec:     (double) Number of sec included in a window
%   winStepInSec:       (double) Number of sec to step for the following window
%                       (a new window is considered at each
%                       "winStepInSec" sec)
%   samplingRate:       (double) Sampling rate in Hz
%
% Outputs:
%   windowsPerSamples:  Matrix of [nWindows, round(winLengthInSec * samplingRate)]
%                       A Matrix of every samples shape to work with window
%                       NaNs fill incomplete windows.
%
% Authors : Karine Lacourse 
% Date    : 2016-10-24
% 
%-------------------------------------------------------------------------
    
    % These variables need to be int
    nWindows = round(nWindows);
    winLengthInSample = round(winLengthInSec * samplingRate);

    windowsPerSamples = nan(nWindows, winLengthInSample);
    % Convert sample event to window event
    for iWind = 1: nWindows
        startIndex = round((iWind-1) * winStepInSec * samplingRate);
        % The number of windows asked could be too high
        % Make sure the index asked is not outside dataPerSamples
        if startIndex >= length(dataPerSamples)
            startIndex = length(dataPerSamples)-1;
            warnings('The number of windows asked retreives samples outside the dataPerSamples');
        end
        % Could have an incomplete window
        % Make sure the index asked is not outside dataPerSamples        
        stopIndex = startIndex + winLengthInSample;
        if stopIndex > length(dataPerSamples)
            stopIndex = length(dataPerSamples);
        end            
        % Shape the data per windows
        nSamples = stopIndex-startIndex;
        windowsPerSamples(iWind,1:nSamples) = dataPerSamples(startIndex+1:stopIndex);
    end
    
end

