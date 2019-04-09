% Convert a data per window in a data per sample.
function [ dataMat ] = windows2SamplesInSec( dataPerWindow, winLengthInSec, ...
    winStepInSec, samplingRate, nDataSamples )
%-------------------------------------------------------------------------
% Convert a data per window in a data per sample.
% If sliding windows are used (winStepInSec < winLengthInSec) 
% the output is an array [overlapped windows X nDataSamples] to support 
% overlapped window
% 
% example to use it:
% If thresholds are computed per window and you want to keep track of them
% per sample (you can have more than on thresold per sample (with the overlap).
% With this function you can easily take the maximum (or the mean) threshold
% through the windows overlapped. 
%
% Input
%   dataPerWindow:      column data vector per window
%   winLengthInSec:     (double) Number of sec included in a window
%   winStepInSec:       (double) Number of sec to step for the following window
%                       (a new window is considered at each
%                       "winStepInSec" sec)
%   samplingRate:       (double) Sampling rate in Hz
%   nDataSamples:       number of samples in the data per sample vector
%   
% Output
%   dataMat:            row data vector per sample
%                       matrix if there is supperposition of window
%                       ex: nSamplesPerWindow = 3; nSamplesPerInterv = 2;                     
%                           dataMat =   [1,  1,  1,NaN,3,3,3...
%                                        NaN,NaN,2,2,  2,NaN,4...]
%                       ex: nSamplesPerWindow = 4; nSamplesPerInterv = 2;                     
%                           dataMat =   [1,1,1,1,3,3,3,3...
%                                        NaN,NaN,2,2,2,2,4...]
%                       ex: nSamplesPerWindow = 4; nSamplesPerInterv = 1;  
%                           dataPerWindow = [1;2;3;4;5;6;7;8;9;10]
%                           dataMat =   [1,  1,  1,  1,5,5,5,5,9, 9, 9,   9,  NaN
%                                        NaN,2,  2,  2,2,6,6,6,6,10,10,   10, 10...
%                                        NaN,NaN,3,  3,3,3,7,7,7, 7, NaN, NaN, NaN...
%                                        NaN,NaN,NaN,4,4,4,4,8,8, 8, 8,   NaN, NaN]                       
% Author : Karine Lacourse
% Date   : 2016-10-24
% -------------------------------------------------------------------------

    % number of windows
    nWindows            = size(dataPerWindow,1);
    % These variables need to be integers
    nSamplesPerWindow   = round(winLengthInSec * samplingRate);  
    nDataSamples        = round(nDataSamples);
    
    % Error check
    if (nDataSamples) > round((nWindows * winStepInSec * samplingRate + ...
            (nSamplesPerWindow-winStepInSec*samplingRate)))
       error(['data per window are missing to generate the data per sample',...
           '\nNb. of samples asked is %d and the number samples possibly generated is %d'],...
           nDataSamples, round((nWindows * winStepInSec*samplingRate + ...
            (nSamplesPerWindow-winStepInSec*samplingRate))));
    end   
    
    % compute how many windows are supperposed
    supperposition  = ceil(winLengthInSec / winStepInSec);
    dataMat         = nan(supperposition, nDataSamples);
    iWin = 1;
    while iWin <= nWindows
        for iOverlapWin = 1 : supperposition
            iStart = round((iWin-1) * winStepInSec * samplingRate);
            iStop  = iStart + nSamplesPerWindow;          
            if iWin <= nWindows
                % The window is complete
                if iStop <= nDataSamples
                    dataMat(iOverlapWin, iStart + 1: iStop) = ...
                        repmat(dataPerWindow(iWin),1,nSamplesPerWindow);
                % The window is incomplete
                else
                    nSamplesIn = nDataSamples-iStart;
                    dataMat(iOverlapWin, iStart + 1: end) = ...
                        repmat(dataPerWindow(iWin),1,nSamplesIn);                
                end
            end
            % watch out iWin is incremented with each overlap
            iWin = iWin + 1;
        end
    end
end

