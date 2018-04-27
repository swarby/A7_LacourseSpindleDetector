function [sigmaCorrInfoTS] = a7subSigmaCorr(timeSeries, DEF_a7)
% Purpose 
%   Detect spindles based on the correlation between the raw data and 
%   the data filtered in the sigma band.
%   A fixed threshold is used (% of correlation r) r can vary from -1 to 1
%   but we don't want any delay between the raw signal and the sigma signal
%   then we threshold on the positive part ex) 0.7 or 70%.
%
% Input
%   dataVector      - input timeseries, detector will run on this data.
%   DEF_a7.standard_sampleRate      - DEF_a7.standard_sampleRate of the timeseries data, needed for the calcPSD.
%   validSampleVect - selection vector of samples to show the BSL histogram
%   channelName     - name of the input channel; transferred to the EVENTS output.
%   DEF_a7 - DEF option for the a7 spindle
%       DEF_a7.corrThresh      = 40;   % percentage of correlation threshold 
%                                      % of the sigma correlation over the raw signal 
%       DEF_a7.relWindLength   = 0.5;    % window length in sec for sigmaCorr
%       DEF_a7.relWindStep     = 0.1;    % window step in sec for sigmaCorr
%       DEF_a7.sigmaFreqHigh   = 16;   % frequency band of the broad band
%       DEF_a7.sigmaFreqLow    = 11;   % frequency band of the sigma
%   
% Output
%   sigmaCorrInfoTS - vector of the information used to detect, same number of datapoints as input dataVector.

%
% Requirements
%   event_StartsEndsDurations.m
% 
% Notes : no minimum or maximum length applied on the detections
%
% Authors
% Karine Lacourse  2016-08-08
% 
% Arrangement
% Jacques Delfrate 2018-02-13
%--------------------------------------------------------------------------
    
    % We dont manage the windows per samples because of rounding error.
    % We want to minimize the difference between different sampling rate.
    % For each window we compute the start index rounded based on window
    % length and window step in sec. The same length in samples is used 
    % for each window (stop = start + window lentgh in samples).  
    % We want to avoid a shift accumulated through the whole recording due 
    % to rounding error (sampling rate).
    
    % Total length of the timeseries; in seconds
    dataLength_sec = length(timeSeries)/DEF_a7.standard_sampleRate ;   
    % Number of complete windows based on the length and step in sec.
    % nWindows = floor((dataLength_sec-DEF_a7.absWindLength)/DEF_a7.absWindStep)+1;  
    % Number of windows (the maximum number of step windows, at least half)
    nWindows = round(dataLength_sec/DEF_a7.absWindStep);    
    % Init the number of samples
    nTotSamples             = length(timeSeries);
    
%     DEF_a7.filterOrderSigCorr             = 20; % The filter order

    % Filter out Delta from the signal
    if DEF_a7.removeDeltaFromRawSigCorr  ==1
        timeSeriesNDelta = butterFiltZPHighPassFiltFilt( ...
            timeSeries, DEF_a7.totalFreqLow, DEF_a7.standard_sampleRate, ... 
            DEF_a7.filterOrderSigCorr);
    else
        timeSeriesNDelta = timeSeries;
    end
    
    % Filter the signal in the sigma band
    timeSeriesFilt  = butterFiltZPHighPassFiltFilt(...
        timeSeries, DEF_a7.sigmaFreqLow, DEF_a7.standard_sampleRate, ...
        DEF_a7.filterOrderSigCorr);
    timeSeriesFilt  = butterFiltZPLowPassFiltFilt(...
        timeSeriesFilt, DEF_a7.sigmaFreqHigh, DEF_a7.standard_sampleRate, ...
        DEF_a7.filterOrderSigCorr);
    
    % Shape the time series (vector of samples) into windows of samples
    tsRawPerWindow = samples2WindowsInSec( timeSeriesNDelta, nWindows, ...
        DEF_a7.relWindLength, DEF_a7.relWindStep, DEF_a7.standard_sampleRate);
    tsSigmaPerWindow = samples2WindowsInSec( timeSeriesFilt, nWindows, ...
        DEF_a7.relWindLength, DEF_a7.relWindStep, DEF_a7.standard_sampleRate);
    
    % Compute the correlation
    correlationWin = zeros(nWindows,1);
    % For each window
    for iWind = 1 : nWindows
        % Compute the correlation on the selected window
        rCoeff = corrcoef(tsRawPerWindow(iWind,:), tsSigmaPerWindow(iWind,:));
        % Store the correlation value
        correlationWin(iWind,1) = rCoeff(1,2);        
    end

    % Error check on missing window to have exactly the length of the
    % time series
    % The number of samples created by windows2Samples with nWindows
    nSamplesFromWin = round(nWindows * DEF_a7.relWindStep*DEF_a7.standard_sampleRate + ...
        (DEF_a7.relWindLength*DEF_a7.standard_sampleRate-DEF_a7.relWindStep*DEF_a7.standard_sampleRate));
    if nSamplesFromWin < length(timeSeries)
        correlationWin = [correlationWin;0];
    end      
    
    % Convert the correlation into a sample vector
    correlationMat = windows2SamplesInSec( correlationWin, DEF_a7.relWindLength, ...
        DEF_a7.relWindStep, DEF_a7.standard_sampleRate, nTotSamples );
    
    % Average all the corr from the overlapped windows
    correlation = mean(correlationMat,'omitnan');
    
    % Output the detection information for analysis of the threshold
    sigmaCorrInfoTS = correlation';    
   
end

