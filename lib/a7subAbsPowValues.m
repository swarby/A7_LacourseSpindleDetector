function absSigmaPow = a7subAbsPowValues( dataVector, DEF_a7)
% Purpose
%   Compute the absolute power for different frequency band through the
%   average sample². Absolute power : absSigmaPow
%
% Input
%   dataVector : time series to compute the abs sigma power (vector)
%   DEF_a7.standard_sampleRate : sampling rate 
%   validSampleVect - selection vector of samples to compute the BSL (only
%       to see histogram, there is no normalization to the baseline for the
%       absSigmaPow.
%   DEF_a7     : struct of options
%               DEF_a7.standard_sampleRate :  sampling rate
%               DEF_a7.absWindLength : window length in seconds
%               DEF_a7.absWindStep : window step in seconds
%               DEF_a7.sigmaFreqLow : the low frequency of the sigma band
%               DEF_a7.sigmaFreqHigh : the high frequency of the sigma band
%
% Output (wide vector, size of dataVector transposed)
%   absSigmaPow
%
% Requirements
%   butterFiltZPHighPassFiltFilt.m
%   butterFiltZPLowPassFiltFilt.m
%   samples2WindowsInSec.m
% 
% Author      : Karine Lacourse  2016-10-24
% Arrangement : Jacques Delfrate 2018-02-13
%-------------------------------------------------------------------------

    % Verification error
    if(iscell(dataVector))
        error('%s: data must be a vector', DEF_a7.eventNameAbsPowValue );
    end
    % Check if dataVector is a column vector
    if ~iscolumn(dataVector)
        dataVector = dataVector';
    end
    lengthWindowInSample    = round(DEF_a7.absWindLength * DEF_a7.standard_sampleRate);
    
    % Total length of the timeseries; in seconds
    dataLength_sec = length(dataVector)/DEF_a7.standard_sampleRate ;   
    % Number of windows (the maximum number of step windows, at least half)
    nWindows = round(dataLength_sec/DEF_a7.absWindStep);     

    %--------------------------------------------------------------------
    %% Compute the RMS (wihtout the root)
    %--------------------------------------------------------------------

    % ------ SIGMA ---------
    % Filter the signal
    timeSeriesFilt  = butterFiltZPHighPassFiltFilt(...
        dataVector, DEF_a7.sigmaFreqLow, DEF_a7.standard_sampleRate, 20);
    timeSeriesFilt  = butterFiltZPLowPassFiltFilt(...
        timeSeriesFilt, DEF_a7.sigmaFreqHigh, DEF_a7.standard_sampleRate, 20);
    % Convert the vector per sample into a matrix 
    % [nWindow x windowLengthInSample]
    sampleMat4E   = samples2WindowsInSec(timeSeriesFilt, nWindows, ...
        DEF_a7.absWindLength, DEF_a7.absWindStep, DEF_a7.standard_sampleRate);
    % The last samples can be NaN if the window is incomplete
    sigma_EValPerWin = sum(sampleMat4E.^2/lengthWindowInSample,2, 'omitnan' );

    %--------------------------------------------------------------------
    %% Convert the values per window into a vector per sample
    %--------------------------------------------------------------------
    % In case of overlap between windows : average the value between window

    % Convert the absolute sigma power into a sample vector
    ESigmaFreqMat = windows2SamplesInSec( sigma_EValPerWin, DEF_a7.absWindLength, ...
        DEF_a7.absWindStep, DEF_a7.standard_sampleRate, length(dataVector) );    

    % The resolution is stepWindowInSample, but we consider only the mean
    % power through all the overlapped window
    absSigmaPow     = mean(ESigmaFreqMat,'omitnan');
    
        % *** Sigma ***
        nMissSamples = length(dataVector) - length(absSigmaPow);
        absSigmaPow = [absSigmaPow, repmat(absSigmaPow(end), 1, nMissSamples)];
        % error check
        % Make sure the converted energy per samples has the same length than dataVector
        if nMissSamples > round((lengthWindowInSample - DEF_a7.absWindStep * DEF_a7.standard_sampleRate)/2)
            warning('Window managment is weird, too many missing samples');
        end

end

