function [detVect, detInfoTS, NREMClass, outputFile] = ...
    a7SpindleDetection(timeSeries, sleepStageVect, ...
    artifactDetectVector, DEF_a7)
% Purpose
%   Compute features related to sigma power to detect spindles.
%
% Input
%   timeSeries           : time series (usually C3 channel) tall vector
%   sleepStageVect       : sleep stage vector tall vector
%   artifactDetectVector : artifact detection vector (logical 0 or 1) tall vector
%   DEF_a7 definition    :
%                          standard_sampleRate : sampling rate 
%                          minimumDuration : minimum accepted length for a spindle (! not used yet!!)
%                          maximumDuration : maximum accepted length for a spindle (! not used yet!!)
%
% Output
%   detVect    : detection vector (logical 0 or 1) size as timeSeries
%   detInfoTS  : matrix of detection info converted into a time series 
%                (with the same size than timeSeries)
%                first 4 rows are the features used to detect spindles
%                1 : PSDSigmaLog (sigma power log10 transformed)
%                2 : relSigPow (z-score of relative sigma power from 30 sec clean around the PSA window)
%                3 : sigmaCov (z-score of the covariance of sigma from 30 sec clean closed to the PSA window)
%                4 : sigmaCorr (correlation of sigma)
%                5 : alphaLog (alpha power log10 transformed)
%                6 : PSDLowFreq (delta + theta power for slow ratio)
%                7 : PSDHighFreq (beta power for slow ratio)
%                8 : artifact vector (logical vector 0: no artifact, 1: artifact)
%   NREMClass  : logical tall vector (0/1 In/Out of NREM spectral context)
%   outputFile : output matrix organised to be saved as a text files
%                It contains start (sample), end (sample), duration 
%                (sample) and spectral context of the detected events
% 
% Requirements
%    a7subAbsPowValues.m
%    a7subRelSigPow.m
%    a7subSigmaCov.m
%    a7subSigmaCorr.m
%    a7subDetWithIntv4featAnd.m
%
% 
% Author      : Karine Lacourse 2016-08-12
% Arrangement : Jacques Delfrate 2018-02-13
%-------------------------------------------------------------------------

    nDetecInfo  = 1; 
    detInfoTS   = zeros(nDetecInfo,length(timeSeries));

    %% Compute the valid samples vector (right sleep stage without artifact)
    % If no sleep staging have been chosen
    if isempty(DEF_a7.bslSleepStaging)
        validSampleVect = ones(size(timeSeries));
    else
        % keep TS in order
        validSampleVect = ismember(sleepStageVect, DEF_a7.bslSleepStaging);
    end
    % Take only samples without any artifact
    validSampleVect = (validSampleVect & artifactDetectVector==0);
    nValidSec       = floor(sum(validSampleVect)/DEF_a7.standard_sampleRate)-DEF_a7.relWindLength;
    nValidWindows   = nValidSec/DEF_a7.relWindStep;
    nWinInBSL       = ceil((DEF_a7.BSLLengthSec-DEF_a7.relWindLength) / DEF_a7.relWindStep);
    
    % warning check
    if nValidWindows < nWinInBSL
        warning('Bsl empty because there is only %i valid windows and we need %i',...
            nValidWindows, DEF_a7.BSLLengthSec);
    end

    %% Compute detection info based on the absolute sigma power
    absSigmaPow    = a7subAbsPowValues(timeSeries, DEF_a7);
    % Output to evaluate the threshold
    detInfoTS(1,:) = log10(absSigmaPow);

    %% Compute detection info based on the relative sigma power
        % relative Sigma Power is converted into z-score
        [relSigPow, PSDLowFreq, PSDHighFreq] = ...
            a7subRelSigPow(timeSeries, validSampleVect, DEF_a7);
        if ~isempty(relSigPow)
            % Output to evaluate the threshold
            detInfoTS(2,:) = relSigPow;
            detInfoTS(5,:) = artifactDetectVector; 
        end

    if ~isempty(relSigPow)
    %% Compute the detection info based on the sigma covariance
        sigmaCovDistInfoTS = a7subSigmaCov(timeSeries, ...
            validSampleVect, DEF_a7);
        % Output to evaluate the threshold
        detInfoTS(3,:)     = sigmaCovDistInfoTS;

    %% Compute the detection info based on the sigma correlation
        sigmaCorrInfoTS = a7subSigmaCorr(timeSeries, DEF_a7);
        % Output to evaluate the threshold
        detInfoTS(4,:)  = sigmaCorrInfoTS; 

    %% Combine detection
    % Detect spindles based on the 4 parameters and a7 features definition 
    [detVect, slowRInfoTS, NREMClass] = a7subDetWithIntv4feat(...
            detInfoTS, DEF_a7, artifactDetectVector, PSDLowFreq, PSDHighFreq);
    % slow ratio output    
    detInfoTS(6,:) = slowRInfoTS;
    detInfoTS      = detInfoTS';
    NREMClass  = cell2mat(NREMClass);
    outputFile = organiseOutputTxtFile(detVect, NREMClass, ...
        sleepStageVect, DEF_a7);
    else
        detVect = [];
        detInfoTS = [];
        NREMClass = [];
    end
end


