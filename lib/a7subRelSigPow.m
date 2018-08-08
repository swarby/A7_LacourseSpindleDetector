
function [relSigmaDistVect, PSDLowFreq, PSDHighFreq] = ...
    a7subRelSigPow(dataVector, validSampleVect, DEF_a7)

% Purpose
%   Detect spindles based on the relative sigma power. 
%   An adaptive threshold is used. The baseline is computed on only clean data.
%   The threshold is specified in STD distance from the mean. 
%
% Input:
%   dataVector      - input timeseries, detector will run on this data.
%   channelName     - name of the input channel; transferred to the EVENTS output.
%   validSampleVect - selection vector of samples to compute the BSL
%   DEF_a7 - DEF option for the a7 spindle
%       DEF_a7.PSAWindLength       = 0.5;    % window length in sec for PSA
%       DEF_a7.PSAZWindLength      = 2;      % window length with zero pad in sec for PSA
%       DEF_a7.PSAWindStep         = 0.1;    % window step in sec for PSA
%       DEF_a7.bslLengthSec        = 30;     % length of the baseline in sec
%       DEF_a7.standard_sampleRate = 100     % sampleRate of the timeseries data, needed for the calcPSD.
%       DEF_a7.relSigPow_Th                  % relative sigma power threshold
%   
% Output
%   relSigmaDistVect - vector of the information used to detect, same number of datapoints as input dataVector.
%   PSDLowFreq - tall vector of the energy in the low band of the slow ratio
%               [nPSDWindow x 1] (info in PSA window, not in sample)
%   PSDHighFreq - tall vector of the energy in the high band of the slow ratio
%               [nPSDWindow x 1] (info in PSA window, not in sample)
%
% Requirements
%   calcPSD.m 
%   event_StartsEndsDurations.m 
%   findFreqIndex.m
%
% Authors:
% Karine Lacourse  2016-08-08
% 
% Arrangement :
% Jacques Delfrate 2018-02-13
%--------------------------------------------------------------------------

% Verification error
if(iscell(dataVector))
    error('%s: data must be a vector', DEF_a7.eventNameRelSigPow);
end

% Check if dataVector is a column vector
if ~iscolumn(dataVector)
    dataVector = dataVector';
end

%-------------------------------------------------------------------------
%% Compute the relative sigma power on every PSA window
%-------------------------------------------------------------------------

    % Compute the PSD - each row is one FFT window, each column is freqBin
    % hard coded
    windowScalingType       = 'hann';
    removeMean              = 1;
    zeroPad                 = 1; % To increase the freq resolution
    normType                = 'IntegSpectPow';
    % Only complete windows are generated for the PSD
    [PSD, freqBins] = calcPSDNormType(dataVector, DEF_a7, ...
        windowScalingType, removeMean, zeroPad, normType);
    nPSDWindow = size(PSD,1);

    % Find the index of the frequency band in the frequency bins of the PSD
    [ iFreqStart, iFreqStop ] = findFreqIndex( DEF_a7.sigmaFreqLow, ...
        DEF_a7.sigmaFreqHigh,  freqBins );
    if iFreqStart == iFreqStop
        error(['The frequency band to compute the relative sigma band is excluded',...
            'of the frequency bins of the PSD']);
    end
    % Calculate the power by summing each row (each FFT window); matrix transposed, as sum works by-column
    PSDSigmaFreq         = sum(PSD(:,iFreqStart:iFreqStop),2); 

    % Find the index of the frequency band in the frequency bins of the PSD
    [ iFreqStart, iFreqStop ] = findFreqIndex( DEF_a7.totalFreqLow, ...
        DEF_a7.totalFreqHigh,  freqBins );
    if iFreqStart == iFreqStop
        error(['The frequency band to compute the relative sigma band is excluded',...
            'of the frequency bins of the PSD']);
    end
    PSDTotalFreq         = sum(PSD(:,iFreqStart:iFreqStop),2); 

    relSigmaPow          = PSDSigmaFreq ./ PSDTotalFreq;
    % error check
    if any(relSigmaPow > 1)
       error('The relative sigma power cannot be higher than 100%, check the code!'); 
    end
   
%----------------------------------------------------------------------
% Extract energy in low band and high band in order to turn off
% detection outside n2-3-4
%----------------------------------------------------------------------
if DEF_a7.inContOn == 1
    % Find the index of the frequency band in the frequency bins of the PSD
    % band of delta + theta
    [ iFreqStart, iFreqStop ] = findFreqIndex( DEF_a7.lowFreqLow, ...
        DEF_a7.lowFreqHigh,  freqBins );
    if iFreqStart == iFreqStop
        error(['The frequency band to compute the low band (slowRatio) is excluded',...
            'of the frequency bins of the PSD']);
    end
    PSDLowFreq         = sum(PSD(:,iFreqStart:iFreqStop),2);     
    % band of beta
    [ iFreqStart, iFreqStop ] = findFreqIndex( DEF_a7.highFreqLow, ...
        DEF_a7.highFreqHigh,  freqBins );
    if iFreqStart == iFreqStop
        error(['The frequency band to compute the high band (slowRatio) is excluded',...
            'of the frequency bins of the PSD']);
    end
    PSDHighFreq         = sum(PSD(:,iFreqStart:iFreqStop),2);
else
    PSDLowFreq  = [];
    PSDHighFreq = [];
end

%-------------------------------------------------------------------------
%% Compute the treshold based on the BSL for each relSigmaPow
%-------------------------------------------------------------------------

    % Convert the valid sample vector per PSA window
    invalidSampleByPSDWin   = samples2WindowsInSec( ~validSampleVect, ...
        size(PSD,1), DEF_a7.PSAWindLength, DEF_a7.PSAWindStep, DEF_a7.standard_sampleRate);
    invalidSampleByPSDWin   = sum(invalidSampleByPSDWin,2);
    validSampleByPSDWin     = (invalidSampleByPSDWin==0);

    % Sum the PSD of windows without any previous detections
    PSDBSL              = relSigmaPow; 
    PSDBSL(validSampleByPSDWin==0) = nan;
    PSDBSLNoNaN         = PSDBSL(validSampleByPSDWin==1);

    % The relative sigma power has to be log10 transform
    PSDBSL      = log10(PSDBSL);
    PSDBSLNoNaN = log10(PSDBSLNoNaN);
    relSigmaPow = log10(relSigmaPow);

    % For each PSD window compute the index vector of the PSD windows to take to create
    % a clean baseline (ex. BLS length is 3 mins and PSD window is 2 sec, then
    % we need 45 clean PSD windows to create a clean baseline)
    iAvailable          = find(validSampleByPSDWin==1);
    nFFTWinInBSL        = round(DEF_a7.bslLengthSec/DEF_a7.PSAWindStep);

    % If there is less valid windows than the number required to compute the
    % baseline (then there less valid windows than 3 mins in the whole
    % recording)
    if nFFTWinInBSL > length(iAvailable)
        % If there is the warningsDir input argument
        if nargin == 4
            % Output error in the warning directory
            warning('Bsl empty because there is only %i valid PSD and we need %i',...
                length(iAvailable), nFFTWinInBSL);
             warning('%s : Bsl empty because there is only %i valid window and we need %i',...
                   eventName, length(iAvailable), nWinInBSL);
        else
            if  nargin > 1
                % Output error in the command window
                warning('%s : Bsl empty because there is only %i valid PSD and we need %i',...
                    DEF_a7.eventNameRelSigPow, length(iAvailable), nFFTWinInBSL);
            else
                % Output error in the command window
                warning(['Spindle with Brunner : Bsl empty because there',...
                    'is only %i valid PSD and we need %i'], length(iAvailable),...
                    nFFTWinInBSL);            
            end
        end
        relSigmaDistVect    = [];
    else
        PSDtotal_threshold  = nan(nPSDWindow,1);
        PSDtotal_SD         = nan(nPSDWindow,1);
        PSDtotal_med        = nan(nPSDWindow,1);
        limits              = percentiles(PSDBSLNoNaN, ... 
            [DEF_a7.lowPerctRelSigPow,DEF_a7.highPerctRelSigPow]);

        % We want to init the baseline threshold for all the PSD windows, 
        % even the ones that includes a previously detected artifact
        % To keep matrixes the same size (NaN are marked anyway)
        for iPSDWinTot = 1 : nPSDWindow

            % Index vector of valid PSA window
            iValPSDWin    = find(iAvailable>=iPSDWinTot,1,'first');
            % If there is no more PSD windows without any artifact (end of the
            % recording is all artifacted) we take the last valid PSD window.
            if isempty(iValPSDWin)
                iValPSDWin = iAvailable(end);
            end

            % If the number of PSD windows to take is even (ex. 52) 
            % one more previous PSD window is taken than the following windows
            %   (ex. 1:26 + 27:27+25)
            % If the number of PSD windows to take is odd (ex. 45) 
            % the same number of previous and following PSD windows are taken
            %   (ex. 1:22 + 23:23+22)

            % At least 1 previous PSD is missing
            if iValPSDWin <= ceil((nFFTWinInBSL-1)/2)
                iPSDWinVector   = iAvailable(1:nFFTWinInBSL);

            % At least 1 following PSD is missing
            elseif iValPSDWin >= (length(iAvailable)-floor((nFFTWinInBSL-1)/2))
                iPSDWinVector   = iAvailable(end-nFFTWinInBSL+1:end);       
            % Enough PSD windows around the current PSD window
            else
                iStart          = iValPSDWin-ceil((nFFTWinInBSL-1)/2);
                iStop           = iStart + nFFTWinInBSL;
                iPSDWinVector   = iAvailable(iStart+1:iStop);       
            end

            % Consider only the baseline included in the percentile selected
            if DEF_a7.useLimPercRelSigPow == 1
                PSDBSLsort  = sort(PSDBSL(iPSDWinVector));
                iStart      = find(PSDBSLsort>=limits(1));
                iStop       = find(PSDBSLsort<=limits(2));
                iInPerc     = intersect(iStart,iStop);
                PSDtotal_SD(iPSDWinTot,1)       = std( PSDBSLsort(iInPerc), 'omitnan' );
                if DEF_a7.useMedianPSAWindRelSigPow==0
                    PSDtotal_med(iPSDWinTot,1)  = mean( PSDBSLsort(iInPerc), 'omitnan' );
                else
                    PSDtotal_med(iPSDWinTot,1)  = median( PSDBSLsort(iInPerc), 'omitnan');
                end            
            else
                PSDtotal_SD(iPSDWinTot,1)       = std( PSDBSL(iPSDWinVector), 'omitnan' );
                if DEF_a7.useMedianPSAWindRelSigPow==0
                    PSDtotal_med(iPSDWinTot,1)  = mean( PSDBSL(iPSDWinVector), 'omitnan' );
                else
                    PSDtotal_med(iPSDWinTot,1)  = median( PSDBSL(iPSDWinVector), 'omitnan');
                end
            end
            % the STD is forced to 0 when the bsl available to compute 
            % the STD is outside the percentile limits
            if isnan(PSDtotal_SD(iPSDWinTot))
                PSDtotal_SD(iPSDWinTot)     = 0;   
            end
            if isnan(PSDtotal_med(iPSDWinTot))
                PSDtotal_med(iPSDWinTot)     = 0;   
            end
            % For validation purpose
            PSDtotal_threshold(iPSDWinTot) = PSDtotal_med(iPSDWinTot) + ...
                (DEF_a7.relSigPow_Th * PSDtotal_SD(iPSDWinTot)); 
        end      
        
        %-------------------------------------------------------------------------
        %% Detection based on the relSigmaPow and the baseline
        %-------------------------------------------------------------------------    

        % We dont care if the event_byWindow is 0 with a PSDtotal_threshold=nan
        % becasue the time series has already been detected with a previous artifact
        % not necessary it could be because its outside the percentile limits!!
        
        % Only for validation purpose
        event_byWindow          = relSigmaPow > PSDtotal_threshold; % (5 > NaN) = 0     

        %--------------------------------------------------------------------------
        % Calculate the standard deviation distance from the threshold based on the baseline
        % To keep track of the detection window information
        %--------------------------------------------------------------------------
        PSDDistWin              = (relSigmaPow-PSDtotal_med)./PSDtotal_SD;

        % Verify that artifacts computed with the std are the
        % same than computed with the fixed threshold from the baseline
        % Only for validation purpose
        sameDetection = (PSDDistWin > DEF_a7.relSigPow_Th) == event_byWindow;
        if ~isempty(find(sameDetection==0,1,'first'))
            error('Computation of the threshold on the baseline is wrong');
        end

        %----------------------------------------------------------------------
        %% Make a timeseries from the detection information
        %----------------------------------------------------------------------

        % Convert the relSigmaPow distance into a sample vector
        % Error check on missing window to have exactly the length of the
        % time series
        if round(size(PSDDistWin,1) * DEF_a7.PSAWindStep * DEF_a7.standard_sampleRate + ...
                ((DEF_a7.PSAWindLength * DEF_a7.standard_sampleRate)-(DEF_a7.PSAWindStep * DEF_a7.standard_sampleRate))) ...
                < length(dataVector)
            PSDDistWin = [PSDDistWin;0];
        end
        PSDDistWinMat = windows2SamplesInSec( PSDDistWin, DEF_a7.PSAWindLength, ...
            DEF_a7.PSAWindStep, DEF_a7.standard_sampleRate, length(dataVector) );            

        % The resolution is DEF_a7.PSAWindStep, but we consider only the maximum
        % relative sigma power through all the overlapped PSA window
        PSDDistWinMean     = mean(PSDDistWinMat,'omitnan');
        PSDDistWinSTD      = std(PSDDistWinMat,'omitnan');

        % This information should be visually inspected on NYX viewer
        % or summarize for every spindle from the gold standard
        relSigmaDistVect   = PSDDistWinMean;

    end

end
