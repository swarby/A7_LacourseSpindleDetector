function [sigmaCovDistInfoTS] = a7subSigmaCov(timeSeries, ...
    validSampleVect, DEF_a7)
% Purpose
%   Detect spindles based on the covariance between the raw data and 
%   the data filtered in the sigma band.
%   An adaptive threshold is used (distance in STD from the mean)
%   For each window we compute the start index rounded based on window
%   length and window step in sec. The same length in samples is used 
%   for each window (stop = start + window lentgh in samples).  
%
% Input 
%   timeSeries      : input timeseries, detector will run on this data.
%   validSampleVect : selection vector of samples to compute the BSL
%   DEF_a7          : structure of a7 settings
%   
% Output
%   sigmaCovDistInfoTS : vector of the information used to detect, same number of datapoints as input dataVector.
% 
% Requirement
%   butterFiltZPHighPassFiltFilt.m
%   butterFiltZPLowPassFiltFilt.m
%   samples2WindowsInSec.m
%
% Author :     Karine Lacourse  2016-08-08
% Arrangment : Jacques Delfrate 2018-02-13
%--------------------------------------------------------------------------


    % ---------------------------- SCRIPT ---------------------------------
    % Total length of the timeseries; in seconds
    dataLength_sec = length(timeSeries)/DEF_a7.standard_sampleRate ;   
    % Number of windows (complete windows)
    nWindows = floor( (dataLength_sec-DEF_a7.winLengthSec)/DEF_a7.WinStepSec) + 1;
    
    %----------------------------------------------------------------------
    %% Filter out Delta from the signal
    %----------------------------------------------------------------------
    if DEF_a7.removeDeltaFromRawSigmaCov==1
        timeSeriesNDelta = butterFiltZPHighPassFiltFilt( ...
            timeSeries, DEF_a7.totalFreqLow, DEF_a7.standard_sampleRate, DEF_a7.filterOrderSigmaCov);
    else
        timeSeriesNDelta = timeSeries;
    end
    
    %----------------------------------------------------------------------
    %% Filter the signal in the sigma band
    %----------------------------------------------------------------------
    timeSeriesFilt  = butterFiltZPHighPassFiltFilt(...
        timeSeries, DEF_a7.sigmaFreqLow, DEF_a7.standard_sampleRate, ...
        DEF_a7.filterOrderSigmaCov);
    timeSeriesFilt  = butterFiltZPLowPassFiltFilt(...
        timeSeriesFilt, DEF_a7.sigmaFreqHigh, DEF_a7.standard_sampleRate, ...
        DEF_a7.filterOrderSigmaCov );
    
    
    %----------------------------------------------------------------------
    %% Compute the covariance
    %----------------------------------------------------------------------
    % Shape the time series (vector of samples) into windows of samples
    tsRawPerWindow = samples2WindowsInSec( timeSeriesNDelta, nWindows, ...
        DEF_a7.winLengthSec, DEF_a7.WinStepSec, DEF_a7.standard_sampleRate);
    tsSigmaPerWindow = samples2WindowsInSec( timeSeriesFilt, nWindows, ...
        DEF_a7.winLengthSec, DEF_a7.WinStepSec, DEF_a7.standard_sampleRate);
    
    covarianceWin = zeros(nWindows,1);
    % For each window
    for iWind = 1 : nWindows
        % Compute the covariance on the selected window
        covResult = cov(tsRawPerWindow(iWind,:), tsSigmaPerWindow(iWind,:), 'partialrows');
        % Store the covariance value
        covarianceWin(iWind,1) = covResult(1,2);        
    end
    
    %----------------------------------------------------------------------
    %% Convert the valid sample vector per window
    %----------------------------------------------------------------------
    % Shape the time series (vector of samples) into windows of samples
    invalidPerWin   = samples2WindowsInSec( ~validSampleVect, nWindows, ...
        DEF_a7.winLengthSec, DEF_a7.WinStepSec, DEF_a7.standard_sampleRate); 
    covInvalidWin   = sum(invalidPerWin,2, 'omitnan');
    validByWin      = (covInvalidWin==0);
    
    % All the covariance windows available to compute the BSL
    % If the log10 values are used, add +1 to avoid log10(0)=-inf
    % Adding +1 to the set of values wont influence the distribution
    if DEF_a7.useLog10ValNoNegSigmaCov==1
        % Create a set of covariance values without negative
        % we dont care about the negative, it should not happen in spindle
        covNoNeg        = covarianceWin;
        covNoNeg(covNoNeg<=0)=0;        
        winBSLNoNaN     = covNoNeg(validByWin==1);
        winBSLNoNaN     = log10(winBSLNoNaN+1);
        covariance      = log10(covNoNeg+1); %update the covariance value 
    else
    	% Sum the windows without any previous detections
        winBSL              = covarianceWin; 
        winBSL(validByWin==0) = nan;
        winBSLNoNaN         = winBSL(validByWin==1);
        covariance          = covarianceWin;
    end
       
    
    %----------------------------------------------------------------------
    %% Compute the treshold based on the BSL for each sigmaCov
    %----------------------------------------------------------------------
    % For each current cov window compute the index vector of the cov windows to take to create
    % a clean baseline (ex. BLS length is 30 sec and cov window step window is 0.1 sec, then
    % we need 300 clean cov windows to create a clean baseline)
    iAvailable      = find(validByWin==1);
    nWinInBSL       = floor( (DEF_a7.bslLengthSec-DEF_a7.winLengthSec) / DEF_a7.WinStepSec) +1;

    % If there is no baseline
    if isempty(iAvailable)
        sigmaCovDistInfoTS  = nan(nWindows,1);
        error('Spindle with cov : Bsl is empty, DEF_a7.bslLengthSec is set to %i sec, no detection possible',...
            DEF_a7.bslLengthSec);            
    else
       % If there is less valid windows than the number required to compute the
        % baseline 
        if nWinInBSL > length(iAvailable)
            % Output error in the command window
            warning(['Spindle with cov sliding windows: Bsl is too short ',...
                'only %i sec valid and DEF_a7.bslLengthSec is set to %i sec'], ...
                round((length(iAvailable)-1)*DEF_a7.WinStepSec+DEF_a7.winLengthSec),...
                DEF_a7.bslLengthSec);  
            nWinInBSL = length(iAvailable);
            nWindows = length(iAvailable);
        end
        % Inits
        total_threshold     = nan(nWindows,1);
        covTotal_SD         = nan(nWindows,1);
        covTotal_med        = nan(nWindows,1);
        limits              = percentiles(winBSLNoNaN,[DEF_a7.lowPerctSigmaCov,DEF_a7.highPerctSigmaCov]);
        
        % Because we need to detect spindles only on the valid sample
        % any sample not valid will have an detection information set to nan

        % We want to init the baseline threshold for all the windows, 
        % even the ones that includes a previously detected artifact
        % To keep matrixes the same size (NaN are marked anyway)
        for iWinTot = 1 : nWindows

            % Index vector of valid PSA window
            iValWin    = find(iAvailable>=iWinTot,1,'first');
            % If there is no more PSD windows without any artifact (end of the
            % recording is all artifacted) we take the last valid PSD window.
            if isempty(iValWin)
                iValWin = iAvailable(end);
            end

            % If the number of PSD windows to take is even (ex. 52) 
            % one more previous PSD window is taken than the following windows
            %   (ex. 1:26 + 27:27+25)
            % If the number of PSD windows to take is odd (ex. 45) 
            % the same number of previous and following PSD windows are taken
            %   (ex. 1:22 + 23:23+22)

            % At least 1 previous PSD is missing
            if iValWin <= ceil((nWinInBSL-1)/2)
                iWinVector   = iAvailable(1:nWinInBSL);
            % At least 1 following PSD is missing
            elseif iValWin >= (length(iAvailable)-floor((nWinInBSL-1)/2))
                iWinVector   = iAvailable(end-nWinInBSL+1:end);       
            % Enough PSD windows around the current PSD window
            else
                iStart       = iValWin-ceil((nWinInBSL-1)/2);
                iStop        = iStart + nWinInBSL;
                iWinVector   = iAvailable(iStart+1:iStop);       
            end

            % Consider only the baseline included in the percentile selected
            if DEF_a7.useLimPercSigmaCov == 1
                BSLsort     = sort(covariance(iWinVector));
                iStart      = find(BSLsort>=limits(1));
                iStop       = find(BSLsort<=limits(2));
                iInPerc     = intersect(iStart,iStop);
                covTotal_SD(iWinTot,1)       = std( BSLsort(iInPerc), 'omitnan' );
                if DEF_a7.useMedianWindSigmaCov==0
                    covTotal_med(iWinTot,1)  = mean( BSLsort(iInPerc), 'omitnan' );
                else
                    covTotal_med(iWinTot,1)  = median( BSLsort(iInPerc), 'omitnan');
                end            
            else
                covTotal_SD(iWinTot,1)       = std( covariance(iWinVector), 'omitnan' );
                if DEF_a7.useMedianWindSigmaCov==0
                    covTotal_med(iWinTot,1)  = mean( covariance(iWinVector), 'omitnan' );
                else
                    covTotal_med(iWinTot,1)  = median( covariance(iWinVector), 'omitnan');
                end
            end
            % the STD is forced to 0 when the bsl available to compute 
            % the STD is outside the percentile limits
            if isnan(covTotal_SD(iWinTot))
                covTotal_SD(iWinTot)     = 0;   
            end
            if isnan(covTotal_med(iWinTot))
                covTotal_med(iWinTot)     = 0;   
            end
            % Only for validtion purpose
            total_threshold(iWinTot) = covTotal_med(iWinTot) + ...
                (DEF_a7.sigCov_Th * covTotal_SD(iWinTot)); 
        end

            
        % Compute the threshold based on threshold * STD + median
        % Only for validation purpose
        covTotal_threshold = covTotal_med + (DEF_a7.sigCov_Th * covTotal_SD);
    
        %-------------------------------------------------------------------------
        %% Detection 
        %-------------------------------------------------------------------------    

        % We dont care if the event_byWindow is 0 with a covTotal_threshold=nan
        % becasue the time series has already been detected with a previous artifact
        % !not necessary it could be because its outside the percentile limits!!

        % Only for validation purpose
        eventPerWin     = covariance > covTotal_threshold; % (5 > NaN) = 0 

        %--------------------------------------------------------------------------
        % Calculate the standard deviation distance from the threshold based on the baseline
        % To keep track of the detection window information
        %--------------------------------------------------------------------------
        covDistWin      = (covariance-covTotal_med)./covTotal_SD;

        % Verify that artifacts computed with the std are the
        % same than computed with the fixed threshold from the baseline
        % Only for validation purpose
        sameDetection = (covDistWin > DEF_a7.sigCov_Th) == eventPerWin;
        if ~isempty(find(sameDetection==0,1,'first'))
            error('Computation of the threshold on the baseline is wrong');
        end    

        %----------------------------------------------------------------------
        %% Make a timeseries from the detection information
        %----------------------------------------------------------------------

        % Convert the cov distance into a sample vector
        % Error check on missing window to have exactly the length of the
        % time series
        % The number of samples created by windows2Samples with nWindows
        nSamplesFromWin = round(nWindows * DEF_a7.WinStepSec*DEF_a7.standard_sampleRate + ...
            (DEF_a7.winLengthSec*DEF_a7.standard_sampleRate-DEF_a7.WinStepSec*DEF_a7.standard_sampleRate));
        if nSamplesFromWin < length(timeSeries)
            covDistWin = [covDistWin;0];
        end       
        distWinMat = windows2SamplesInSec( covDistWin, DEF_a7.winLengthSec, ...
            DEF_a7.WinStepSec, DEF_a7.standard_sampleRate, length(timeSeries) );            

        % Average all the cov from the overlapped windows
        distWinMean     = mean(distWinMat,'omitnan');
        distWinMean     = fillmissing(distWinMean,'previous');

        % Output the detection information for analysis of the threshold
        sigmaCovDistInfoTS   = distWinMean;
    
    end

end

