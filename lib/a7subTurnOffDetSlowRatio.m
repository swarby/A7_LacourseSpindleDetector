function [possSpindle, slowRValidTS, slowRInfoTS] = a7subTurnOffDetSlowRatio(...
            possSpindle, PSDLowFreq, PSDHighFreq, artifactDetectVector, DEF_a7)
% Purpose: 
%   Turn off detections when they occur in an unexpected context which can
%   be estimated through a spectral profile.
%
% Input: 
%   possSpindle          : detection vector (double) of the length of the time
%                          series which the a7 detector is run on.
%   PSDLowFreq           : tall vector of the energy in the low band of the slow ratio
%                          [nPSDWindow x 1] (info in PSA window, not in sample)
%   PSDHighFreq          : tall vector of the energy in the high band of the slow ratio
%                          [nPSDWindow x 1] (info in PSA window, not in sample)
%   artifactDetectVector : vector of artifact in sample (tall vector)
%                           (1:artifact;  0:valid)
%   DEF_a7               : structure of a7 settings
%   
% Output:
%   possSpindle - spindle detection vector (same number of datapoints as
%                   input dataVector)
%   slowRValidTS - logical tall vector (0: means not in the NREM spectral
%                   context) same number of datapoints as input dataVector
%   slowRInfoTS - slow ratio information vector (same number of datapoints as
%                   input dataVector)
% 
% Notes : no minimum or maximum length applied on the detections
%
% Authors:
%   Karine Lacourse 2016-08-08
% 
% Arrangement:
%   Jacques Delfrate 2018-02-13
%--------------------------------------------------------------------------
    
    % ---------------------------- INIT -----------------------------------  
    % Total length of the timeseries; in seconds
    dataLength_sec = length(possSpindle)/DEF_a7.standard_sampleRate ;   
    % Number of complete windows based on the length and step in sec.
    % The PSD is performed only on complete window
    nWindows = floor((dataLength_sec- DEF_a7.winLengthSec)/DEF_a7.WinStepSec)+1;   
    % Init the number of samples
    nTotSamples         = length(possSpindle);    
    
    
    % ---------------------------- SCRIPT ---------------------------------
    % Error check
        if abs(nWindows - length(PSDLowFreq))>1
            error('The number of window in PSDLowFreq is unexpected');
        end
        if abs(nWindows - length(PSDHighFreq))>1
            error('The number of window in PSDLowFreq is unexpected');
        end
    
    %----------------------------------------------------------------------
    %% Convert the valid sample vector per PSA window
    %----------------------------------------------------------------------
    invalidWin = samples2WindowsInSec( artifactDetectVector, nWindows, ...
         DEF_a7.winLengthSec, DEF_a7.WinStepSec, DEF_a7.standard_sampleRate);
    invalidWin = sum(invalidWin,2);
    validByWin = (invalidWin==0);    
      
    %----------------------------------------------------------------------
    %% Compute the treshold based on the BSL for each slowRatio
    %----------------------------------------------------------------------
    
    % For each current window compute the index vector of the windows to take to create
    % a clean baseline 
    %   We consider only complete windows, window length should be considered
    %   BLS length is 30 sec - the length of the window (0.3 sec) and the
    %   rest can be splitted into window step window of 0.1 sec
    
    iAvailable      = find(validByWin==1);
    nWinInBSL       = floor( (DEF_a7.bslLengthSec-DEF_a7.winLengthSec)/DEF_a7.WinStepSec ) + 1;
    
    % If there is less valid windows than the number required to compute the
    % baseline (then there less valid windows than 3 mins in the whole
    % recording)
    if isempty(iAvailable)
        error(['Spindle with slow ratio sliding windows: Bsl is empty, ',...
            'DEF_a7.bslLengthSec is set to %i sec, no detection possible'],...
            DEF_a7.bslLengthSec);            
    else
        if nWinInBSL > length(iAvailable)
            if  nargin > 4
                % Output error in the command window with the event name
                warning('%s : Bsl is too short it is %d sec and DEF_a7.bslLengthSec is set to %d sec',...
                    DEF_a7.eventNameSlowRatio, ...
                    round((length(iAvailable)-1)*DEF_a7.WinStepSec+DEF_a7.winLengthSec),...
                    DEF_a7.bslLengthSec);
            else
                % Output error in the command window without the event name
                warning(['Spindle with cov : Bsl is too short, ',...
                    'only %i valid sec and DEF_a7.bslLengthSec is set to %i'], ...
                    round((length(iAvailable)-1)*DEF_a7.WinStepSec+DEF_a7.winLengthSec),...
                    DEF_a7.bslLengthSec);            
            end
            nWinInBSL = length(iAvailable);
        end        
        slowRTotal_med        = nan(nWindows,1);
        
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

            % Use the mean of the energy through all the PSA windows
            if DEF_a7.useMedianWindSlowRatio == 0
                lowTotal_med  = mean( PSDLowFreq(iWinVector), 'omitnan' );
                highTotal_med  = mean( PSDHighFreq(iWinVector), 'omitnan' );
                slowRTotal_med(iWinTot,1) = lowTotal_med ./ highTotal_med;
            % Use the median of the energy through all the PSA windows
            else
                lowTotal_med  = median( PSDLowFreq(iWinVector), 'omitnan');
                highTotal_med  = median( PSDHighFreq(iWinVector), 'omitnan');
                slowRTotal_med(iWinTot,1) = lowTotal_med ./ highTotal_med;
            end
        end

        % All the windows available to compute the BSL
        % If the log10 values are used, add +1 to avoid log10(0)=-inf
        % Adding +1 to the set of values wont influence the distribution
        if DEF_a7.useLog10ValNoNegSlowRatio == 1
            % Create a set of covariance values without negative
            %  we dont care about the negative, it should not happen in spindle
            slowRNoNeg      = slowRTotal_med;
            slowRNoNeg(slowRNoNeg<=0)=0;        
            winBSLNoNaN     = slowRNoNeg(validByWin==1);
            winBSLNoNaN     = log10(winBSLNoNaN+1);
            slowRTotal_med  = log10(slowRNoNeg+1); %update the covariance value     
        end       
        
        %-------------------------------------------------------------------------
        %% Detection 
        %-------------------------------------------------------------------------    

        % We dont care if the event_byWindow is 0 with a covTotal_threshold=nan
        % becasue the time series has already been detected with a previous artifact
        
        % Apply threshold to select NREM samples (estimated)
        slowRValidPerWin     = slowRTotal_med > DEF_a7.slowRat_Th; % (5 > NaN) = 0 

        %----------------------------------------------------------------------
        %% Make a timeseries from the detection information
        %----------------------------------------------------------------------

        % Convert into a sample vector
        % Error check on missing window to have exactly the length of the
        % time series
        if round(size(slowRValidPerWin,1) * DEF_a7.WinStepSec * DEF_a7.standard_sampleRate + ...
                ( DEF_a7.winLengthSec-DEF_a7.WinStepSec) * DEF_a7.standard_sampleRate) < nTotSamples
            slowRValidPerWin = [slowRValidPerWin;0];
            % To output the slowRatio information to debug (learn thresholds)
            slowRTotal_med = [slowRTotal_med;0];
        end
        % We have an overlap between window
        slowRValidTS = windows2SamplesInSec( slowRValidPerWin,  DEF_a7.winLengthSec, ...
            DEF_a7.WinStepSec, DEF_a7.standard_sampleRate, nTotSamples );          
        slowRInfoTS = windows2SamplesInSec( slowRTotal_med,  DEF_a7.winLengthSec, ...
            DEF_a7.WinStepSec, DEF_a7.standard_sampleRate, nTotSamples );  
        
        % Take the max of the overlap values 
        % (if any of the PSA window is an estimation of the NREM)
        slowRValidTS = max(slowRValidTS, [], 1, 'omitnan');
        slowRValidTS = slowRValidTS'; % function is column wise
        slowRInfoTS = mean(slowRInfoTS, 1, 'omitnan');  
        
        % Update the possible spindle with the selection of the NREM2
        % estimated
        possSpindle = possSpindle & slowRValidTS;
    
    end

end

