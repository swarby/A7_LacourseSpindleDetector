function [detectionVectFix, slowRInfoTS, NREMClass] = ...
    a7subDetWithIntv4feat(detInfoTS, DEF_a7, artifactDetectVector, ...
    PSDLowFreq, PSDHighFreq)
%
% Purpose : Detect spindles based on the 4 features, detections are merged
% if they are too close to each others, too short and too long
% detections are removed.
%
% Input
%   - detInfoTS            : matrix [nSamples x nFeatures] of the detection features
%   - DEF_a7               : structure of a7 detection settings
%       - absSigPow_Th      : thresholds vector for abs sigma power
%       - relSigPow_Th      : thresholds vector for rel sigma power
%       - sigCov_Th         : thresholds vector for sigma covariance
%       - sigCorr_Th        : thresholds vector for sigma correlation
%   - artifactDetectVector : vector of artifact in sample (tall vector)
%                            (1:artifact;  0:valid)
%   - PSDLowFreq           : tall vector of the energy in the low band of the slow ratio
%                            [nPSDWindow x 1] (info in PSA window, not in sample)
%   - PSDHighFreq          : tall vector of the energy in the high band of the slow ratio
%                            [nPSDWindow x 1] (info in PSA window, not in sample)
%
% Output 
%   - detectionVectFix : tall vector [ nSamples ] of the detection
%   - slowRInfoTS - slow ratio information vector (same number of datapoints as
%                   input dataVector)
%   - NREMClass : tall vector of the number of events detected
%               (0: means not "IN" the NREM spectral context)
%
% Author : Karine Lacourse, 
% Date : 2017-01-04
% Arrangement : Jacques Delfrate
% Date : 2018-02-13
% 
%--------------------------------------------------------------------------
   
    %----------------------------------------------------------------------
    % INIT
    %----------------------------------------------------------------------
    % Thresholds
    absSigmaThresh      = DEF_a7.absSigPow_Th;
    relSigmaThresh      = DEF_a7.relSigPow_Th;
    sigmaCovThresh      = DEF_a7.sigCov_Th;
    sigmaCorrThresh     = DEF_a7.sigCorr_Th;
    % Spindles
    maxSpindleLengthSec     = DEF_a7.maxDurSpindleSec;
    minSpindleLengthSec     = DEF_a7.minDurSpindleSec;
    % minimum duration in sec to merge spindles
    % Should be kept to 0
    minSpindleStepSec2Merge = 0;

    %----------------------------------------------------------------------
    % Script
    %----------------------------------------------------------------------
    nSamples        = length(detInfoTS(:,1));
       
    %------------------------------------------------------------------
    %% Find possible spindle
    %------------------------------------------------------------------
    possSpindleAbsSig   = detInfoTS(:,1) > absSigmaThresh;
    possSpindleRelSig   = detInfoTS(:,2) > relSigmaThresh;
    possSpindleSigCov   = detInfoTS(:,3) > sigmaCovThresh;
    possSpindleSigCorr  = detInfoTS(:,4) > sigmaCorrThresh;

    possSpindle  = possSpindleSigCorr & possSpindleSigCov & ...
        possSpindleAbsSig & possSpindleRelSig;

    % Possible spindles are turn off if there is an artifact during the
    % detection of a spindle    
    if DEF_a7.spindleNoArt==1
        possSpindle = possSpindle & ~artifactDetectVector;
    end 

    % Fixed the length of the spindle based on absSigmaPow and sigmaCov
    % Create a list of events
    [startEventSmp, endsEventSmp, ~] = event_StartsEndsDurations(possSpindle); 

    % Create the vector for possible length
    possLength = possSpindleAbsSig & possSpindleSigCov;      

    % For every events extend the length based on the absolute features
    for iEvt = 1 : length(startEventSmp)

        % Adjust the start of the detection
        i = 1;
        % While the asbSigma overcomes the threshold the start of the 
        % detection is extended 
        % (check to not ask for possSpindleAbsSig(0))
        while ( (startEventSmp(iEvt)-1 > 0) && ...
                (possLength(startEventSmp(iEvt)-1)==1) )
            startEventSmp(iEvt) = startEventSmp(iEvt)-1;
            i = i +1;
        end

        % Adjust the end of the detection
        i = 1;
        % While the asbSigma overcomes the threshold the end of the 
        % detection is extended
        % (check to not ask for possSpindleAbsSig outside the matrix)
        while (endsEventSmp(iEvt) < length(possLength)) &&...
                (possLength(endsEventSmp(iEvt)+1)==1)
            endsEventSmp(iEvt) = endsEventSmp(iEvt)+1;  
            i = i +1;                
        end

    end    

    % Remove any duplicated spindles
    eventListStartEnd = [startEventSmp, endsEventSmp];
    eventListStartEndUnique = unique(eventListStartEnd,'rows');
    startEventSmp = eventListStartEndUnique(:,1); 
    endsEventSmp = eventListStartEndUnique(:,2); 

    % Create duration in sample tab
    durationsEventSmp = endsEventSmp - startEventSmp + 1;

    % If there is at least one event
    if ~isempty(startEventSmp)

        if minSpindleStepSec2Merge>0
            % Fix spindle list by merging any too close spinldles
            eventLst2Merge  = ((startEventSmp(2:end) - endsEventSmp(1:end-1))...
                / DEF_a7.standard_sampleRate) <= minSpindleStepSec2Merge;
            % Add an event at the end to analyze all the set of events
            % Even if the last event has to be merge, it works
            eventLst2Merge = [eventLst2Merge; 0];
            % Verify the duration of the event to merge (must be shorter
            % than the minimum duration)
            iEventLst2Merge = find(eventLst2Merge==1);
            for iEvt2Merge = 1 : length(iEventLst2Merge)
                durSecFstEvt2merge = durationsEventSmp(iEventLst2Merge(iEvt2Merge))...
                    / DEF_a7.standard_sampleRate;
                durSecSecEvt2merge = durationsEventSmp(iEventLst2Merge(iEvt2Merge)+1)...
                    / DEF_a7.standard_sampleRate;
                if ( durSecFstEvt2merge < minSpindleLengthSec && ...
                        durSecSecEvt2merge < minSpindleLengthSec )
                    eventLst2Merge(iEventLst2Merge(iEvt2Merge))=1;
                    fprintf('event %i : 2 spindles are merged\n\n', ...
            iEvt2Merge);
                else
                    eventLst2Merge(iEventLst2Merge(iEvt2Merge))=0;
                end
            end
            iEventLstOK = find(eventLst2Merge==0);            

            iEvntKept = 0;
            updateEnd = 0;
            startEventSmpMerge      = zeros(length(iEventLstOK),1);
            endsEventSmpMerge       = zeros(length(iEventLstOK),1);
            durationsEventSmpMerge  = zeros(length(iEventLstOK),1);
            for iEvt = 1 : length(eventLst2Merge)
                % No modif, event does not have to be merge
                if eventLst2Merge(iEvt)==0 && updateEnd==0
                    iEvntKept = iEvntKept +1;
                    startEventSmpMerge(iEvntKept)       = startEventSmp(iEvt);
                    endsEventSmpMerge(iEvntKept)        = endsEventSmp(iEvt);
                    durationsEventSmpMerge(iEvntKept)   = durationsEventSmp(iEvt);
                % Merge (completed) : end and duration have to be updated 
                % It is already done for the start
                elseif eventLst2Merge(iEvt)==0 && updateEnd==1
                    updateEnd = 0; 
                    endsEventSmpMerge(iEvntKept)  = endsEventSmp(iEvt);
                    durationsEventSmpMerge(iEvntKept) = ...
                        endsEventSmpMerge(iEvntKept)-startEventSmpMerge(iEvntKept);
                % Merge : end is removed
                else 
                    % First event of the merge of events
                    if updateEnd==0
                        iEvntKept = iEvntKept +1;
                        startEventSmpMerge(iEvntKept) = startEventSmp(iEvt);              
                    end
                    % We accumulate every events that have to be merge until the
                    % last to know the final end and duration
                    % The start does not have to be updated
                    updateEnd = 1; % It is not more the first event of a merge of events    
                end
            end
        else
            startEventSmpMerge = startEventSmp;
            endsEventSmpMerge = endsEventSmp;
            durationsEventSmpMerge = durationsEventSmp;
        end

        %------------------------------------------------------------------
        %% Fix spindle list by removing too short spindle
        %------------------------------------------------------------------
        iSS2Keep = find(durationsEventSmpMerge >= (minSpindleLengthSec * DEF_a7.standard_sampleRate));
        startsEventSmpFixLength = startEventSmpMerge(iSS2Keep);
        endsEventSmpFixLength   = endsEventSmpMerge(iSS2Keep);
        durEventSmpFixLength    = durationsEventSmpMerge(iSS2Keep);

        %------------------------------------------------------------------
        %% Fix spindle list by removing too long spindle
        %------------------------------------------------------------------
        iSS2Keep = find(durEventSmpFixLength <= (maxSpindleLengthSec * DEF_a7.standard_sampleRate));
        startsEventSmpFixLength = startsEventSmpFixLength(iSS2Keep);
        endsEventSmpFixLength   = endsEventSmpFixLength(iSS2Keep); 

        % Write a warning if no spindle is found
        if isempty(startsEventSmpFixLength)
            warning('No spindles detected');   
        end

        % Go back to the detection vector in order to evaluate performance
        detectionVectFix = eventSmpList2DetectVect( [...
            startsEventSmpFixLength,endsEventSmpFixLength], size(detInfoTS,1) );

        %------------------------------------------------------------------
        %% Spindle context classifier
        %------------------------------------------------------------------        
        % if PSDLowFreq empty --> context classifier Off
        if ~isempty(PSDLowFreq)
            % Keep track of IN or OUT the NREM spectral context
            [~, slowRValidTS, slowRInfoTS] = ...
                a7subTurnOffDetSlowRatio(possSpindle, PSDLowFreq, ...
                PSDHighFreq, artifactDetectVector, DEF_a7);
            % Verify for each detection the NREM classifier flag
            NREMClass = zeros(length(startsEventSmpFixLength),1);
            for iEvt = 1 : length(startsEventSmpFixLength)
                % The context is considered "IN" if at least one sample of
                % the spindle is "IN" the NREM spectral context
                NREMClass(iEvt,1) = any( slowRValidTS(startsEventSmpFixLength(iEvt) ...
                    : endsEventSmpFixLength(iEvt)));
            end             
        else
            slowRInfoTS  = [];
        end            


    else
        % No detections possible
        detectionVectFix = zeros(nSamples,1);
        slowRInfoTS  = [];
        % Write a warning if no spindle is found
        warning('No spindles detected');        
    end


end

