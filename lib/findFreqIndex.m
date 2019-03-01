function [ iFreqStart, iFreqStop ] = findFreqIndex( lowFreq, highFreq, freqBins )
%
% Function to find the nearest index in a frequency vector that corresponds to a frequency band.
%
% Input
%      lowFreq  : low frequency cut off for band pass filter
%      highFreq : high frequency cut off for band pass filter
%      freqBins : single row indicating the header of the frequency bin
%
% Output
%      iFreqStart : start index for frequency band
%      iFreqStop  : end index for frequency band
%      
% Author : Karine Lacourse
%
% Arrangement : Jacques Delfrate
%

    % LOW FREQ
    iFreqStart  = find(freqBins >=  lowFreq, 1, 'first');
    % The freqBin are always lower than the lowFreq
    if isempty(iFreqStart)
        iFreqStart = length(freqBins);
    end
    % Check if the previous one is closer
    if iFreqStart > 1
        [val, indice] = min([abs(freqBins(iFreqStart)-lowFreq),...
            abs(freqBins(iFreqStart-1)-lowFreq)]);
        iFreqStart = iFreqStart-indice+1;
    end

    % HIGH FREQ
    iFreqStop   = find(freqBins >=  highFreq, 1, 'first');
    % The freqBin are always lower than the highFreq
    if isempty(iFreqStop)
        iFreqStop = length(freqBins);
    end    
    % Check if the previous one is closer
    if iFreqStop > 1
        [val, indice] = min([abs(freqBins(iFreqStop)-highFreq),...
            abs(freqBins(iFreqStop-1)-highFreq)]);
        iFreqStop = iFreqStop-indice+1;
    end    

    if iFreqStart == iFreqStop
        error('The selected frequency band is excluded of the frequency bins');
    end

end

