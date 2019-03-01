function [PSD, freqBins] = calcPSDNormType(dataVector, DEF_a7, windowScalingType,...
    removeMean, zeroPad, normType)

% calcPSD.m
% 
% Purpose:
% calcPSD calculates the Power Spectral Density (PSD) of a time signal. 
% 4 different normalization can be used 'NoNorm', 'ScalFact', 'EnergyCons'
% and 'IntegSpectPow'.
%
% Input: 
%   dataVector        : input timeseries, should be a tall vector.
%   DEF_a7            : structure of a7 settings
%   windowScalingType : window function ie 'hann', 'rectwin'
%   removeMean:       : To remove the mean (the DC value is written at the
%                       freq=0, but it is not taken into account in the
%                       normalization
%   zeroPad:          : To pad data with zeros to a length factor of 2 (usually for speed).
%                       The resolution of the FFT is inscreased in that
%                       case.
%   normType:         : normalization method to use: 
%                       'NoNorm', 'ScalFact', 'EnergyCons', 'IntegSpectPow'
%
% Output:
%   PSD:                output matrix of power spectral density, 
%                       rows are each FFT window, cols are frequency bins
%   freqBins:           single row indicating the header of the frequency bin
%
%  
% Authors:
% Modified from calcPSD scripts
% from Emil Munk and Hyatt Moore 2013
% Simon Warby 2014-06-11
% Karine lacourse 2015-05
% Jacques Delfrate 2018-02

    if DEF_a7.WinStepSec == 0
        error('PSD computation: the overlap of the fft window is 100%, a minimum slide of 1 sample is needed');
    end

    %--------------------------------------------------------------------------
    % General init before the FFT computation
    %--------------------------------------------------------------------------
    if nargin==5
        if zeroPad==0
            nFFTSamples = round( DEF_a7.winLengthSec*DEF_a7.standard_sampleRate);
        else
            % or use next highest power of 2 greater than or equal to length(x) to calculate FFT.
            nFFTSamples = 2^(nextpow2( DEF_a7.winLengthSec*DEF_a7.standard_sampleRate));
        end
    elseif nargin==6
        if zeroPad==0
            nFFTSamples = round( DEF_a7.winLengthSec*DEF_a7.standard_sampleRate);
        else 
            % or use the specified length in the input argument
            nFFTSamples = round( DEF_a7.ZeroPadSec*DEF_a7.standard_sampleRate);
        end        
    else
        error('calcPSDNormType needs 8 or 9 input arguments');
    end
    
    dataLength_sec = length(dataVector)/DEF_a7.standard_sampleRate ;   % total length of the timeseries; in seconds
    nFFTWind = floor((dataLength_sec- DEF_a7.winLengthSec)/DEF_a7.WinStepSec)+1 ;  % number of fft windows needed

    % Number of samples from the time seires included in the FFT window
    nDataSampleInFFTWind = round( DEF_a7.winLengthSec * DEF_a7.standard_sampleRate);
    % Window Scaling Function - weight the values in the analysis window (ie hann)
    wt = (0:nDataSampleInFFTWind-1)/(nDataSampleInFFTWind-1);
    % Some scaling windows are not computed with the matlab function "eval"
    % in order to be able to modfify the factor
    if strcmp(windowScalingType,'blackman')
        windowScalingWeights = 0.4266 - 0.4965*cos(2*pi*wt) + 0.0768*cos(4*pi*wt);
        windowScalingWeights = windowScalingWeights';
    elseif strcmp(windowScalingType,'hann')
        defaultFact = 0.5; %  defaultFact = 0.8; % less attenuation
        windowScalingWeights = defaultFact - defaultFact*cos(2*pi*wt);  % Hann
        windowScalingWeights = windowScalingWeights';
    else
        windowScalingWeights = eval([windowScalingType '(' num2str(nDataSampleInFFTWind) ')']) ;
    end                                  


    % Coherent gain
    CG = sum(windowScalingWeights)/length(windowScalingWeights);
    % Noise gain
    NG = sum(windowScalingWeights.^2)/length(windowScalingWeights); 

    % Calculate the number of frequency bins (= number of cols in the output)
    if ( ~ rem(nFFTSamples,2) )    % one-sided, Nfft is even
        numFreqBins = round(nFFTSamples/2+1);
    else                            % one-sided, Nfft is odd
        numFreqBins = round((nFFTSamples+1)/2);
    end

    % Create vector with evenly spaced frequency bins, numFreqBins = numCols in the output
    freqBins = (0:numFreqBins-1)*DEF_a7.standard_sampleRate/nFFTSamples ;


    % The zero padding normalization factor
    % The zeros padding increases the frequency resolution (more bins) 
    % for the same sampling rate. The FFT output un-normalized are the same
    % in average but not for the sum.
    %
    % zeroPadNorm only exist to view the FFT with or without zeros padding
    % supperposed. The scale will be kept even if the sum of the energy is
    % not the same.
    zeroPadNorm         = nDataSampleInFFTWind/nFFTSamples;    

    % Pre-allocate the matrix for the PSD data (rows = results from each FFT window, columns = frequencyBin) 
    PSD = zeros(nFFTWind, numFreqBins); 

    % loop over the sliding FFT windows
    for iA = 1:nFFTWind

        % Extract the FFT window from data, the start position moves with each loop 
        start = round((iA-1)*DEF_a7.WinStepSec*DEF_a7.standard_sampleRate);
        dataWindow = dataVector(start+1 : start + round( DEF_a7.winLengthSec*DEF_a7.standard_sampleRate)) ;  

        %----------------------------------------------------------------------
        % (1) Remove the mean of the signal, dc offset
        %----------------------------------------------------------------------
        if(removeMean)
            dataWindow_mean = mean(dataWindow) ;
            dataWindow = dataWindow-dataWindow_mean ;
        end

        %----------------------------------------------------------------------
        % (2) Apply the window scaling function (ie hanning window)
        %----------------------------------------------------------------------
        % The window shape the data and affects the PSD
        % The zeros padding shapes the data and affects the PSD
        % The scale of the fft is kept if the padding is added after the hann
        dataWindowScaled = dataWindow .* windowScalingWeights ;

        %----------------------------------------------------------------------
        % (3) Pad the data after applying the scaling window
        %----------------------------------------------------------------------
        if(zeroPad)
            dataWindowZeroPad = zeros(nFFTSamples,1);
            nZeroPad = round(nFFTSamples- DEF_a7.winLengthSec*DEF_a7.standard_sampleRate);
            dataWindowZeroPad( round(nZeroPad/2) + 1 : round(nZeroPad/2) + ...
                length(dataWindowScaled) ) = dataWindowScaled;
            dataWindowScaled = dataWindowZeroPad;
        end

        %----------------------------------------------------------------------
        % (4) Perform the fft
        %----------------------------------------------------------------------
        % Apply the correcting factor from the scaling window
        if strcmp(normType,'ScalFact')
            dataWindowScaled = dataWindowScaled / CG;
        end


        % Even with a dyadic length for the FFT window, the length of the
        % dataWindowScaled is always the same than the nFFTSamples
        % The zero padding is added around the data.   
        fft_dataWindow = fft(dataWindowScaled, nFFTSamples);  
        FFTModule = abs(fft_dataWindow);            % Calculate the PSD (un-normalized)



        if strcmp(normType,'EnergyCons')
            %----------------------------------------------------------------------
            % (5) fft is symmetric, throw away second half of psd, which are redundant
            %
            % Frequency resolution stays the same through different sampling rate,
            % but the length is longer for the higher sampling rate to reach the nyquist frequency.
            %
            % The FFT is divided by the (number of samples)^2 in the single side FFT
            % to make the scale of the FFT independant of the sampling rate.
            % Then the sum of the energy of the FFT is the same for a sampling 
            % rate of 100 Hz or 256 Hz.
            %
            %----------------------------------------------------------------------
            FFTModule2              = FFTModule.^2; 
            FFTModule2_1Sided(1)    = FFTModule2(1);
            FFTModule2_1Sided(2:numFreqBins) = 2 * FFTModule2(2:numFreqBins);
            FFTSingleSide           = FFTModule2_1Sided/((nFFTSamples/2).^2);  

            %----------------------------------------------------------------------
            % (6) Normalization to satisfy Parseval's theorem of conservation of power
            %----------------------------------------------------------------------
            % (a) The average energy of the original signal, before window scaling and padding 
            dataWinFFT_EngAvg   = sum(dataWindow.^2/length(dataWindow));

            % (b) The energy of the PSD (only normalized for the sampling rate)
            PSD_EngAvg          = sum(FFTSingleSide);


            % (d) The Factor to normalize the energy
            % (e) To normalize the energy of the single sided FFT
            % The energy of the time series is zero (flat line)
            if dataWinFFT_EngAvg < eps
                PSDNormalized       = zeros(size(FFTSingleSide));
            % The energy of the time series is valid
            else
                power2NormFactor    = PSD_EngAvg*zeroPadNorm/dataWinFFT_EngAvg ; 
                % (e) To normalize the energy of the single sided FFT
                PSDNormalized       = FFTSingleSide/power2NormFactor;              
            end   

        elseif strcmp(normType,'IntegSpectPow')

            FFTModule2      = (FFTModule/nFFTSamples).^2;
            FFTModule2_1Sided(1) = FFTModule2(1);
            FFTModule2_1Sided(2:numFreqBins) = 2 * FFTModule2(2:numFreqBins);
            PSDNormalized = FFTModule2_1Sided/NG;

        elseif strcmp(normType,'EEGLAB')

            FFTModule2      = FFTModule.^2;
            FFTModule2_1Sided(1) = FFTModule2(1);
            FFTModule2_1Sided(2:numFreqBins) = 2 * FFTModule2(2:numFreqBins);
            PSDNormalized = FFTModule2_1Sided/(NG * nFFTSamples * DEF_a7.standard_sampleRate);      


        elseif strcmp(normType,'NoNorm') || strcmp(normType,'ScalFact')
            % zeroPadNorm only exist to view the FFT with or without zeros padding
            % supperposed. The scale will be kept even if the sum of the energy is
            % not the same.
            FFTModule2              = (FFTModule/nFFTSamples).^2; 
            FFTModule2_1Sided(1)    = FFTModule2(1);
            FFTModule2_1Sided(2:numFreqBins) = 2 * FFTModule2(2:numFreqBins);        
            PSDNormalized           = FFTModule2_1Sided / zeroPadNorm.^2;
        else
            error('Unsupported PSD normalization method : %s', normType);
        end

        % Add the fft results (row) to the psd matrix
        PSD(iA,:) = PSDNormalized';

        % The first frequency is the DC offset
        if(removeMean)
            PSD(iA,1) = dataWindow_mean;
        end

    end
end


