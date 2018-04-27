
function [output, errorMess] = butterFiltZPLowPassFiltFilt(input, Fpass, Fs, filterOrder)
% Lowpass Butterworth IIR filter.
% Filter constructed using butterworth
% Filter constructed with zero-pole-gain from  
% converted into a second order section (SOS)
% and the non-linear phase is corrected by the filtfilt.
%
% Input: 
%   input is an EEG timeseries, tall vector
%   Fpass : cut off frequency (last frequency to pass)
%   Fs : frequency sampling
%   (optional) filterOrder : the order of the filter (default is 10)
%
% Output:
%   output : the filtered EEG signal, tall vector.
%   errorMess : cell of string, empty if everyting is ok
%
% Author     : Simon Warby      2014-04-30 
%              (from filtfilt_butter_lowpass.m) 
% Arrangment : Karine Lacourse  2015-10-13
%              Jacques Delfrate 2018-02-13
%-------------------------------------------------------------------------

    if nargin < 4
        filterOrder = 10; % the higher the order, the steeper the slope of the filter and longer processing
    end

    nyq             = Fs*0.5;   % Nyquist frequency
    Wn              = Fpass/nyq;% Window; fraction of 0-nyq being filtered
    showFreqResp    = 0;        % To show the frequency response of the filter

    errorMess = cell(0,0);

    % Check that input is numeric
    if ~isnumeric(input)
        errorMess{length(errorMess)+1} = 'lowpass - input is not numeric'; 
    end

    % Check that Fpass is not less than zero or greater than, or too close to nyquist
    if Fpass < 0
        errorMess{length(errorMess)+1} = 'lowpass - cutoff freq must be > 0 Hz';
    end

    if Fpass>=nyq
        errorMess{length(errorMess)+1} = 'lowpass - frequency cannot be > fs/2';
    end

    % Check that the Fs is not too low for the Fpass   
    if Wn>=1
        errorMess{length(errorMess)+1} = 'lowpass - frequency is higher than Nyquist';      
    end
    % If an error occurs dont apply the filter (xml file has to be fixed)
    if isempty(errorMess)
        
        % Design the filter with zero pole
        [z,p,k] = butter(filterOrder, Wn, 'low');
        
        % Convert it into the second order section(SOS)
        % Note : Embedding the gain in the first section when scaling a 
        % direct-form II structure is not recommended and may result in 
        % erratic scaling. To avoid embedding the gain, use ss2sos with two outputs.          
        [sos,scaleFactG] = zp2sos(z,p,k);
        
        % Warning on the scaling factor, because the magnitude response
        % will have a gain, but it does not affect the filtered time
        % series. I guess that the gain is known and the time series is
        % corrected.  
        if showFreqResp==1
            if scaleFactG<0.9
                warning('The scaling factor is %f',scaleFactG);
            end
        end
        
        % Plot Freq response
        if showFreqResp == 1
            fvtool(sos,'Fs',Fs); xlim([0 nyq/2]);
        end
        
        % Apply filtering with FilFilt to remove the phase and group delay
        % It is applied 2 times on the eeg signal (forward and backward)
        output = filtfilt(sos,scaleFactG,input);
    else
        output = input;
    end

end


