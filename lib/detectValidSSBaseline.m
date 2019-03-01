function detectValidSSBaseline(sleepStageVect, DEF_a7)
%detectValidSSBaseline display warning for sleep staging
%  If a sleep stage baseline defined in the initA7_DEF is not member of the
%  input sleep stage vector a warning is displayed
%
% Input 
%      sleepStageVect : input sleep stage vector - tall vector (cell or tab of char)
%      DEF_a7         : structure of a7 settings
%
% Author : Jacques Delfrate 2018-08-13

    % Make sure the bslSleepStaging is a cell of char
    DEF_a7.bslSleepStaging = cellstr(DEF_a7.bslSleepStaging);
    % Look for every stage included in the baseline
    for iBsl=1:length(DEF_a7.bslSleepStaging)
        valideBsl = unique(ismember(sleepStageVect, DEF_a7.bslSleepStaging{iBsl}));
        if ~any(valideBsl)
            warning('%s is not member of the input sleep stage vector \n', ...
                DEF_a7.bslSleepStaging{iBsl});
        end
    end

end

