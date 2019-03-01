function cell2tab(filename, cellArray, permission)
% cell2tab.m
%
% Purpose:
% Writes 2D cell array content into a tab delimited file. 
% Tab delim is default.  All cell contents are written as strings.
% 
% Usage:
% cell2tab(filename, cellArray, permission)
% 
% Variables:
% filename     : Name of the file to save. [ i.e. 'text.txt' ]
% cellArray    : Name of the cell to write
% permission   : fopen permission 'w' is default (overwrite); change to 'a' for append
%
% Examples:
% % append WORK.data cell to the text.txt file
% cell2tab('text.txt.', WORK.data, 'a')
% % output the list of filenames; list is transposed vertically using []'
% cell2tab(fullfile(WORK.outputDir, 'edf.list.txt'), cellstr([WORK.inputFiles]')) ;
%
% Authors:
% Simon Warby 2012-11-23
% 
% Changelog:
% 2014-04-17    Fixed tab delimited output; removed delimeter option


if nargin<3
    permission = 'w';
end
 

FID = fopen(filename,permission);
for z=1:size(cellArray,1)
    for s=1:size(cellArray,2)
        
        var = eval(['cellArray{z,s}']);
        
        if size(var,1) == 0
            var = '';
        end
        
        if iscell(var)
            var = cell2mat(var);
        end
        
        if isnumeric(var)
            var = num2str(var);
        end
                
        fprintf(FID,'%s\t', var);  % sting + tab delim
        
    end
    fprintf(FID,'\r\n');  % windows line ending

end
    fclose(FID);
    
