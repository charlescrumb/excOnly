function writeNetworkScoresInCSV(filename,scores,networkID,separator,header)
%writeNetworkScoresInCSV(filename,scores,networkID,separator,header)
% will write the scores in Kaggle format as a CSV file when given the network ID.  
% scores should be an nxn matrix, whose (i,j)th entry indicates the strength of 
% connection from neuron i to neuron j.

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Authors: Mehreen Saeed
% Date: Jan 2014
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

if nargin<4
    separator=',';
end

if nargin<5
    %construct the header
    header = {'NET_neuronI_neuronJ', 'Strength'};
end

%check the file extension
fileExtension = filename(end-3:end);

if ~strcmp(fileExtension, '.csv'), 
    filename=[filename '.csv']; 
end
%----------------------------------------------------------------------
if ~exist(filename, 'file')
    fp=fopen(filename,'w');

    %write the header
    if ~isempty(header)
        for k=1:length(header)-1
            fprintf(fp, '%s%s', header{k}, separator);
        end
        fprintf(fp, '%s\r\n', header{end});
    end
else
	fp=fopen(filename,'a');
end

%write the rest of the scores
for i=1:size(scores, 1)
    for j=1:size(scores, 2)
        firstColumn = [networkID '_' num2str(i) '_' num2str(j) ];
        fprintf(fp,'%s%s',firstColumn,separator);
        fprintf(fp,'%f\n', full(scores(i, j)));
    end
end

fclose(fp);