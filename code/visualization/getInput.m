function [neuron_i neuron_j] = getInput(maxValue)
%[neuron_i neuron_j] = getInput(maxValue)

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Authors: Mehreen Saeed
% Date: Jan 2014
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

neuron_i=1;
neuron_j=2;

prompt = {'Enter i value:','Enter j value:'};
dlg_title = 'Input neuron (i,j) pair';
num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines);

if isempty(answer)
    return;
end

neuron_i = str2num(answer{1});
if (neuron_i > maxValue || neuron_i < 1)
    neuron_i=1;
end

neuron_j = str2num(answer{2});
if (neuron_j > maxValue || neuron_j < 1)
    neuron_j=1;
end

