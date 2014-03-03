function samples = getSampleInput(maxValue)
%samples = getSampleInput(maxValue)

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Authors: Mehreen Saeed
% Date: Jan 2014
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

samples = 1:1000;

prompt = {'Enter start value:','Enter end value:','Enter interval'};
dlg_title = 'Input samples to plot';
num_lines = 1;
defaults = {'1','1000','1'};
answer = inputdlg(prompt,dlg_title,num_lines,defaults);

if isempty(answer)
    return;
end

starting = str2num(answer{1});
ending = str2num(answer{2});
interval = str2num(answer{3});

%make sure we have a valid input
if starting > 0 & starting <= maxValue & ending >0 & ending <=maxValue ...
        & ending > starting & interval < maxValue
    samples = starting:interval:ending;
end


