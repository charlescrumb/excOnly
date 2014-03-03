function offset = getOffset(offset)
%offset = getOffset(offset)

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Authors: Mehreen Saeed
% Date: Jan 2014
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

prompt = {'Enter offset value:'};
dlg_title = 'Input offset';
num_lines = 1;
defaultanswer = {'0.2'}
answer = inputdlg(prompt,dlg_title,num_lines,defaultanswer);

if isempty(answer)
    return;
end

offset = str2num(answer{1});
%%check if its valid
if (offset<0|offset>4)
    return;
end

