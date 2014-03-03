function M=load_data(filename)
% M=load_data(filename)
% Load the data in filename as a Matlab matrix and saves it in Matlab
% format for faster re-load.

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Authors: Isabelle Guyon
% Date: Jan 2014
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

% Check the extension
ext=filename(end-3:end);
if ~strcmp(ext(1), '.'), 
    ext='.mat';
else
    filename=filename(1:end-4);
end

if exist([filename, '.mat'], 'file')
    ext = '.mat';
end

if strcmp(ext, '.txt') || strcmp(ext, '.csv') || strcmp(ext, '.mat')
    fprintf('Loading %s%s...', filename, ext);
    M=load([filename ext]);
    if(isstruct(M))
        M=M.M;
    end    
    fprintf(' done\n');
else
    fprintf('Extension %s not supported\n', ext);
    M=[];
    return;
end

if ~exist([filename, '.mat'], 'file')
    fprintf('Saving %s.mat...', filename);
    save(filename, 'M');
    fprintf(' done\n');
end

