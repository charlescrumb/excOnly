function concatenateFiles(filenames, outfilename)
%concatenateFiles(filenames, outfilename)
% Concatenate all the files in filenames, removing duplicate headers.
% The result is in outfilename.
% Goal: make a suitable submission.

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Author: Isabelle Guyon 
% Date: Jan 2014
% Last modified: 01/30/2014
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

copyfile(filenames{1}, outfilename);
fout=fopen(outfilename, 'a');
fprintf('Concatenating files...\n');
tic
for k=2:length(filenames)
    fid=fopen(filenames{k});
    tline = fgetl(fid); %get rid of header
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        fprintf(fout, '%s\n', tline);
    end
    fclose(fid);
end
fclose(fout);
toc
fprintf('...done\n');