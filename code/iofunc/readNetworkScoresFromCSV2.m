function [scores,header,networkID]=readNetworkScoresFromCSV2(filename,networkID,n)
%[scores,header,networkID]=readNetworkScoresFromCSV(filename,networkID,n)
% Read the scores from kaggle CSV format anf fills a  
% sparse "scores" nxn matrix, whose (i,j)th entry indicates the strength of 
% connection between neuron i and j.
% n is the dimension of the scores matrix....if it is not specified it will
% be inferred from the CSV file.

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Authors: Isabelle Guyon
% Date: Jan 2014
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

if nargin < 3, n=-1; end
if nargin < 2, networkID=''; end
scores=[]; header=[];

fileExtension = filename(end-3:end);

if ~(strcmp(fileExtension, '.csv') || strcmp(fileExtension, '.txt')), 
    filename=[filename '.csv']; 
end

if ~exist(filename, 'file')
    fprintf([filename ' not found']);
    return; 
end

% Get Header and number of columns
fid=fopen(filename,'r');
tline = fgetl(fid);
comma=[strfind(tline, ','), length(tline)+1];
header{1}=tline(1:comma(1)-1);
header{2}=tline(comma(1)+1:comma(2)-1);
header{3}=tline(comma(2)+1:end);

% Get network name
tline = fgetl(fid);
if isempty(networkID), 
    uns=strfind(tline, '_');
    networkID = tline(1:uns(end-1)-1);
end

% Get the matrix dimension
if n==-1
    while 1
        tl = fgetl(fid);
        if ~ischar(tl), break, end
        tline =tl;
    end
    uns=strfind(tline, '_');
    n = str2double(tline(uns(1)+1: uns(2)-1));
end
fclose(fid);

scores = sparse(n, n);

% Second pass
fid=fopen(filename,'r');
tline = fgetl(fid);
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    uns=strfind(tline, '_');
    net=tline(1:uns(1)-1);
    if strcmp(net, networkID)
        comma=[strfind(tline, ','), length(tline)+1];
        i = str2double(tline(uns(1)+1:uns(2)-1));
        j = str2double(tline(uns(2)+1:comma(1)-1));        
        strength=str2double(tline(comma(1)+1:comma(2)-1));
        if strength>0
            scores(i, j)=strength;
        end
    end
end

fclose(fid);