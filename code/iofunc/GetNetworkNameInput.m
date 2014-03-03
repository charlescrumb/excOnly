function [F, D, G, network, totalSamples, totalNeurons, samples] = GetNetworkNameInput(names,F,D,G, network,dataDirectory,h,selection)
%%% [F, D, G, network, totalSamples, totalNeurons, samples] = GetNetworkNameInput(names,F,D,G, network,dataDirectory,h,openDlg)
%%% Function to input network name from a given list and then read the
%%% corresponding network and return fluorescence signal F, discretized
%%% signal D, target values of network (if file exists) and index of hte
%%% network name selected by the user.

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Author: Mehreen Saee 
% Date: Feb 2014
% Last modified: Isabelle Guyon 02/02/2014
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

if nargin<8
    [selection ok] = listdlg('PromptString','Select a network',...
                            'SelectionMode','single','ListString', names);
else
    ok = 1;
end

if ~ok
    [totalSamples totalNeurons] = size(F);
    samples = 1:1000;
    return;
end

extension = '.txt';

%% define file name
networkId = names{selection};
fluorescenceFile = [dataDirectory filesep 'fluorescence_' networkId extension];
%% Load the Fluorescence signal
F = load_data(fluorescenceFile);
%% name of the network file that contains truth values (may not be present)
% (network architecture provided only for training data to the participants)
networkFile = [dataDirectory filesep 'network_' networkId extension];

%% load the truth values if the file exists
if exist(networkFile,'file')
    fprintf('Loading %s\n', networkFile);
    network = readNetworkScores(networkFile);
else
    network = [];
end
%% Discretize the fluorescence signal
fprintf('Discretizing the fluorescence signal\n');
[D, G] = discretizeFluorescenceSignal(F,  'conditioningLevel', 0.25, 'bins',3,'debug',false);
index = selection;
%% size of data and set sample value
[totalSamples totalNeurons] = size(F);
samples = 1:1000;
if (totalSamples<1000)
    samples = 1:totalSamples;
end
%%update the title of the figure
if ~isempty(h)
    set(h,'Name', ['Connectomics Neuron Pair Browser for network : ' names{index}]);
end
