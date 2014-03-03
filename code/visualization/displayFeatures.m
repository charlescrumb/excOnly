function displayFeatures(F,D,neuron_i,neuron_j,G,network,h)
%%%displayFeatures(F,D,neuron_i,neuron_j,G,network,h)
%%%Helper function for the challenge browser
%%this will first compute the features for a neuron pair (i,j) given their
%%indices and then display them
%%% D is the discretized fluorescence signal
%%% G is the global conditioning level
%%% h is the figure handle where the features are to be displayed
%%% network contains the actual i,j scores
%%% the truth value

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Authors: Mehreen Saeed
% Date: Jan 2014
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

x=600;
y=600;
width = 240;
height = 30;
nxt = 240;
nxty = 30;
%% display the truth values if network is given 
if (length(network))
    truth_ij = network(neuron_i,neuron_j);
    truth_ji = network(neuron_j,neuron_i);

    truthstr1 = uicontrol('Parent',h,'Style','text','Position', [x y width height], 'FontSize', 10,...
                        'String',['truth ' num2str(neuron_i) '->' num2str(neuron_j) ' = ' num2str(truth_ij)]);
    truthstr2 = uicontrol('Parent',h,'Style','text','Position', [x+nxt y width height], 'FontSize', 10,...
                        'String',['truth ' num2str(neuron_j) '->' num2str(neuron_i) ' = ' num2str(truth_ji)]);
else
    truthstr1 = uicontrol('Parent',h,'Style','text','Position', [x y width height], 'FontSize', 10,...
                        'String',['         ']);
    truthstr2 = uicontrol('Parent',h,'Style','text','Position', [x+nxt y width height], 'FontSize', 10,...
                        'String',['         ']);
end

%% Calculate the GTE from the joint PDF
%% Calculate the joint PDF
P = calculateJointPDFforGTE(D(:,[neuron_i neuron_j]), G,'debug',false);
GTE = calculateGTEfromJointPDF(P,'debug','false');

%% display GTE
GTEstr1 = uicontrol('Parent',h,'Style','text','Position', [x y-nxty width height], 'FontSize', 10,...
                    'String',['GTE (' num2str(neuron_i) ',' num2str(neuron_j) ') = ' num2str(GTE(1,2),'%3.3e')]);
GTEstr2 = uicontrol('Parent',h,'Style','text','Position', [x+nxt y-nxty width height], 'FontSize', 10,...
                    'String',['GTE (' num2str(neuron_j) ',' num2str(neuron_i) ') = ' num2str(GTE(2,1),'%3.3e')]);
                
%% Compute Cross Correlation 
xcorre = computeCrossCorrelation(F(:,[neuron_i neuron_j]));

%% display cross correlation
xcorrstr1 = uicontrol('Parent',h,'Style','text','Position', [x y-nxty*2 width height], 'FontSize', 10,...
                    'String',['xcorr(' num2str(neuron_i) ',' num2str(neuron_j) ') = ' num2str(xcorre(1,2),'%3.3e')]);
xcorrstr2 = uicontrol('Parent',h,'Style','text','Position', [x+nxt y-nxty*2 width height], 'FontSize', 10,...
                    'String',['xcorr(' num2str(neuron_j) ',' num2str(neuron_i) ') = ' num2str(xcorre(2,1),'%3.3e')]);
                
%% Compute Information Gain Entropy 
IGentropy = computeIGEntropy(F(:,[neuron_i neuron_j]));

%% display Information Gain Entropy
IGentropystr1 = uicontrol('Parent',h,'Style','text','Position', [x y-nxty*3 width height], 'FontSize', 10,...
                    'String',['IGEntropy(' num2str(neuron_i) ',' num2str(neuron_j) ') = ' num2str(IGentropy(1,2),'%3.3e')]);
IGentropystr2 = uicontrol('Parent',h,'Style','text','Position', [x+nxt y-nxty*3 width height], 'FontSize', 10,...
                    'String',['IGEntropy(' num2str(neuron_j) ',' num2str(neuron_i) ') = ' num2str(IGentropy(2,1),'%3.3e')]);

%% Compute Information Gain Gini Index 
IGgini = computeIGGini(F(:,[neuron_i neuron_j]));

%% display Information Gain Gini Index
IGginistr1 = uicontrol('Parent',h,'Style','text','Position', [x y-nxty*4 width height], 'FontSize', 10,...
                    'String',['IGGini(' num2str(neuron_i) ',' num2str(neuron_j) ') = ' num2str(IGgini(1,2),'%3.3e')]);
IGginistr2 = uicontrol('Parent',h,'Style','text','Position', [x+nxt y-nxty*4 width height], 'FontSize', 10,...
                    'String',['IGGini(' num2str(neuron_j) ',' num2str(neuron_i) ') = ' num2str(IGgini(2,1),'%3.3e')]);
