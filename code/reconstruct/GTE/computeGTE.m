function GTE = computeGTE(F, debug_)
%GTE = computeGE(F)
%GTE = computeGE(F, debug_)
% Baseline method to compute scores based on 
%    (Stetter 2013) Stetter, O., Battaglia, D., Soriano, J. & Geisel, T. 
%    Model-free reconstruction of excitatory neuronal connectivity from 
%    calcium imaging signals. PLoS Comput Biol 8, e1002653 (2012).
% Set flag debug_ to true to visualize histograms.

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Authors: Javier Orlandi
% Date: Dec 2013
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

if nargin<2,
    debug_ = false; 
end

%% Discretize the fluorescence signal
[D, G] = discretizeFluorescenceSignal(F, 'debug', debug_, 'conditioningLevel', 0.25, 'bins', [-10,0.12,10]);

%% Calculate the joint PDF
P = calculateJointPDFforGTE(D, G);

%% Calculate the GTE from the joint PDF
GTE = calculateGTEfromJointPDF(P);
