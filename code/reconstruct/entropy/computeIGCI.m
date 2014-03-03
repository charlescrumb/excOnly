function scores = computeIGCI(F, debug_)
% scores = computeIGCI(F, debug_)
% Baseline method to compute scores based on 
% Information Geometry Causal Inference,
% Inspired by:
% [1]  P. Daniušis, D. Janzing, J. Mooij, J. Zscheischler, B. Steudel,
%      K. Zhang, B. Schölkopf:  Inferring deterministic causal relations.
%      Proceedings of the 26th Annual Conference on Uncertainty in Artificial 
%      Intelligence (UAI-2010).  
%      http://event.cwi.nl/uai2010/papers/UAI2010_0121.pdf
% It boils down to computing the difference in entropy between pairs of
% variables:
% scores(i, j) = H(j) - H(i)

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Author: Isabelle Guyon
% Date: Jan 2014
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

if nargin<2,
    debug_ = false; 
end

%% Discretize the fluorescence signal
if debug_, fprintf('discretize...\n'); end
D = discretizeFluorescenceSignal(F,'bins',3);
if debug_,  fprintf('...done\n'); end

%% Compute the entropy
if debug_, fprintf('compute entropy...\n'); end
H = entropy(D);
if debug_,  fprintf('...done\n'); end

%% Compute the scores as entropy differences (vectorized :-))
if debug_, fprintf('compute entropy differences...\n'); end
n=length(H);
scores = H(ones(n,1),:);
scores = scores - scores';
if debug_,  fprintf('...done\n'); end
