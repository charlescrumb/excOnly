function scores = computeMI(F, debug_)
% scores = computeMI(F, debug_)
% Baseline method to compute scores based on 
% Mutual Information
% scores(i, j) = H(j) + H(i) - H(i, j)
% Note: this should be the same as computeIGEntropy (except perhaps
% diagonal elements...)

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
tic
%% Compute the entropy
if debug_, fprintf('compute entropy...\n'); end
Hi = entropy(D);
if debug_,  fprintf('...done\n'); end

%% Compute the join entropies
if debug_, fprintf('compute joint entropies...\n'); end
Hij = jointEntropy(D);
if debug_,  fprintf('...done\n'); end

%% Compute the scores as H(j) + H(i) - H(i, j) (vectorized :-))
if debug_, fprintf('compute mutual information...\n'); end
n=length(Hi);
scores = Hi(ones(n,1),:);
scores = scores + scores' - Hij;
if debug_,  fprintf('...done\n'); end
toc

end
