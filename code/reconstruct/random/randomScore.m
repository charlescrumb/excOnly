function score = ramdomScore(F, arg)
%score = ramdomScore(F)
% Baseline method to compute scores assigning random values.

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Authors: Isabelle Guyon
% Date: Jan 2014
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

N = size(F, 2);
score=sprand(N,N,0.1);
