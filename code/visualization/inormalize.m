function Y = inormalize(X)
%Y = inormalize(X)
% Inputs:
% X -- matrix of reals.
% Retruns:
% Y -- matrix X normalized between 0 and 1.
% If all values are equal, return all ones.

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Authors: Isabelle Guyon
% Date: April 2002
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

mini = min(min(X));
maxi = max(max(X));
delta = maxi-mini;
if delta==0,
   [n, p] = size(X);
   Y = ones(n, p);
else
   Y = (X - mini)/delta;
end