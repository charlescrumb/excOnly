function H = entropy(X, dim)
%Z = entropy(X)
%Z = entropy(X, dim)
% Compute entropy H(X) of a discrete variable X.
% This works even if X has non integer discrete levels.
% If X is a matrix, compute the entropy of each column vector.
% if dim is given, compute the entropy of the columns (dim =1) or lines
% (dim =2).
% see also jointEntropy.

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Author: Mo Chen
% Date: 13 Mar 2012
% Last modified: 01/31/2014 by Isabelle Guyon
% Contact: mochen80@gmail.com
% Source: http://www.mathworks.com/matlabcentral/fileexchange/35625-information-theory-toolbox
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

[nlin, ncol]=size(X);
if nargin<2, 
    if ncol==1, 
        dim=1; 
    elseif nlin==1,
        dim=2;
    else
        dim=1;
    end
end

if dim==1 % The propose of not transposing is not to touch X to save memory
    H=zeros(1,ncol);
    for k=1:ncol
        [u,~,label] = unique(X(:,k));
        p = full(mean(sparse(1:nlin,label,1,nlin,numel(u),nlin),1));
        H(k) = -dot(p,log2(p+eps));
    end
else
    H=zeros(nlin,1);
    for k=1:nlin
        [u,~,label] = unique(X(k,:));
        p = full(mean(sparse(1:ncol,label,1,ncol,numel(u),ncol),1));
        H(k) = -dot(p,log2(p+eps));
    end
end

end