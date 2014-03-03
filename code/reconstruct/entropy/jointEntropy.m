function H = jointEntropy(X, Y, dim)
%H = jointEntropy(X)
%H = jointEntropy(X, Y)
%H = jointEntropy(X, [], dim)
% Compute joint entropy H(X,Y) of two discrete variables X and Y.
% This works only for integer discrete levels. It avoids calling "unique".
% If X is a matrix (p,n), compute the (n,n) matrix joint entropies
% of all pairs of columns. 
% if Y={}, and dim is given, dim =1 => use columns; dim =2 => use lines 
% and return a (p,p) matrix.
% See also: entropy.
% Note: jointEntropy(X,X)=entropy(X).

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

debug_=0;
    
[nlin, ncol]=size(X);
if nlin>1 && ncol>1, % matrix case
    if nargin>1,
        if ~isempty(Y),
            help jointEntropy
            return
        end
    end
else
    if nargin<2, 
        Y=[]; 
    else
        assert(all(size(X)==size(Y)));
    end
end

if nargin<3, 
    if ncol==1, 
        dim=1; 
        X=[X,Y];
    elseif nlin==1,
        dim=2;
        X=[X;Y];
    else
        dim=1;
    end
end

% Here we assume that dim=1 is the preferred dimension and we won't have to
% transpose. Otherwise, duplicate the code for the transposed version to
% get a memory saving.
if dim==2, X=X'; end

% This assumes discrete integer values
[p,n]=size(X);
m = min(min(X));
if m~=1, X = X-m+1; end % one base the data

M = max(max(X));
idx = 1:nlin;
H=zeros(ncol,ncol);
tot=ncol*(ncol+1)/2;
k=0;
for i=1:ncol
    for j=i:ncol
        p = nonzeros(sparse(idx,X(:,i),1,nlin,M,nlin)'*...
                     sparse(idx,X(:,j),1,nlin,M,nlin)/nlin); %joint distribution of x and y
        H(i,j) = -dot(p,log2(p+eps));
        H(j,i) = H(i,j);
        k=k+1;
        if mod(k/tot*100,10)==0 && debug_
            fprintf('\n percentage done : %d%%',i/ncol*100);
        end
    end
end
fprintf('\n');
end
