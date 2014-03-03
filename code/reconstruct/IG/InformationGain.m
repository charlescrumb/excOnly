function IG = InformationGain(A,B,method )
%IG = InformationGain(A,B,method )
%Compute the information gain using the following three possible methods
%max classification error, gini or entropy.
%Here we assume that A is the attribute (input) and B is the label
%(output), so the value IG returned will be a measure of A->B (larger
%values for A->B, smaller values for B->A.
% The function returns
% IG(B, A) = I(B) - I(B|A)
% where I is the "method".
%Set "method" as 'entropy','gini' or 'maxError', default is 'entropy'.

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Author: Mehreen Saeed
% Date: Jan 2014
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

if nargin<3
    method = 'entropy';
end

% Compute the "impurity" in B I(B)
dataSubset = B;
str = [method '(dataSubset);'];
measure = eval(str);

% Compute 
IG = measure;
vals = unique(A);
for i=1:1:length(vals)
    ind = find(A==vals(i));
    prob = length(ind)/length(A);
    dataSubset = B(ind);
    IG = IG - prob*eval(str);
end


%-------------------------------------------------------------------
%compute entropy
function e = entropy(X)

total = length(X);
vals = unique(X);
e=0;

%is it a constant
if length(vals)==1
    return; 
end

for i=1:1:length(vals)
    ind = find(X==vals(i));
    prob = length(ind)/total;
    e = e - (prob*log2(prob));
end


%-------------------------------------------------------------------
%compute gini
function e = gini(X)

total = length(X);
vals = unique(X);
e=1;


for i=1:1:length(vals)
    ind = find(X==vals(i));
    prob = length(ind)/total;
    e = e - (prob*prob);
end



%-------------------------------------------------------------------
%compute max classification error
function e = maxError(X)

total = length(X);
vals = unique(X);

prob = [];

for i=1:1:length(vals)
    ind = find(X==vals(i));
    prob = [prob length(ind)/total];
end

e = 1 - max(prob);