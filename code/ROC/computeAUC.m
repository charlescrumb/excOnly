function auc = computeAUC(Target, Scores)
% auc = computeAUC((Target, Scores)
%
% Calculates the area under the ROC for a given set
% of predictions scores and target labels.  Currently limited to two classes.
%
% Scores: n*1 matrix of posterior probabilities for class 1
% Target: n*1 matrix of categories {0,1}
% auc: Area under the curve

% Algorithm found in
% A  Simple Generalisation of the Area Under the ROC
% Curve for Multiple Class Classification Problems
% David Hand and Robert Till
% http://www.springerlink.com/content/nn141j42838n7u21/fulltext.pdf
% Should be consistent with the values shown on the Kaggle leaderboard.

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Author: Ben Hamner
% Date: Jan 2013
% Last modified: NA
% Contact: ben@benhamner.com
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

category=Target(:);
posterior=Scores(:);

if exist('tiedrank')
    r = tiedrank(posterior);
else
    r = tiedrank_metrics(posterior);
end
auc = (sum(r(category==1)) - sum(category==1)*(sum(category==1)+1)/2) / ...
    ( sum(category<1)*sum(category==1));

function r = tiedrank_metrics(x)

[~,I] = sort(x);

r = 0*x;
cur_val = x(I(1));
last_pos = 1;

for i=1:length(I)
    if cur_val ~= x(I(i))
        r(I(last_pos:i-1)) = (last_pos+i-1)/2;
        last_pos = i;
        cur_val = x(I(i));
    end
    if i==length(I)
        r(I(last_pos:i)) = (last_pos+i)/2;
    end
end