function scores = computeIGEntropy(F, debug_)
% scores = computeIGGini(F, debug_)
% Baseline method to compute scores based on information Gain measure using
% Entropy

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Author: Mehreen Saeed
% Date: Jan 2014
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

if nargin<2,
    debug_ = false; 
end

%% Discretize the fluorescence signal
D = discretizeFluorescenceSignal(F,'bins',3);

%% compute the information gain measure
col = size(D, 2);
arr = [];
for i=1:1:col
    for j=i+1:1:col
        Values_ij = InformationGain(D(:,i),D(:,j),'entropy');
        Values_ji = InformationGain(D(:,j),D(:,i),'entropy');
        
        arr = [arr;i j Values_ij];
        arr = [arr;j i Values_ji];
    end
    
    if debug_
        if mod(i/col*100,10)==0
            fprintf('\n percentage done : %d',i/col*100);
        end
    end
end
if debug_
	fprintf('\n');
end

scores =  sparse(arr(:,1), arr(:,2), arr(:,3), col,col);

end
