function scores = computeCrossCorrelation(F, debug_, timeLag)
%%% scores = computeCrossCorrelation(F)
%%% Baseline method to compute scores based on cross correlation measure 
%%% Cross correlation is normalized
%%% default value for timeLag is zero
%%% if other timeLag value is specified then crossCorrelation values are
%%% averaged for connection strength i->j as average
%%% crossCorrelation(-timeLag:0) and for j->i are averaged for
%%% crossCorrelation(0:timeLag)

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Authors: Mehreen Saeed
% Date: Jan 2014
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

if nargin<2,
    debug_ = false; 
end

if nargin<3
    timeLag = 0;
end

%% Discretize the fluorescence signal
F = discretizeFluorescenceSignal(F,'bins',3);

%% compute cross correlation
[row col] = size(F);
arr = [];
for i=1:1:col
    for j=i+1:1:col
        crossCorr = myxcorr(F(:,i),F(:,j),timeLag);
        Values_ij = mean(crossCorr(1:timeLag+1));        
        Values_ji = mean(crossCorr(timeLag+1:2*timeLag+1));
        
        arr = [arr;i j Values_ij];
        arr = [arr;j i Values_ji];
    end
    
    if debug_
        if mod(i/col*100,10)==0
            fprintf('Percentage done : %d\n',i/col*100);
        end
    end
end

scores =  sparse(arr(:,1), arr(:,2), arr(:,3), col,col);



function crossCorr = myxcorr(x,y,timeLag)
%%% function to compute cross correlation between two time series
    
index = 1;
n = length(x);
for i=-timeLag:1:-1
    j = -i;
    crossCorr(index) = y(j+1:n)'*x(1:n-j);
    index = index+1;
end

for i=0:timeLag
    
    crossCorr(index) = x(i+1:n)'*y(1:n-i);
    index = index+1;
end
 
% normalize crossCorrelation
crossCorr = crossCorr/sqrt(x'*x)/sqrt(y'*y);