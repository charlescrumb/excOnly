function [scoresPval scoresFstat] = computeGranger(F, debug_,timeLag)
% scores = computeGranger(F)
% Baseline method to compute scores based on Granger causality index.

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
    timeLag = 1;
end


%% compute Granger causality index for all pairs
[row col] = size(F);

scoresPval = zeros(col,col);
scoresFstat = zeros(col,col);
for i=1:1:col
    for j=i+1:1:col
        [pij fstatij] = calculateIndex(F(:,i),F(:,j),timeLag);        
        [pji fstatji] = calculateIndex(F(:,j),F(:,i),timeLag);        

        scoresPval(i,j) = pij;
        scoresPval(j,i) = pji;
        
        scoresFstat(i,j) = fstatij;
        scoresFstat(j,i) = fstatji;
        
    end
    
    if debug_
        if mod(i/col*100,10)==0
            fprintf('Percentage done : %d\n',i/col*100);
        end
    end
end

function [pVal Fstatistic] = calculateIndex(y,x,p)
%%% function to compute granger causality index when given two time series
%%% x and y.  x and y should be the same length.  p are the total time lag values
%%% of x and y used for predicting x.
%%% The F statistic is computed from the residual sum of square values
%%% Web reference: http://support.sas.com/rnd/app/examples/ets/granger/
%%% following their notations
%%% returns the p value of the Ftest that y Granger causes x

%% data size
timeSteps = length(x);
%% step 1: compute sum of square of residuals when regressing x using past p
%values of x
currentX = x(1+p:timeSteps,1);
pastX = [ones(timeSteps-p,1)];
for i=1:p
    pastX = [pastX x(i:timeSteps-i,1)];
end
%% find the residual sum of squares from regression
[b1 bint1 residuals1] = regress(currentX,pastX);
RSS0 = residuals1'*residuals1;

%% perform regression by including values of y
for i=1:p
    pastX = [pastX y(i:timeSteps-i,1)];
end
[b2 bint2 residuals2] = regress(currentX,pastX);
RSS1 = residuals2'*residuals2;

%% compute the fstatistic
Fstatistic = (RSS0-RSS1)/p;
Fstatistic = Fstatistic/RSS1*(timeSteps-2*p-1);
%% this is approximately and F distribution with DOF (p,T-2p-1)
pVal = 1-fcdf(Fstatistic,p,timeSteps-2*p-1);