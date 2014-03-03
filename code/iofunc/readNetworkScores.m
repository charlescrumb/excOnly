function scores = readNetworkScores(filename,N)
%scores = readNetworkScores(filename,N)
% Read the network architecture from the original format as a spase nxn matrix, 
% whose (i,j)th entry indicates neuron i is connected to neuron j 
% with a certain strength.
% N is the dimension of the scores matrix....if it is not specified it will
% be inferred.

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Authors: Isabelle Guyon
% Date: Jan 2014
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

networkData = load(filename);
if nargin<2
    N = max(max(networkData(:,1:2)));
end

% Keep only 0/1 weights, ignore blocked connections
networkData(networkData(:,3)>0, 3) = 1;
networkData(networkData(:,3)<0, 3) = 0;

scores = sparse(networkData(:,1), networkData(:,2), networkData(:,3), N, N);