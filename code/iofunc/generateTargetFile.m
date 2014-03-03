%script for generating target files in CSV
%-------------------INPUTS
targetFile = 'E:\connectomics\GTE-Challenge-master\MATLAB\data\network_train1'
networkFile = 'E:\connectomics\GTE-Challenge-master\MATLAB\challenge\network_iNet1_Size100_CC03.txt'
networkID = 'targets_train1';
%--------------------------------------------------------------
networkData = load(networkFile);
N = max(max(networkData(:,1:2)));
targets = sparse(networkData(:,1), networkData(:,2), networkData(:,3), N, N);

writeNetworkScoresInCSV(targetFile,targets,networkID);