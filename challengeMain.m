%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     CHALLENGE MAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example generating a challenge submission. Loads the provided fluorescence file
%%% and performs the reconstruction (based on sample algorithms). As output it generates
%%% a scoring matrix whose ROC is computed against the true network (if provided).
%%%
%%% See also: challengeFastBaseline.m
%%% See also: challengeTrain.m
%%% See also: challengeVisualization.m
%%%

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Authors: Javier Orlandi, Mehreen Saeed, Isabelle Guyon
% Date: Jan 2014
% Last modified: 01/30/2014
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

clear all;
close all;
%% Begin -- User defined options --
challengeFolder = [pwd filesep 'krumby/'];
dataDirectory = [challengeFolder 'data'];
submissionDirectory = [challengeFolder 'results']; % Where ready-to-submit result files are found
modelDirectory = [challengeFolder 'models'];       % Where trained predictive models end up
networkIdNames = {'mocktest_small'};%'mockvalid', 'mocktest'};        % IMPORTANT: these are the base names of your data files in the data directory
%networkIdNames = {'iNet1_Size100_CC01inh', 'iNet1_Size100_CC02inh', 'iNet1_Size100_CC03inh', 'iNet1_Size100_CC04inh', 'iNet1_Size100_CC05inh', 'iNet1_Size100_CC06inh'};
scoringMethods = {@computeRW};%@randomScore, @computeGTE}; % Other methods include: @computeIGGini, @computeIGEntropy, @computeCrossCorrelation, @computeGranger, @trainedPredictor};     
                                 % Use:
                                 % 1) @randomScore to rapidly generate random results.
                                 % 2) @computeGTE to generate the baseline result with the GTE method to detect causality in time series. 
						 % @computeGranger to generate the Granger causality statistic (similar the GTE but with linear model).
                                 % 3) The following methods to treat the time series observations as iid points:
                                 % @computeIGGini calculates information gain using the Gini measure.
                                 % @computeIGEntropy calculates information gain using Entropy.
                                 % @computeCrossCorrelation computes the cross correlation.
                                 % 4) A trained model (see challengeTrain):
                                 % @trainedPredictor
modelName = 'sample_model';% Name of the model used by the trainedPredictor scoring method
concatenateScores = 1;     % 1 if all scores from various networks are appended to the same file (for each method)
% End -- User defined options --

%% Initializations
addpath(genpath(challengeFolder));
if ~exist(submissionDirectory,'dir') mkdir(submissionDirectory); end
extension='.txt';
the_date = datestr(now, 'yyyymmddTHHMMSS');
logfile = [submissionDirectory filesep 'logfile.txt'];
flog=fopen(logfile, 'a');
fprintf('==========================================================\n');
fprintf('\n ChaLearn connectomics challenge, sample code version %s\n', this_version);
fprintf(' Date and time started: %s\n', the_date);
fprintf(' Saving AUC results in %s\n\n', logfile);
fprintf('==========================================================\n\n');

metNum=length(scoringMethods);
netNum=length(networkIdNames);
scores=cell(netNum,metNum);
    
%% Loop over all methods you want to try
for j=1:metNum
    
    scoringMethod = scoringMethods{j};
    scoreFile = [submissionDirectory filesep func2str(scoringMethod) '_' sprintf('%s_', networkIdNames{:}) the_date];
    
    %% Loop over all networks you want to process   
    for i=1:netNum

        networkId = networkIdNames{i};
        fprintf('*** %s on %s ***\n\n', func2str(scoringMethod), networkId);

        %% Load the Fluorescence signal
        fluorescenceFile = [dataDirectory filesep 'fluorescence_' networkId extension];
        positionFile = [dataDirectory filesep 'networkPositions_' networkId extension];
        F = load_data(fluorescenceFile);
        P = load_data(positionFile);

        %% Compute the scores for all pairs of neurons i, j
        tic
        fprintf('Computing scores with %s\n', func2str(scoringMethod));
        if strcmp(func2str(scoringMethod), 'trainedPredictor')
            arg =  [modelDirectory filesep modelName '.mat'];
        else
            arg = false;
        end
        if strcmp(func2str(scoringMethod),'computeHash')
            scores{i,j} = scoringMethod(F, arg, P);
        else
            scores{i,j} = scoringMethod(F, arg);
        end
        % Note: these scoring methods do not make use of available neuron
        % positions (in their 2-D layout simulating neuron cultures).
        execution_time = toc

        %% Write the scores to the submission directory, ready to be submitted
        resuFile = [submissionDirectory filesep func2str(scoringMethod) '_' networkId '_' the_date];
        if ~concatenateScores, scoreFile = resuFile; end
        fprintf('Writing %s\n', scoreFile);
        writeNetworkScoresInCSV(scoreFile, scores{i,j}, networkId);

        %% If we have the network architecture... compute/plot the ROC curve:
        % (network architecture provided only for training data to the participants)
        networkFile = [dataDirectory filesep 'network_' networkId extension];
        fprintf('Computing ROC with using network %s\n', networkFile);
        if exist(networkFile,'file')
            network = readNetworkScores(networkFile);
            h = figure('Name', [func2str(scoringMethod) ' ' networkId ' ROC curve']);
            [AUC, FPR, TPR, TPRatMark, raw] = computeROC(network, scores{i,j}, 'plot', h);
            saveas(h, resuFile, 'png');
            fprintf('\n==> AUC = %5.4f\n\n', AUC);
            fprintf(flog, '\n%s\t%s\t%s\t%5.4f\t%5.0f', the_date, func2str(scoringMethod), networkId, AUC, execution_time);
        end
    end %loop over networks
end %loop over methods

fclose all;
MSG = 'Challenge solved.';
disp([datestr(now, 'HH:MM:SS'), ' ', MSG]);