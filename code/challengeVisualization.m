%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 CHALLENGE VISUALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example of the various functions available to visualize the Challenge
%%% data. The figures, the network in Gelphi format suitable for 
%%% visualization and the movie in AVI format are exported to the
%%% graphicDirectory.
%%%
%%% See also: challengeMain.m
%%% See also: challengeFastBaseline.m
%%% See also: challengeTrain.m

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
challengeFolder = [pwd filesep];
dataDirectory = [challengeFolder 'data'];
graphicDirectory = [challengeFolder 'graphic'];
networkIdNames = {'mocktrain', 'mockvalid','mocktest'}; % IMPORTANT: these are the base names of your data files in the data directory
% End -- User defined options --

%% Keyboard input preferences
fprintf('List of networks:\n');
for k=1:length(networkIdNames), fprintf('%d: %s\n', k, networkIdNames{k}); end
selection = input('Choose your network: ');  
networkID = networkIdNames{selection};

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initializations and data loading %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Loading data, please be patient...\n');

addpath(genpath(challengeFolder));
extension='.txt';
%create a movie folder if it does not exist
if ~exist(graphicDirectory,'dir') mkdir(graphicDirectory); end

% Load the network and fluorescence values:
[F, D, G, network, totalSamples, totalNeurons, samples] = GetNetworkNameInput(networkIdNames,[],[],[],[],dataDirectory,[],selection);
networkFile = [dataDirectory filesep 'network_' networkID extension];

% Load the neuron positions, create a Network structure:
Network.RS = network;
networkPositionsFile = [dataDirectory filesep 'networkPositions_' networkID extension];
fprintf('Loading %s\n', networkPositionsFile);
temp = load(networkPositionsFile);
Network.X = temp(:,1);
Network.Y = temp(:,2);

% Define the initial pair of neurons
neuron_i = 10;
neuron_j = 14;

% Define the maximum number of samples
maxSamples = min(1000,totalSamples); % Show only part of the time series

% Define the visualization offset (so time series on the same axis do not
% overlap)
offset=0.2;

fprintf('... data loaded.\n');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  DATA VISUALIZATION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

option=true;
while option
    fprintf('Options:\n1: View fluorescence data.\n2: View network data.\n3: View neuron movie.\n4: Browse through the time series.\n0: Quit.\n');
    option = input('Choose your option: ');


    switch option
    
%% VISUALIZATION OF FLUORESCENCE DATA
    case 1 

        fprintf(2, 'Plotting fluorescence data.\nSaved to directory %s\n', graphicDirectory);

        % Plot the average of the whole time series
        figure('Name', 'visualizeFluorescenceTimeSeries: Average', 'Position',[440 378 560 420]);
        h = visualizeFluorescenceTimeSeries(F, 'samples', 1:size(F,1));
        saveas(gcf, [graphicDirectory filesep 'averageSeries_' networkID], 'jpg');

        % Plot all the individual traces separated
        figure('Name', 'visualizeFluorescenceTimeSeries: All series slightly offset', 'Position',[450 358 560 420]);
        h = visualizeFluorescenceTimeSeries(F, 'type', 'single','offset',offset/10, 'neuronList', 1:size(F,2));
        saveas(gcf, [graphicDirectory filesep 'allSeries_' networkID], 'jpg');

        % Visualize a single neuronal pair; compare the original and discretized signals
        % Note: these functions expect a structure Network, with fields RS, X, and Y.
        figure('Name', 'visualizePair: Compare original and discretized', 'Position',[460 338 560 420]);
        subplot(2,1,1); 
        h = visualizePair(Network, F, neuron_i, neuron_j, 'offset', offset,'samples',1:maxSamples);
        subplot(2,1,2);
        h = visualizePair(Network, D, neuron_i, neuron_j, 'offset', offset,'samples',1:maxSamples);
        saveas(gcf, [graphicDirectory filesep 'pair_' num2str(neuron_i) '_' num2str(neuron_j) '_' networkID], 'jpg');

    
%% VISUALIZATION OF NETWORK DATA
    case 2

        % Plot the network connections as a matrix
        h=figure('Name', 'Network connections', 'Position',[470 318 560 420]);
        gmat_display(full(network),[],[],h);
        saveas(gcf, [graphicDirectory filesep 'networkConnections_' networkID], 'jpg');

        % Save the network as a GEXF file for visualization with Gephi
        fprintf(2, 'Saved the network as a GEXF file for visualization with Gephi\n in directory: %s\n', graphicDirectory);
        if exist(networkFile,'file')
            gephiFile = [graphicDirectory filesep 'network_' networkID  '.gexf'];
            networkToGEXF(Network, gephiFile);
        end
        % Now you can open the GEXF file with Gephi to visualize the network.
        % Download and install Gephi from: https://gephi.org/users/download/.

%% MOVIE
    case 3
        % Fluorescence movie, 1 in every 10 samples and no saving
        if strfind(computer, 'PCWIN') % Works only under windows
            cmap = hot(256);
            visualizeFluorescenceMovie(Network, F, 'colorMode', 'absolute','samples',1:10:maxSamples, 'fileName',[],'cmap', cmap);
        end

        % Save the movie, but only the fluorescence data
        if strfind(computer, 'PCWIN') % Works only under windows
            fprintf(2,'Saving movie to directory %s\n', graphicDirectory);
            visualizeFluorescenceMovie(Network, F, 'colorMode', 'absolute','movieMode','fast','samples',1:maxSamples, 'fileName',[graphicDirectory filesep 'movie_' networkID]);
        else
            fprintf(2,'Cannot create movie: not running Windows.\n');
        end
        %you can use implay to play the movie or any AVI player

%% DATA BROWSER
    case 4

        h=figure('Name', ['Connectomics Neuron Pair Browser for network: ' networkID]);
        set(h,'OuterPosition',[150 50 1100 750]);

        samples = 1:maxSamples; % Show only part of the time series

        % plot the data
        plotData(F,D,neuron_i,neuron_j,samples,offset);
        displayFeatures(F,D,neuron_i,neuron_j,G,network,h);

        % Network
        NetworkButton=uicontrol('Parent', h, 'Position', [5 400 80 40], 'FontSize', 12, 'BackgroundColor', [1 0.9 0.1], 'String', 'Network', ...
                    'Callback','[F D G network  totalSamples totalNeurons samples] = GetNetworkNameInput(networkIdNames,F,D,G,network,dataDirectory,h);plotData(F,D,neuron_i,neuron_j,samples,offset);displayFeatures(F,D,neuron_i,neuron_j,G,network,h);');
        % Prev
        PrevButton=uicontrol('Parent', h, 'Position', [5 360 80 40], 'FontSize', 12, 'BackgroundColor', [1 0.9 0.1], 'String', '< Prev', ...
                    'Callback','if neuron_i~=1 neuron_i=neuron_i-1;elseif neuron_j~=1 neuron_j=neuron_j-1;end;plotData(F,D,neuron_i,neuron_j,samples,offset);displayFeatures(F,D,neuron_i,neuron_j,G,network,h);');
        % Next
        NextButton=uicontrol('Parent', h, 'Position', [5 320 80 40], 'FontSize', 12, 'BackgroundColor', [1 0.9 0.1], 'String', 'Next >',...
                    'Callback', 'if (neuron_i~=totalNeurons) neuron_i=neuron_i+1; elseif neuron_j~=totalNeurons neuron_j = neuron_j+1;end;plotData(F,D,neuron_i,neuron_j,samples,offset);displayFeatures(F,D,neuron_i,neuron_j,G,network,h);');
        % Jump to another pair
        JumpButton=uicontrol('Parent', h, 'Position', [5 280 80 40], 'FontSize', 12, 'BackgroundColor', [1 0.7 0.1], 'String', 'Jump to', ...
                    'Callback', '[neuron_i neuron_j] = getInput(totalNeurons);plotData(F,D,neuron_i,neuron_j,samples,offset);displayFeatures(F,D,neuron_i,neuron_j,G,network,h);');
        % get different samples
        SamplesButton=uicontrol('Parent', h, 'Position', [5 240 80 40], 'FontSize', 12, 'BackgroundColor', [1 0.7 0.1], 'String', 'Samples', ...
                    'Callback', 'samples=getSampleInput(totalSamples);plotData(F,D,neuron_i,neuron_j,samples,offset);displayFeatures(F,D,neuron_i,neuron_j,G,network,h);');
        % get the offset
        OffsetButton=uicontrol('Parent', h, 'Position', [5 200 80 40], 'FontSize', 12, 'BackgroundColor', [1 0.7 0.1], 'String', 'Offset', ...
                    'Callback', 'offset=getOffset(offset);plotData(F,D,neuron_i,neuron_j,samples,offset);displayFeatures(F,D,neuron_i,neuron_j,G,network,h);');
        break
        
        case 0
            %do nothing here...as the option is for quit
        otherwise
            fprintf(2, 'No such option!\n');
    end  % Switch option
end  % While not_quit_me

