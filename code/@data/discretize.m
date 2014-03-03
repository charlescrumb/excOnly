function discretize(this, varargin)
% DISCRETIZE discretizes the fluorescence signal so it
% can be used to compute the joint PDF. If conditioning is applied, the
% entries above the conditioning level are returned in the G vector.
% Fills in data members D (same dimension as F) and G (Vector defining the global conditioning level of the signal at
% that given time (for now 1 and 2 for below and above the level).
%
% USAGE:
%    discretize(this, varargin)
%
% optional arguments ('key' followed by its value): 
%    'bins' - Number of bins to use in the discretization (before
%    conditioning). If the entry is a vector, it will define the bin edges
%    (default 3).
%
%    'relativeBins' - (true/false). If true the bins are defined based on
%    the min and max values of each sample. If false, they are defined
%    based on the absolute limits (default false).
%
%    'conditioningLevel' - Value used for conditioning. If the
%    value is 0 the level is guessed (for now, the peak of the histogram
%    plus 0.05). Set it to inf to avoid conditioning (default 0).
%
%   'highPassFilter' - (true/false). Apply a high pass filter to the
%   fluorescence signal, i.e., work with the derivative (default true).
%
%    'debug' - true/false. Show additional partial information (default
%    false).
%
% EXAMPLE:
%    discretize(this, 'bins', 3, 'debug', true);
%
%    (Stetter 2013) Stetter, O., Battaglia, D., Soriano, J. & Geisel, T. 
%    Model-free reconstruction of excitatory neuronal connectivity from 
%    calcium imaging signals. PLoS Comput Biol 8, e1002653 (2012).
%
% SEE ALSO: [D, G] = discretizeFluorescenceSignal(F, varargin)

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Author: Javier Orlandi
% Date: Dec 2013
% Last modified: 01/29/2014 by Isabelle Guyon
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

if ~isempty(this.D) && ~this.rediscretize, return; end % IG: Do not recompute if already done!

%%% Assign default values
params.bins = [-10,0.12,10]; % default for GTE = [-10,0.12,10], default original: 3
params.relativeBins = false;
params.conditioningLevel = 0.25; % Old default = 0
params.debug = false;
params.highPassFilter = true;
params = parse_pv_pairs(params,varargin); 

epsilon = 1e-3; % To avoid weird effects at bin edges
if length(params.bins)==1
    numBins=params.bins;
else
    numBins=length(params.bins)-1;
end

%%% Get the conditioning level
avgF = mean(this.F,2); % Average time series
% IG: with 'conditioningLevel', 0.25 we don't go throught this
if(params.conditioningLevel == 0)
    [hits, pos] = hist(avgF, 100);
    [~, idx] = max(hits);
    CL = pos(idx)+0.05;
    fprintf('Best guess for conditioning found at: %.2f\n', CL);
else
    CL = params.conditioningLevel;
end

%%% Apply the conditioning
this.G = (avgF >= CL)+1;  % two levels: 1=av. activity below threshold; 2=av. activity above threshold.

%%% Show the result of conditioning
if(params.debug)
    figure;
    h = plotFluorescenceHistogram(F,'bins',100);
    hold on;
    yl = ylim;
    hCL = plot([1, 1]*CL, yl,'k');
    legend([h, hCL], 'Average F histogram', 'Conditioning Level');
    xlabel('Fluorescence');
    ylabel('hits');
end

%%% Apply the high pass filter
[numSamp, numNeuron]=size(this.F);
if(params.highPassFilter)
%    F = diff(F); %IG: to save memory, we do not do that
    %G = G(1:end-1);
    this.G = this.G(2:end);
    numSamp=numSamp-1;
end

%%% Discretize the signal
%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Starting to discretize [numBins=%d, relativeBins=%d, conditioningLevel=%g, highPassFilter=%d]...\n', ...
    numBins, params.relativeBins, CL, params.highPassFilter);
tic
this.D = NaN(numSamp, numNeuron);
if(length(params.bins) > 1)
    params.relativeBins = false; % Just in case
end

if(params.relativeBins)
    % IG: by default we do not go through this
    for i = 1:size(this.F, 2) % loop over neurons
        phi=this.F(:, i);
        if(params.highPassFilter)
            phi = diff(phi); %IG: we do the high pass filter 1 neuron at a time to save memory
        end
        binEdges = linspace( min(phi)-epsilon, max(phi)+epsilon, params.bins+1);
        for j = 1:(length(binEdges)-1)
            hits = phi >= binEdges(j) & phi < binEdges(j+1);
            D(hits, i) = j;
        end
    end        
else
    binEdges = params.bins; % IG: binEdges in the case of length(params.bins) == 1 is now inside the loop
    for i = 1:size(this.F, 2) % IG: Loop over neurons to save memory
        phi=this.F(:, i);
        if(params.highPassFilter)
            phi = diff(phi); 
        end
        if(length(params.bins) == 1)
            % IG: with option 'bins', [-10,0.12,10] we do not go through this
            binEdges = linspace(min(phi)-epsilon, max(phi)+epsilon, params.bins+1);
        end
        for j = 1:(length(binEdges)-1)
            hits = phi >= binEdges(j) & phi < binEdges(j+1);
            this.D(hits, i) = j;
        end
    end
    if(params.debug)
        fprintf('Global bin edges set at: (');
        for j = 1:(length(binEdges)-1)
            fprintf('%.2f,', binEdges(j));
        end
        fprintf('%.2f).\n', binEdges(end));
    end
end
toc
fprintf('...done\n');

end