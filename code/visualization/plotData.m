function plotData(F,D,neuron_i,neuron_j,samples,offset)
%%% plotData(F,D,neuron_i,neuron_j,samples)
%%% function to plot fluorescence signals F and its discretized version D
%%% for the (neuron_i,neuron_j) pair using samples specified as parameter

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Authors: Mehreen Saeed
% Date: Jan 2014
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

subplot('Position',[0.13 .6 0.4 .35])
plot(samples,F(samples,neuron_i),'r',samples,F(samples,neuron_j)+offset);
legend(['neuron ' num2str(neuron_i)],['neuron ' num2str(neuron_j)]);
xlabel('sample nunmber');
ylabel('Fluorscence signal');

%%plot the discretized signal
minival = min(min(D(samples,[neuron_i neuron_j])));
maxival = max(max(D(samples,[neuron_i neuron_j])))+offset;

subplot('Position',[0.13 .15 0.4 .35])
plot(samples,D(samples,neuron_i),'r',samples,D(samples,neuron_j)+offset);
legend(['neuron ' num2str(neuron_i)],['neuron ' num2str(neuron_j)]);
xlabel('sample nunmber');
ylabel('Discrete Fluorscence signal');
axis([min(samples) max(samples) minival-.5 maxival+.5]);

