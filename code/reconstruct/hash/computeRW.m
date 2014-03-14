function [E] = computeHash(F, debug_, P)
%% Uses econometrics toolbox random walk function

clear HashTable
global HashTable
%workaround initialization of map structure
HashTable=containers.Map(uint32(1), [1 1]);
remove(HashTable,1);

%% parameters
%search window in frames
searchWindow=2;

%% Discretize the fluorescence signal
%this is what was used in starter kit. i dont understand it much
[D, G] = discretizeFluorescenceSignal(F, 'debug', false, 'conditioningLevel', 0.25, 'bins', [-10,0.12,10]);

%oopsi params. it went a lot slower than starter kit method, abandoned but
%left in case i was doing it wrong
%V=struct('fast_iter_max',4,'smc_iter_max',0,'dt',1/50,'preprocess',1);

%initialize spike index matrix. discretize takes 1 frame away
%booleanSpikes=logical(sparse([],[],[],size(F,1)-1,size(F,2)));



%much slower oopsi attempt, maybe i did it wrong  
%send each channel to oopsi
% for i=1:size(F,2)
%     fast=run_oopsi(F(:,i),V);
%     n_fast=fast.n/max(fast.n);
%     spikeThresh=mean(n_fast)+4.5*std(n_fast);
%     %spts=find(n_fast>spikeThresh);
%     tic
%     booleanSpikes(:,i)=n_fast>spikeThresh;
%     toc
% end

%2's in D were around spike times, not sure how it worked though. it had 1
%less row than F, maybe it looks at past times and frame 1 doesn't have them?
booleanSpikes=logical(D==2);

nSamples=size(booleanSpikes,1);
nNeurons=size(F,2);
weightsM = zeros(nNeurons);

tic
E=zeros(nNeurons);
timeStacked=zeros(nNeurons, nNeurons, searchWindow+1);
for i=1:nSamples-searchWindow
    
    spikeSamp1 = find(booleanSpikes(i,:))';

    if spikeSamp1 ~= 0
        
        spikeSamp2 = find(booleanSpikes(i+1,:))';
        
        if spikeSamp2 ~= 0

            spikeSamp3 = find(booleanSpikes(i+2,:))';
            
            if spikeSamp3 ~=0

                % Random walk, requires at least 3 spikes in row
                % outputs X which is list of nodes visited

                % reset parameters
                X = [];
                j = 1;
                stepProb = 0.999;

                for i=1:500

                    %record node
                    switch j
                        case 1
                            rSamp = datasample(spikeSamp1,1,'Replace',false);
                            sSize = numel(spikeSamp1);
                        case 2
                            stepProb = 0.95;
                            rSamp = datasample(spikeSamp2,1,'Replace',false);
                            sSize = numel(spikeSamp2);
                        case 3
                            rSamp = datasample(spikeSamp2,1,'Replace',false);
                            stepProb = 0.8;
                            sSize = numel(spikeSamp3);
                    end
                    
                    X=[X;rSamp];
                    
                    if sSize==1
                        % adjust stepProb so 0% chance of self connection 
                        % [ ] can look up real prob later
                        stepProb = 0;
                    end

                    %flip coin to determine if move to next frame
                    if rand>stepProb

                        %move to next frame
                        j = j+1;

                        if j==4
                            break    
                        end

                    end

                end

                for ii = 2:numel(X)
                        weightsM(X(ii),X(ii-1))=weightsM(X(ii),X(ii-1))+1; %edge count
                end
            end
        end
    end
end


toc
disp('sie')
tic

        figure
        imagesc(weightsM)

        % Draw points on circle
A      = 1;  %starting point
p      = 101; %number of points
angles = linspace(0,2*pi,p); 
f      = @(theta) [cos(theta),sin(theta)]; 
r      = arrayfun(@(i)f(angles(i)),1:p,'UniformOutput',false); 
coordsM     = reshape(cell2mat(r),[2 p]); 
% figure
% line(coordsM(1,:),coordsM(2,:)); 
% hold on 
% plot(coordsM(1,:),coordsM(2,:),'ro');
coordsM(:,end)=[]; % remove last coordinate (which ties together circle)
coordsM = coordsM';


        % Using wgplot
        figure;
        [he,hv]=wgPlot(weightsM,coordsM);


E = weightsM;
E=E./max(max(max(E)));

%setting neuron1=neuron2 counts to zero improved the score a lot.
E(logical(eye(size(E)))) = 0;

end

