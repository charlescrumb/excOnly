function [E] = computeHash(F, debug_, P)
%% Uses econometrics toolbox random walk function
%% DOESN'T WORK!

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
    [r, c]=find(booleanSpikes(i:i+searchWindow,:));

    if r ~= 0
        % connM = sparse(r,c,1,3,100);

        % Set up added noise matrix parameters.
        randVec = rand(300,1)<0.5;
        m = reshape(randVec,3,100);
        
        % create connection matrix
%         connM = sparse(m);

        % bug fix; add random noise to matrix
        % fullN = full(connN);
        m = m | booleanSpikes(i:i+searchWindow,:);
        connM = sparse(m);

        
        
        walks = numel(c);

        % Random walk
        Ind = 1;
        F = @(t,X) 0;
        G = @(t,X) 1;
        S = sde(F,G,'startState',Ind);
        X = S.simByEuler(walks,'ntrials',1,'Z',@(t,X) RandomGraphMove(X,connM)-X);

    %     % Figure for connectivity matrix
    %     spos = get(0,'screensize');
    %     figure('position',[spos(3)/10 spos(4)/10 spos(3)*.8 spos(4)*.8],'color','w');
    %     gplot(connM,coordsM); hold on
    %     view([60 -30 ]);

        % Draw number of visits at each node
    %     Visited = unique(X);
    %     xdata = zeros(size(Visited)); ydata = xdata; zdata = xdata;
    %     for ii = 1:numel(Visited)
    %         zdata(ii) = sum(X == Visited(ii));
    %         xdata(ii) = coordsM(Visited(ii),1);
    %         ydata(ii) = coordsM(Visited(ii),2);
    %         line(xdata(ii)*[1;1],ydata(ii)*[1;1],[0 zdata(ii)],'color','b',...
    %             'linewidth',2);
    %     end
    %     set(gca,'xtick',[],'ytick',[]);
    %     hold on

        % Draw # of edge traversals

        L = [];
        for ii = 2:numel(X)

                weightsM(X(ii),X(ii-1))=weightsM(X(ii),X(ii-1))+1; %edge count

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

