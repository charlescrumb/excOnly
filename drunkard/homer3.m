


%% Random walk visualizations
% 20140312 - using kaggle contest data
load('spikeSamp.mat');
network = csvread('/Users/idaniboy/Documents/MATLAB/data/small/network_iNet1_Size100_CC01inh.txt');
inhibs = network(:,3)==-1;
network(inhibs,:)=[]; %remove inhibitory connections
nS = sparse(network(:,1),network(:,2),1,100,100);


%% Code for creating spike sample
%# spikeSamp.mat : spikeSamp1-3
%# in each sample, col1 = node, col2 = a possible (validated) outgoing edge
% nRows = size(network,1); % number of rows
% nSample = 90; % number of samples
% rndIDX = randperm(nRows); 
% spikeSamp = network(rndIDX(1:nSample), :); 
% save('spikeSamp.mat','spikeSamp');

% draw network edges
figure
spy(nS)
title('Existent Network Edges')

% Random walk
% outputs X which is list of nodes visited

% reset parameters
X = [];
j = 1;
stepProb = 0.7; % for now this is = to chance of self connection

% randomize nodes for walk
rw1 = randperm(size(spikeSamp1,1));
rw2 = randperm(size(spikeSamp2,1));
rw3 = randperm(size(spikeSamp3,1));

for i=1:10
    
    %record node
    switch j
        case 1
            X=[X;spikeSamp1(rw1(i),1)];
        case 2
            X=[X;spikeSamp2(rw2(i),1)];
            stepProb = 0.6;
        case 3
            X=[X;spikeSamp3(rw3(i),1)];
            stepProb = 0.5;
    end
    
    if numel(rw1)==1
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


%% Random walk visualizations


% get coords on circle for visualization
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

% Figure for connectivity matrix
spos = get(0,'screensize');
figure('position',[spos(3)/10 spos(4)/10 spos(3)*.8 spos(4)*.8],'color','w');
gplot(nS,coordsM); hold on
view([60 -30 ]);

% Draw number of visits at each node
Visited = unique(X);
xdata = zeros(size(Visited)); ydata = xdata; zdata = xdata;
for ii = 1:numel(Visited)
    zdata(ii) = sum(X == Visited(ii));
    xdata(ii) = coordsM(Visited(ii),1);
    ydata(ii) = coordsM(Visited(ii),2);
    line(xdata(ii)*[1;1],ydata(ii)*[1;1],[0 zdata(ii)],'color','b',...
        'linewidth',2);
end
set(gca,'xtick',[],'ytick',[]);
hold on

% Draw # of edge traversals
weightsM = zeros(100,100);

L = [];
for ii = 1:numel(X)
    xdata = coordsM(X(ii),1);
    ydata = coordsM(X(ii),2);
    if isempty(L)
        L = line(xdata,ydata,'color','r','marker','o','markersize',12,...
            'markeredgecolor','k','markerfacecolor','r');
    else
        set(L,'xdata',xdata,'ydata',ydata);
        weightsM(X(ii),X(ii-1))=weightsM(X(ii),X(ii-1))+1; %edge count
    end
%     pause(0.0001);
%     drawnow
end
hold on

% 
figure
imagesc(weightsM)

% Using wgplot
figure;
[he,hv]=wgPlot(weightsM,coordsM);
