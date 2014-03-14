


%% Random walk visualizations
% 20140312 - using kaggle contest data
load('spikeSamp.mat');
network = csvread('/Users/idaniboy/Documents/MATLAB/data/small/network_iNet1_Size100_CC01inh.txt');
inhibs = network(:,3)==-1;
network(inhibs,:)=[]; %remove inhibitory connections
networkSp = sparse(network(:,1),network(:,2),1,100,100);


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
spy(networkSp)
title('Existent Network Edges')

% create self-connectivity matrix (edge-list)for spikeSamps1 (n=10)
eS1 = combnk(spikeSamp1(:,1),2); %list of edges

% create self-connectivity matrix (edge-list) for spikeSamps2 (n=100; however limit n
% to 12 b/c of combinatorial explosion)
nRows = size(spikeSamp2,1); % number of rows
nSample = 12; % number of samples
rndIDX = randperm(nRows); 
miniSpikeSamp2 = spikeSamp2(rndIDX(1:nSample), 1);
eS2 = combnk(miniSpikeSamp2(:,1),2); %list of edges

% create self-connectivity matrix (edge-list) for spikeSamps3 (n=20)
nRows = size(spikeSamp3,1); % number of rows
nSample = 10; % number of samples
rndIDX = randperm(nRows); 
miniSpikeSamp3 = spikeSamp3(rndIDX(1:nSample), 1);
eS3 = combnk(miniSpikeSamp3(:,1),2); %list of edges

% create connectivity matrix (edge-list) for spikeSamps1-->spikeSamps2 (n=10x100)
[p,q] = meshgrid(spikeSamp1(:,1), miniSpikeSamp2(:,1));
s1c2 = [p(:) q(:)];

% create connectivity matrix (edge-list)for spikeSamps2-->spikeSamps3 (n=100x20)
[p,q] = meshgrid(miniSpikeSamp2(:,1), miniSpikeSamp3(:,1));
s2c3 = [p(:) q(:)];

edges = [eS1;eS2;eS3;s1c2;s2c3];

% create connection matrix
edges = unique(edges,'rows');
connN = sparse(edges(:,1),edges(:,2),1,100,100);
%connM = sparse(edges(:,1),edges(:,2),1,100,100);
figure; spy(connN);

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



%% Random walk visualizations
% Set up parameters.
rows = 100;
columns = 100;
onesPerColumn = 2;
% Initialize matrix.
m = zeros(rows, columns);
for col = 1 : columns
	% Get random order of rows.
	randRows = randperm(rows);
	% Pick out "onesPerColumn" rows that will be set to 1.
	rowsWithOne = randRows(1:onesPerColumn);
	% Set those rows only to 1 for this column.
	m(rowsWithOne, col) = 1;
end
% create connection matrix
connM = sparse(m);
figure; spy(connM);title('added edge noise');

% bug fix; added random noise
fullN = full(connN);
m = m + fullN;
connM = sparse(m);

% Random walk
Ind = 1;
F = @(t,X) 0;
G = @(t,X) 1;
S = sde(F,G,'startState',Ind);
X = S.simByEuler(1000,'ntrials',1,'Z',@(t,X) RandomGraphMove(X,connM)-X);

% Figure for connectivity matrix
spos = get(0,'screensize');
figure('position',[spos(3)/10 spos(4)/10 spos(3)*.8 spos(4)*.8],'color','w');
gplot(connM,coordsM); hold on
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
weightsM = zeros(size(m));

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
