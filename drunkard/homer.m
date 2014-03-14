%% Random walk visualizations

% Set up parameters.
rows = 10;
columns = 10;
onesPerColumn = 4;
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
% Display m

% set one node to default to 1 as starting point
m(1,1)=1;

% create connection matrix
connM = sparse(m);
figure; spy(connM);

% Draw points on circle
A      = 1;  %starting point
p      = 11; %number of points
angles = linspace(0,2*pi,p); 
f      = @(theta) [cos(theta),sin(theta)]; 
r      = arrayfun(@(i)f(angles(i)),1:p,'UniformOutput',false); 
coordsM     = reshape(cell2mat(r),[2 p]); 
figure
line(coordsM(1,:),coordsM(2,:)); 
hold on 
plot(coordsM(1,:),coordsM(2,:),'ro');
coordsM(:,end)=[];
coordsM = coordsM';

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
    pause(0.0001);
    drawnow
end
hold on

% 
figure
imagesc(weightsM)




% Using wgplot
figure;
[he,hv]=wgPlot(weightsM,coordsM);
