function [E] = computeHash(F, debug_, P)
clear HashTable
global HashTable
%workaround initialization of map structure
HashTable=containers.Map(uint32(1), [1 1]);
remove(HashTable,1);

%% parameters
%search window in frames
searchWindow=3;

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

tic
E=zeros(nNeurons);
timeStacked=zeros(nNeurons, nNeurons, searchWindow+1);
for i=1:nSamples-searchWindow
    [r, c]=find(booleanSpikes(i:i+searchWindow,:));
    firstNeuron=find(booleanSpikes(i,:));
    sampleIndex=i;
    for j=1:length(firstNeuron)
        for secondNeuron=1:length(c)
            %the -1 is to not skip hash values 1 through nNeurons
            dt=r(secondNeuron);
            hash=(firstNeuron(j)-1)*(nNeurons)^2+(c(secondNeuron)-1)*(nNeurons)+dt;
            E(c(secondNeuron),firstNeuron(j))=E(c(secondNeuron),firstNeuron(j))+1/dt/(1+norm(P(firstNeuron(j),:)-P(c(secondNeuron),:)))^(1/3);
            timeStacked(c(secondNeuron),firstNeuron(j),dt)=timeStacked(c(secondNeuron),firstNeuron(j),dt)+1;
%             try 
%                 HashTable(hash) = [HashTable(hash); [sampleIndex dt]];
%             catch me
%                 HashTable(hash)= [sampleIndex dt];
%             end
        end
    end
end
toc
disp('sie')
tic
keySet=keys(HashTable);
C=zeros(1,size(keySet,2));
for i=1:size(keySet,2)
    A=values(HashTable, keySet(i));
    B=cellfun(@size,A,'uni',false);
    C(i)=B{1}(1);
end
toc
save('mrtable.mat','HashTable')
D=C./max(C);
E=E./max(max(max(E)));
% A=diag(diag(E));

%setting neuron1=neuron2 counts to zero improved the score a lot.
E(logical(eye(size(E)))) = 0;
for i=1:searchWindow
    timeStackedTemp=timeStacked(:,:,i);
    timeStackedTemp(logical(eye(size(timeStacked(:,:,1))))) = 0;
    timeStacked(:,:,i)=timeStackedTemp;
end


% add_timeSeries(datas,sr, 1,handles);

end

