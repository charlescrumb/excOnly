function [L,S,T,maxes] = find_landmarks(D,SR,handles,ID)
%%
% Includes:
%   - weightedCentroid based landmarking
%   - peakfinding based on fastPeakFinder
%

if nargin < 4
    ID=0;
end
%SETTINGS 

%prevent landmark peaks within this many blocks
barrier=10;

% Set initial threshold envelope based on peaks in first x frames. consider
% increasing or randomizing to have this more based on average signal source energy 
% changed from default of 10 -->
numFramesThresh = 15; 


% Limit the number of pairs that we'll accept from each peak
maxpairsperpeak=6;   % moved to front by DAn

%% overmasking factor?  1.0 --> 0.95; increases sensitivity
s_sup = 0.995;


%%%%%%%%%%%%%%%%%%%%%%%%%%

% L = find_landmarks(D,SR)
%   Make a set of spectral feature pair landmarks from some audio data
%   D is an audio waveform at sampling rate SR
%   L returns as a set of landmarks, as rows of a 4-column matrix
%   {start-time-col start-freq-row end-freq-row delta-time]
%
%  REVISED VERSION FINDS PEAKS INCREMENTALLY IN TIME WITH DECAYING THRESHOLD
% 
% 2008-12-13 Dan Ellis dpwe@ee.columbia.edu

% The scheme relies on just a few landmarks being common to both
% query and reference items.  The greater the density of landmarks,
% the more like this is to occur (meaning that shorter and noisier
% queries can be tolerated), but the greater the load on the
% database holding the hashes.
%
% The factors influencing the number of landmarks returned are:
%  A.  The number of local maxima found, which in turn depends on 
%    A.1 The spreading width applied to the masking skirt from each
%        found peak (gaussian half-width in frequency bins).  A
%        larger value means fewer peaks found.
%30 to
f_sd = 10;

%    A.2 The decay rate of the masking skirt behind each peak
%        (proportion per frame).  A value closer to one means fewer
%        peaks found.
%.99 to
a_dec = 0.93;

%    A.3 The maximum number of peaks allowed for each frame.  In
%        practice, this is rarely reached, since most peaks fall
%        below the masking skirt
maxpksperframe = 5;

%    A.4 The high-pass filter applied to the log-magnitude
%        envelope, which is parameterized by the position of the
%        single real pole.  A pole close to +1.0 results in a
%        relatively flat high-pass filter that just removes very
%        slowly varying parts; a pole closer to -1.0 introduces
%        increasingly extreme emphasis of rapid variations, which
%        leads to more peaks initially.
%from .98 to
hpf_pole = 0.93;

%  B. The number of pairs made with each peak.  All maxes within a
%     "target region" following the seed max are made into pairs,
%     so the larger this region is (in time and frequency), the
%     more maxes there will be.  The target region is defined by a
%     freqency half-width (in bins)
targetdf = 31;  % +/- 50 bins in freq (LIMITED TO -32..31 IN LANDMARK2HASH)

%     .. and a time duration (maximum look ahead)
targetdt = 30;  % (LIMITED TO <64 IN LANDMARK2HASH)

%     The actual frequency and time differences are quantized and
%     packed into the final hash; if they exceed the limited size
%     described above, the hashes become irreversible (aliased);
%     however, in most cases they still work (since they are
%     handled the same way for query and reference).

%to translate time to EEG time
global totalFFTtime
D=double(D);
verbose = 0;

% % Convert D to a mono row-vector
% [nr,nc] = size(D);
% if nr > nc
%   D = D';
%   [nr,nc] = size(D);
% end
% if nr > 1
%   D = mean(D);
%   nr = 1;
% end

%Resample to target sampling rate
targetSR = SR;

%# Check to see if multivariate data
if length(handles.EEG.chanlocs) > 1

    %# Determine if fingerprint or not
    if size(D,2)==size(handles.EEG_fp.data,2)
        datas = handles.EEG_fp.data;
    else %# If not, load sample data
        
    datas = handles.EEG.data;
    end
    
    nmaxes3 = 0;
    maxes3 = [0 0 0]';
    
    %# Loop through channels and add maxes to nmaxes2
    for h=1:length(handles.EEG.chanlocs)
        disp('Channel:');
        disp(h);
        [S,freq] = morFinger(datas(h,:),SR,handles,ID);

        % for highRez
        rezBoxes = 32; % resolution; how many boxes in spectra. changed from 32-->128

        %% DOWNSAMPLE MATRIX
        %% newlen should be relative to size of the data sample!... so fingerprint and source will be different
        % make multiples of 100
        lengthD = size(D,2); %length of sample
        % disp('Old length of D:');
        % disp(lengthD);
        lengthSig = lengthD/SR; %length of sample rate in sec
        matrix = S';
        newlen = lengthSig*rezBoxes; %downsampled from 2400 --> 100 (factor of ~32)

        [rows cols]=size(matrix);
        newmatrix=[];
        for i=1:cols
        len=rows;
        x=1:len;x=x';y=matrix(:,i);
        xx=1:(len-1)/(newlen-1):len;
        xx=xx';
        yy=interp1(x,y,xx,'linear');
        newmatrix=[newmatrix yy];
        end

        S = round(newmatrix'); %store downsampled sig ALSO /8 to reduce frequencies back to prior range

        %%move to positive domain
            %strategy one -- remove all negatives
        %     zeroS = S<=0;
        %     S(zeroS) = 0;

            %strategy two -- add lowest negative number
            Smin = min(min(S));
            S = S + abs(Smin);

        %%

        if ID==1 && lengthD==length(handles.EEG.data)
            totalFFTtime=size(S);
        %     disp('totalFFTtime assigned:');
        %     disp(totalFFTtime);
        end

        %%%%% FILTERS

        %% OLD HPF 
        % S=abs(specgram(double(D),60,targetSR,20));
        % % convert to log domain, and emphasize onsets %(needed for linear CWT)
        % Smax = max(S(:));                              %(needed for linear CWT)
        % % Work on the log-magnitude surface    %(needed for linear CWT)
        % S = log(max(Smax/1e6,S));               %(needed for linear CWT)
        % Make it zero-mean, so the start-up transients for the filter are
        % minimized
        % S = S - mean(S(:)); % removed b/c morLet already zero-meaned

        %     S = (filter([1 -1],[1 -hpf_pole],S')');
        %     sThresh = 7; %Removes values below this number (like fill-flood)
        %     S(S<=sThresh) = 0;

        %% 
        % GAUSIAN LP FILTER
        %     filt = (fspecial('gaussian', [2 2],0.5));
        %     S = imfilter(S,filt);

        %% 2d FFT Filter
        %
        %Determine good padding for Fourier transform
        PQ = paddedsize(size(S));
        %Create a Gaussian Highpass filter 5% the width of the Fourier transform
        D0 = 0.10*PQ(1);
        H = hpfilter('gaussian', PQ(1), PQ(2), D0);
        % Calculate the discrete Fourier transform of the image
        F=fft2(double(S),size(H,1),size(H,2));
        % Apply the highpass filter to the Fourier spectrum of the image
        HPFS_S = H.*F;
        % convert the result to the spacial domain.
        HPF_S=real(ifft2(HPFS_S));
        % Crop the image to undo padding
        HPF_S=HPF_S(1:size(S,1), 1:size(S,2));

        % Store HPF_S
        S = HPF_S;
        % REMOVE values less than a threshold (<0)
        zeroS = S<=0;
        S(zeroS) = 0;

        % REMOVE LOW FREQS < 5 Hz
        S(1:10,:) = 0;

        % REMOVE LOW FREQS > 100 Hz
        S(54:64,:) = 0;


%         %%
%         %%% OLD PEAK FINDER
%         %Estimate for how many maxes we keep - < 30/sec (to preallocate array)
% 
% 
%         maxespersec = 30;
% 
%         ddur = length(D)/targetSR;
%         nmaxkeep = round(maxespersec * ddur);
%         maxes = zeros(4,nmaxkeep); % added column to maxes for Channel
%         nmaxes = 0;
%         maxix = 0;
% 
%         % find all the local prominent peaks, store as maxes(i,:) = [t,f];
% 
%         % initial threshold envelope based on peaks in first 10 frames
%         sthresh = s_sup*spread(max(S(:,1:numFramesThresh),[],2),f_sd)';
% 
%         % T stores the actual decaying threshold, for debugging
%         T = 0*S;
% 
%         for i = 1:size(S,2)-1
%           s_this = S(:,i);
%           sdiff = max(0,(s_this - sthresh))';
%           % find local maxima
%           sdiff = locmax(sdiff);
%           % (make sure last bin is never a local max since its index
%           % doesn't fit in 8 bits)
%           sdiff(end) = 0;  % i.e. bin 257 from the sgram
%           % take up to 5 largest
%           [vv,xx] = sort(sdiff, 'descend');
%           % (keep only nonzero)
%           xx = xx(vv>0);
%           % store those peaks and update the decay envelope
%           nmaxthistime = 0;
%           for j = 1:length(xx)
%             p = xx(j);
%             if nmaxthistime < maxpksperframe
%               % Check to see if this peak is under our updated threshold
%               if s_this(p) > sthresh(p)
%                 nmaxthistime = nmaxthistime + 1;
%                 nmaxes = nmaxes + 1;
%                 maxes(2,nmaxes) = p;
%                 maxes(1,nmaxes) = i;
%                 maxes(3,nmaxes) = s_this(p);
%                 eww = exp(-0.5*(([1:length(sthresh)]'- p)/f_sd).^2);
%                 sthresh = max(sthresh, s_this(p)*s_sup*eww);
%               end
%             end
%           end
%           T(:,i) = sthresh;
%           sthresh = a_dec*sthresh;
%         end
% 
%         % Backwards pruning of maxes
%         maxes2 = [];
%         nmaxes2 = 0;
%         whichmax = nmaxes;
%         sthresh = s_sup*spread(S(:,end),f_sd)';
%         for i = (size(S,2)-1):-1:1
%           while whichmax > 0 && maxes(1,whichmax) == i
%             p = maxes(2,whichmax);
%             v = maxes(3,whichmax);
%             if  v >= sthresh(p)
%               % keep this one
%               nmaxes2 = nmaxes2 + 1;
%               maxes2(:,nmaxes2) = [i;p;h];
%               eww = exp(-0.5*(([1:length(sthresh)]'- p)/f_sd).^2);
%               sthresh = max(sthresh, v*s_sup*eww);
%             end
%                % if  v >= sthresh(p) && (abs(i-maxes2(1,nmaxes2))>barrier || abs(p-maxes2(2,nmaxes2))>barrier)
%             whichmax = whichmax - 1;
%           end
%           sthresh = a_dec*sthresh;
%         end
% 
%         maxes2 = fliplr(maxes2);
%         %# concatenate maxes2 into maxes3
%         
%         maxes3 = horzcat(maxes3,maxes2);
%         nmaxes3 = nmaxes3 + nmaxes2;
% 
%         %% END OF OLD PEAKFINDER
%         
         %NEW Centroid finder    
        filt = (fspecial('gaussian', [3 3],0.5));
        p=FastPeakFind(S',18,filt,10,2);  
%         %For Debugging morlet code
%         if h==20;
%             disp('min S');
%             disp(min(min(S)));
%             disp(max(max(S)));
%             disp(mean(mean(S))*0.1);
%           
%             figure
%               imagesc([],[],(40*log(abs(S))));  % set contrast; adjusted 25-->
%               colormap(1-hot); %bluescale
%               axis xy
%               axis tight
%               pause
%         end
        
        nmaxes2=size(p,1)/2;
        maxes2=zeros(2,nmaxes2);
        maxes2=reshape(p,2,nmaxes2)';   

        %Swap columns
        maxes2 = fliplr(maxes2);
        %add column with channel info
        maxes2 = horzcat(maxes2,repmat(h,size(maxes2,1),1))';  
        
        %# concatenate maxes2 into maxes3
        maxes3 = horzcat(maxes3,maxes2);
        nmaxes3 = nmaxes3 + nmaxes2;
        
        
        
    end %end of multichannel for loop

else
        nmaxes3 = 0;
    
end%end of multichannel for loop



%% New peakfinder; 
% 
% %Centroid finder    
% filt = (fspecial('gaussian', [3 3],0.5));
% p=FastPeakFind(S',17,filt,10,2);
% 
% disp('show landmarks');
% disp(size(p));
% 
% nmaxes2=size(p,1)/2;
% maxes2=zeros(2,nmaxes2);
% maxes2=reshape(p,2,nmaxes2)';   
%
% Swap columns
% maxes2 = fliplr(maxes2)';

%% For Debugging
% %for debugging landmarking peakFinding system compared to output
%   set(handles.figure1,'CurrentAxes',handles.axes3);
%   tbase=1;
%   [nr,nc] = size(S);
%   tt = [1:nc]*tbase;
%   
%   imagesc(tt,[],(50*log10(S)));
%   axis tight
%   hold on
%   filt = (fspecial('gaussian', 2,4));
%     q=FastPeakFind(S,25,filt,5,2);
%     
%     nmaxesT=size(p,1)/2;
%     maxesT=zeros(2,nmaxesT);
%     maxesT=reshape(q,2,nmaxesT)';    
%     
%     plot(q(2:2:end)*tbase,q(1:2:end),'r.');  % red dots are what packingMaxes looks for!!  ; it's matrix is rotated 
%     hold on
%     plot(maxesT(:,1),maxesT(:,2),'g+');   
%     hold off
%   set(handles.figure1,'CurrentAxes',handles.axes2);
% 
%     disp(maxes2);
%     disp(nmaxes2);
%     disp(q)
%     disp(maxes2);
%     disp(nmaxes2);

%% Pack the maxes into nearby pairs = landmarks
  
% Limit the number of pairs that we'll accept from each peak
% maxpairsperpeak=3;   % moved to front by DAn

% Landmark is <starttime F1 endtime F2>
L = zeros(nmaxes3*maxpairsperpeak,4);

% Store maxes3 into maxes2;
maxes2 = maxes3;
clear maxes3

nlmarks = 0;

for i =1:nmaxes3
  startt = maxes2(1,i);     
  F1 = maxes2(2,i);  
  maxt = startt + targetdt;
  minf = F1 - targetdf;
  maxf = F1 + targetdf;
  matchmaxs = find((maxes2(1,:)>startt)&(maxes2(1,:)<maxt)&(maxes2(2,:)>minf)&(maxes2(2,:)<maxf)&(abs(maxes2(2,:)-F1)>barrier|(maxes2(1,:)>startt+barrier)));
  if length(matchmaxs) > maxpairsperpeak
    % limit the number of pairs we make; take first ones, as they
    % will be closest in time
    matchmaxs = matchmaxs(1:maxpairsperpeak);
  end
  for match = matchmaxs
    nlmarks = nlmarks+1;
    L(nlmarks,1) = startt;
    L(nlmarks,2) = F1; 
    L(nlmarks,3) = maxes2(2,match);  % frequency row
    
    %if time column difference is >10, then store time column difference
    %instead as this pair is more likely within the same channel...
    if (maxes2(1,match)-startt)>10
        L(nlmarks,4) = maxes2(1,match)-startt;  % time column difference
    else
        L(nlmarks,4) = maxes2(3,match);  % Channel number  
    end
    
  end
end

L = L(1:nlmarks,:);

if verbose
  disp(['find_landmarks: ',num2str(length(D)/targetSR),' secs, ',...
      num2str(size(S,2)),' cols, ', ...
      num2str(nmaxes),' maxes, ', ...
      num2str(nmaxes3),' bwd-pruned maxes, ', ...
      num2str(nlmarks),' lmarks']);
end
  
% for debug return, return the pruned set of maxes
maxes = maxes2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = locmax(X)
%  Y contains only the points in (vector) X which are local maxima

% Make X a row
X = X(:)';
nbr = [X,X(end)] >= [X(1),X];
% >= makes sure final bin is always zero
Y = X .* nbr(1:end-1) .* (1-nbr(2:end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = spread(X,E)
%  Each point (maxima) in X is "spread" (convolved) with the
%  profile E; Y is the pointwise max of all of these.
%  If E is a scalar, it's the SD of a gaussian used as the
%  spreading function (default 4).
% 2009-03-15 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; E = 4; end
  
if length(E) == 1
  W = 4*E;
  E = exp(-0.5*[(-W:W)/E].^2);
end

X = locmax(X);
Y = 0*X;
lenx = length(X);
maxi = length(X) + length(E);
spos = 1+round((length(E)-1)/2);
for i = find(X>0)
  EE = [zeros(1,i),E];
  EE(maxi) = 0;
  EE = EE(spos+(1:lenx));
  Y = max(Y,X(i)*EE);
end


