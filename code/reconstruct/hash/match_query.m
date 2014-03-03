function [R,L] = match_query(D,SR,handles,fp)
% [R,L] = match_query(D,SR)
%     Match landmarks from an audio query against the database.
%     Rows of R are potential maxes, in format
%      songID  modalDTcount modalDT
%     i.e. there were <modalDTcount> occurrences of hashes 
%     that occurred in the query and reference with a difference of 
%     <modalDT> frames.
%     L returns the actual landmarks that this implies.
% 2008-12-29 Dan Ellis dpwe@ee.columbia.edu
global HashTable
global totalFFTtime
global recurrencePlot
clear Lmarks
global Lmarks
global fpLm
global fpLmbox
%Rt = get_hash_hits(landmark2hash(find_landmarks(D,SR)));
%Lq = find_landmarks(D,SR,handles);
if nargin < 4
    fp=0;
end

%Lq = fuzzify_landmarks(Lq);

% Augment with landmarks calculated half-a-window advanced too
% landmarks_hopt = (handles.EEG.times(end)-handles.EEG.times(1))/totalFFTtime(1,2)/1000;
% landmarks_hopt=

%max number of matches to find per channel/and for recurrence plot
perChan=10;
rperChan=10;

%how many dt blocks to meld close matches
meldWindow=1;

%allow all landmarks in 2*lmBox+1 side length box
lmBox=0;

%this many rows counting up from the bottom will not be made into boxes (box creates too many false matches at lower frequencies)
boxCutoff=13;

%check block around initial 2nd landmark point, except for fingerprint (already done)
% if fp==1 && size(fpLm,1)~=0
%     Lq=fpLm;
% else
Lq=[];
preLq= find_landmarks2(D,SR,handles);
for h=1:size(preLq,1)
    if preLq(h,4)~=1
        for j=-lmBox:lmBox
            for k=-lmBox:lmBox
                Lq(end+1,1)=preLq(h,1);
                if preLq(h,2)+j>boxCutoff && preLq(h,2)+j<54
                    Lq(end,2)=preLq(h,2)+j;
                else
                    Lq(end,2)=preLq(h,2);
                end
                if preLq(h,3)+k>boxCutoff && preLq(h,3)+k<54
                    Lq(end,3)=preLq(h,3)+k;
                else
                    Lq(end,3)=preLq(h,3);
                end                
                Lq(end,4)=preLq(h,4);
                Lq(end,5)=preLq(h,5);
                Lq(end,6)=preLq(h,6);
            end
        end
    end
end
    %     preLq=[preLq;find_landmarks(D(i,round(landmarks_hopt/4*SR):end),SR,handles)];
    %     preLq=[preLq;find_landmarks(D(i,round(landmarks_hopt/2*SR):end),SR,handles)];
    %     preLq=[preLq;find_landmarks(D(i,round(3*landmarks_hopt/4*SR):end),SR,handles)];
    %     preLq=[preLq;find_landmarks(D(i,round(landmarks_hopt*SR):end),SR,handles)];
    %     preLq=[preLq;find_landmarks(D(i,round(5*landmarks_hopt/4*SR):end),SR,handles)];
    %     preLq=[preLq;find_landmarks(D(i,round(6*landmarks_hopt/4*SR):end),SR,handles)];
    %     preLq=[preLq;find_landmarks(D(i,round(7*landmarks_hopt/4*SR):end),SR,handles)];
% end
%Lq=unique(Lq,'rows');
%Lq = [Lq;find_landmarks(D(round(landmarks_hopt/4*SR):end),SR)];
%Lq = [Lq;find_landmarks(D(:,round(landmarks_hopt/2*SR):end),SR,handles)];
%Lq = [Lq;find_landmarks(D(round(3*landmarks_hopt/4*SR):end),SR)];
% add in quarter-hop offsets too for even better recall

if size(Lq,1)==1
    Rt=[];
else
    Hq = landmark2hash(Lq);
    Rt = get_hash_hits(Hq);
end
nr = size(Rt,1);

if nr > 0
tkR=Rt;
offsets=[];
j=1;
if recurrencePlot==1
    perChan=rperChan;
end
while j<=perChan
    if j>1
        rows_to_remove = any(tkR==dts(xx), 2);
        tkR(rows_to_remove,:) = [];
    end
    [dts,xx] = unique(sort(tkR(:,2)),'first');
    dtcounts = 1+diff([xx',size(tkR,1)]);
    [vv,xx] = max(dtcounts);
    R(j,:) = [1,vv,dts(xx)];
    doubles=0;
    for k=1:j-1
        if abs(R(j-doubles,3)-R(k,3))<meldWindow
            offsets=[offsets; R(k,3) R(j-doubles,3)];
            R(k,2)=R(k,2)+R(j-doubles,2);
            R(j-doubles,:)=[];
            doubles=doubles+1;
        end
    end
    j=j+1-doubles;
    if length(tkR)<1
        j=perChan+1;
    end
end

  % Sort by descending match count
  [vv,xx] = sort(R(:,2),'descend');
  R = R(xx,:);

  % Extract the actual landmarks
  for m=1:perChan
      if size(R,1)>=m
          H = Rt(Rt(:,2)==R(m,3),:);
          try
            extras=offsets(offsets(:,1)==R(m,3),2);
          catch me
              extras=[];
          end
          for o=1:length(extras)
              H=[H;Rt(Rt(:,2)==extras(o),:)];
          end
          % Restore the original times
          for i = 1:size(H,1)
            hix = find(Hq(:,3)==H(i,3));
            hix = hix(1);  % if more than one...
            H(i,2) = H(i,2)+Hq(hix,2);
            L(i,:) = hash2landmark(H(i,:));
            if R(m,3)==0
                R(m,3)=1;
            end
            Lmarks(i,:,m)=[L(i,:) R(m,3) H(i,1)];
           end


      % Return no more than 10 hits, and only down to half the #hits in
      % most popular
%       if size(R,1) > 10
%         R = R(1:10,:);
%       end
          maxhits = R(m,2);
          nuffhits = R(:,2)>(maxhits/2);
      end
      %R = R(nuffhits,:);
  end
else
  R = [];
  L = [];
  disp('*** NO HITS FOUND ***');
end
recurrencePlot=0;
end
