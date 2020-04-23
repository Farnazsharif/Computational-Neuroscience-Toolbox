%% function allpower=(filename,speedthreshold,durationthreshold,ch)
ch=115;
speedthreshold=5;
durationthreshold=2;
filename='FM05_1';


load([filename '.mat']);
eeg= readmulti([filename '.lfp'], 128,ch);
EEGsamplerate=(spkinfo.samplerate./25);
Teeg = (1:length(eeg))/EEGsamplerate;
[trial]=Trialfinder(filename);
% [restT,runT]=findrunrestT(behav.TXDTS(:,1),behav.TXDTS(:,3),speedthreshold,durationthreshold);
T=behav.runT;
%T=behav.restT;
%T=[behav.TXDTS(trial(:,1)) behav.TXDTS(trial(:,2))];
Ndx=[];
for j=1: length(T)
  [~, Ndx(j,1)]=min(abs(T(j,1)-Teeg));
  [~, Ndx(j,2)]=min(abs(T(j,2)-Teeg));
end
Ndx;
%% Making Trial vector from Time
for i=1:length(T) 
I(i) = find(behav.TXDTS(trial(:,1),1)< T(i,2)& T(i,2)<=behav.TXDTS(trial(:,2),1));
end
I=I+5;
