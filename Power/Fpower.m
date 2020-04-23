%% power 
%close all
%1-make correlation btween time and trials, 2-make correlation for
%frequency
clc
clear
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
%Senv = envelop(S)

%% ca=[];
lfpf=[];
for e=1:length(Ndx)
    eegf=eeg(Ndx(e,1):Ndx(e,2));
 lfpf=cat(1,eegf,lfpf);
end

b1=floor(length(lfpf)./(2.9*EEGsamplerate))
b2=floor(2.9*EEGsamplerate)
seg=b1*b2
lfpf1=lfpf(1:seg);
B=reshape(lfpf1,[b1,b2]);
%%
%figure(2)
[s1 s2]=size(B)
Ndxl=1:b2:seg;
power=[];
for g=1:length(Ndxl)-1 
[Pxx F]= pwelch(lfpf1(Ndxl(g):Ndxl(g+1)),0.5*s2,[],[],EEGsamplerate);
power(:,g)=Pxx;
end
%%
% %figure(2)
power=[];
for g=2:2%length(Ndx)
    L = floor(0.5*(Ndx(g,2)- Ndx(g,1)))
[Pxx F]= pwelch(eeg(Ndx(g,1):Ndx(g,2)),L,[],[],EEGsamplerate);
power(:,g)=Pxx;
end
%%
imagesc(matnorm(power,1)')
xlim([ 0 170])
%%
% bb=[behav.TXDTS(trial(:,1),1) behav.TXDTS(trial(:,2),1) trial(:,3) trial(:,4)];

%%
% figure(1)
% imagesc(F,(1:length(Ndx)),matnorm(power,1)')
% xlim([ 0 80])
%%
%[trial]=Trialfinder(filename);
for i=1:length(behav.restT)%length(runT)   
%I(i) = find(behav.TXDTS(trial(:,1),1)< behav.restT(i,2)& behav.restT(i,2)<=behav.TXDTS(trial(:,2),1));
I(i) = find(behav.TXDTS(trial(:,1),1)< behav.restT(i,2)& behav.restT(i,2)<=behav.TXDTS(trial(:,2),1));
end

figure(1)
imagesc(F,I,matnorm(power,1)')

%%
figure(2)
plot (F , matnorm(power(:,160),1) , 'b')
xlim([ 0 100])
hold on
%%
hold on
plot (F , matnorm(power1,1) , 'k')
xlim([ 0 100])
hold on
% for i=1:length(runT)
%     tr(i)=find(behav.TXDTS(:,1)==runT(:,2));
% end
% bb=[behav.TXDTS(trial(:,1),1) behav.TXDTS(trial(:,2),1) trial(:,3) trial(:,4)];
% for i=1:length(runT)   
% I(i) = find(behav.TXDTS(trial(:,1),1)< runT(i,2)& runT(i,2)<=behav.TXDTS(trial(:,2),1));
% end







