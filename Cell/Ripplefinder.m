% clc
% clear
% filename='LAMP8_013_190619_151810';
eeg= readmulti([filename '.lfp'], 64);  
%% Speed check For maze
% eeg= readmulti([filename '.lfp'], 128);
% EEGsamplerate=(spkinfo.samplerate./25);
Ref_ch=5;
rippleT = detectripples(eeg(:,Ref_ch),spkinfo.samplerate/30,500);

for i=1:length(rippleT)
[~,b(i)]=min(abs(rippleT(i)-TXYVC(:,1)));
end

rippleT=[rippleT TXYVC(b,4)];
rippleT(rippleT(:,2)>5,:)=[];
%% Speed check for TRD
% rippleT=ripplet.peaks;
speedsm=x2v(behav.TXDTS(:,3),behav.TXDTS(:,1),1);
TXDTSV=[behav.TXDTS speedsm];

for i=1:length(rippleT)
[~,b(i)]=min(abs(rippleT(i)-behav.TXDTS(:,1)));
end
rippleT=[rippleT TXDTSV(b,6)];
rippleT(rippleT(:,2)>5,:)=[];

%%

%  lfp = lfp for ripple channel. One column with V values (integers)
%  thresholds = [low high] in SD. typical [2 5]
% lfp=eeg(:,Ref_ch);
% eegsamplingrate=spkinfo.samplerate/30;
% timestamps = (0:(1/eegsamplingrate):(length(lfp)-1)/eegsamplingrate)'; % assumes sampling rate of  =1250 Hz
% ripplet = bz_FindRipples(lfp,timestamps,'thresholds',[2 5],'durations',[20 300],'passband',[100 250],'EMGThresh',0);

%%

[t]= Scan_multi([],[],[],eeg(:,:),spkinfo.samplerate/30,rippleT(:,1));
% [t]= Scan_multi(behav.TXDTS,[],[],eeg(:,:),spkinfo.samplerate/30,ripplet.peaks);

%%

a=[195.13 253.4 287.65 308.25];
 t=[t a];
% oredreing 17 24 18 23 19 22 20 21
%%
% Ripple=[RippleT_Gamma RippleT_multiple RippleT_Single];
Ripple=[];
save('Ripple','Ripple')
Ripple_Time=t;
save('Ripple','-append','Ripple_Time','Ref_ch')
% save('Ripple','-append','RippleT_Single','RippleT_multiple','RippleT_Gamma','Ref_ch')
save('Ripple','-append','Ripple_Time','Ref_ch')

%%

% eeg= readmulti([filename '.lfp'], 64);

load('Ripple.mat')
Rt=[];
channelN=64;
for k=1:length(Ripple_Time)
%     close all


k
 Rt(k,:)=Ripple_t(filename,Ripple_Time(k),80,spkinfo.samplerate/30,channelN);

% plotripples('FM05_1',Rt,40,k);
end

save(['Ripple.mat'],'-append','Rt')
toc


%% find bigest ripple per shank and plot ripple laypout

eeg= readmulti([filename '.lfp'], 64);              
ShankN=6;
ChN=10;
Winrange=100;
chorder=[1 2 3 4 5 6 7 8 9 10]; 
load('Ripple.mat')
EEGsamplerate=spkinfo.samplerate/30;

[Rf_ch_shank,Rf_ch_shank_pow,bigR]=Ref_ch_finder(eeg(:,1:60),Winrange,Rt,spkinfo.samplerate/30,ShankN,ChN,[1 2 3 4 5 6 7 8 9 10]);

save(['Ripple.mat'],'-append','Rf_ch_shank','Rf_ch_shank_pow','bigR')

Rf_ch_shank_pow(:,6)=Rf_ch_shank_pow(:,4);
% Rf_ch_shank_pow(2,6)=7;
% Rf_ch_shank_pow(6,6)=5;
save(['Ripple.mat'],'-append','Rf_ch_shank_pow')

%%
figure
Ripplelayout(eeg(:,1:60),spkinfo.samplerate/30,[10 9 8 7 6 5 4 3 2 1 ],Ripple_Time(40),3)

%%
channelN=64;
eegsamplingrate=spkinfo.samplerate/30;sample_win=50;rippleT2=Ripple_Time(k);
Ripple=[RippleT_Single,RippleT_multiple,RippleT_Gamma];
%%

%%
filename='FM05_1';
s=find(G==903);
binrange=100;
cellN=Cellinfo.Cellchn_T(s,1)
ch=Cellinfo.Cellchn_T(s,2)
Cellripple(filename,ch,cellN,Rt4,binrange )
 %% ripple plots
 filename='FM01_2';
   dur=100;
  ch=34;
% N=1;
nshk_bgn=1;
nshk_end=64;
for N=1:length(Ripple_Time)

plotripples(filename,Rt,dur,N,ch,nshk_bgn,nshk_end);

 close all
end


%%
% TracelayoutFar(eeg2,EEGsamplerate,[5 4 6 3 7 2 8 1],[-0.05 0.05]+rippleT2(10),2,1000) 

 
