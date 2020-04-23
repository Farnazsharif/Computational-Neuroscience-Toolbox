function [restT,runT]=allpower(filename,speedthreshold,durationthreshold,a,b)
%  allpower('FM05_1',5,2,1,150)

% speedthreshold=5;
% durationthreshold=2;
% filename='FM05_1';
% a=1;
% b=150;

pyrch=[115];
% pyrch=[66,67,70,71,73,74,79,80,81,82,88,89,90,96,105,106,111,112,113,114,115,119,120,122,123,126,127];
% pyrch=[69,77,78,83,84,85,86,87,91,92,93,94,95,98,99,100,103,107,108,109,110,116,117,118,124,125];
% pyrch=[65:128];
load([filename '.mat']);
eeg= readmulti([filename '.lfp'], 128);
[row,culomn]=size(eeg);
EEGsamplerate=(spkinfo.samplerate./25);
Teeg = (1:row)/EEGsamplerate;
[trial]=Trialfinder1(filename);
% speedthreshold=5;
% durationthreshold=2;
d=1;
% T=behav.runT;
% T=behav.restT;
% T=[behav.TXDTS(trial(:,1)) behav.TXDTS(trial(:,2))];
[restT,runT]=findrunrestT(behav.TXDTS(:,1),behav.TXDTS(:,3),speedthreshold,durationthreshold,d);
T=restT;
Ndx=[];
for j=1: length(T)
  [~, Ndx(j,1)]=min(abs(T(j,1)-Teeg));
  [~, Ndx(j,2)]=min(abs(T(j,2)-Teeg));
end
Ndx;

%% Method1 sacling based on frequency ( to study the theta we need to have a , b)
% tic
% 
% power1=[];
% PowerT1=0;
% 
% for channel=pyrch
% eegch=eeg(:,channel);
% 
%     for g=1:length(Ndx)
%     L = floor(0.5*(Ndx(g,2)- Ndx(g,1)));
%     [Pxx F1]= pwelch(eegch(Ndx(g,1):Ndx(g,2)),L,[],[a:b],EEGsamplerate);
%     power1(:,g)=Pxx;
%     end
%     PowerT1=power1+PowerT1;
%     
% end
% save  (['F1_115_all'],'F1')
% PowerT1=PowerT1./length(pyrch);
% save  (['Power1_115_all'],'PowerT1')
% toc
%% Method 2 wrong power ( here a ,b are not important, it is so fast and and lim is working)
% tic
% 
% power2=[];
% PowerT2=0;
% 
% for channel=pyrch
% eegch=eeg(:,channel);
% 
%     for g=1:length(Ndx)    
% %     [Pxx F2]= pwelch(eegch(Ndx(g,1):Ndx(g,2)),floor(0.5*EEGsamplerate),[],[],EEGsamplerate);
%     [Pxx F2]= pwelch(eegch(Ndx(g,1):Ndx(g,2)),floor(0.5*EEGsamplerate),[],[a:b],EEGsamplerate);
%     power2(:,g)=Pxx;
%     end
%     PowerT2=power2+PowerT2;
%     
% end
% save  (['F2_115_all'],'F2')
% PowerT2=PowerT2./length(pyrch);
% save  (['Power2_115_all'],'PowerT2')
% toc
%% Method 3  wavelet
% freqband=[1,150];
tic
PowerT3=[];
S=[];
freq = a:b;
scale = frq2scal(freq,EEGsamplerate);
eegch=eeg(:,115);
for g=1:length(Ndx) 
S = cwt(eegch(Ndx(g,1):Ndx(g,2)),scale,'morl');
PowerT3(g,:) = mean(envelop(S.*S),2)';
end

save  (['PowerT3_n_150_rest'],'PowerT3')
toc
%% Method 4 Spectrogrm
% [yo,fo,to,do]=spectrogram(eegch(Ndx(g,1):Ndx(g,2)),50,25,1:100,EEGsamplerate);


end
