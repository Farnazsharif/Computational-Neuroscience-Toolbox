function [Rf_ch_shank]=Ref_ch_finder(filename,Winrange,Rt)
% filename='FM05_1';
% Winrange=100;
% Rt=Rt4;

%%
load([filename '.mat']);
eeg= readmulti([filename '.lfp'], 128);
EEGsamplerate=(spkinfo.samplerate./25);
%%
[R_n,~]=size(Rt);
for chn=1:128
    eegF=eegfilt(eeg(:,chn)',EEGsamplerate,120,180); 
for i=1:R_n
    eg(i,:)=eegF(Rt(i,chn)-Winrange:Rt(i,chn)+Winrange);
end
    eegT(chn,:)=mean(eg);
    
end

a=max(abs(eegT),[],2);
for i=1:16
[s(i),Rf_ch_shank(i)]=max(a(i*8-7:i*8));
end

%% ploting ripple per shanks
% chorder=[1 8 2 7 3 6 4 5];
% nsite=8;
% nshk_bgn=1;
% nshk_end=128;
% reorder=[];
% for ii = 1:16
%     reorder = [reorder chorder+(ii-1)*nsite];
% end
% M=reorder(nshk_bgn:nshk_end);
% W=reshape(M,8,16);
% 
% %%
% eegTF=eegT(reorder,:);
% 
% for i=1:8
%     subplot(1,8,i)
% imagesc(eegTF(i*8-7:i*8,:))
% hold on
% end












