function plotripples(filename,Rt,dur,N,ch,nshk_bgn,nshk_end)
%   dur=100;
%   ch=39;
% k=1;
% nshk_bgn=65;
% nshk_end=128;
% Rt=Rt4;
% filename='FM05_1';
% Rtseb=round(rippleT2*(spkinfo.samplerate/25));
Rt_Number=Rt(N,:);
Rt_chs=Rt(:,ch);
% Rt_chs=Rtseb;
%%
chorder=[1 8 2 7 3 6 4 5];
nsite=8;

reorder=[];
for ii = 1:16
    reorder = [reorder chorder+(ii-1)*nsite];
end
M=reorder(nshk_bgn:nshk_end);
W=reshape(M,8,8);
%% plot one ripple time for all channels
load([filename '.mat']);
eegsamplingrate=(spkinfo.samplerate./25);
eeg= readmulti([filename '.lfp'], 128);
% 
Ndx=[];
Ndx(:,1)=(Rt_Number-dur)';
Ndx(:,2)=(Rt_Number+dur)';
eegS=[];
eegF=[];
% 
for i=1:nshk_end
eegS(:,i)=eeg(Ndx(i,1):Ndx(i,2),i);
[B,A]=butter(4,[120 180]/eegsamplingrate*2);
eegF(:,i)=filtfilt(B,A,eegS(:,i));
end
% 
%%
figure('position',[20,0,3000,2000])  
for kk=1:8

for ll=1:8
site=W(kk,ll);

subplot(8,8,ll+(kk-1)*8)
% plot(eegS(:,site))
plot(eegF(:,site))
% hold on 
% plot(dur+1,eegF(dur+1,ii),'r*')
end
end
set(gcf,'color','w')

print('-depsc',['Ripple_' num2str(N) ])
%% plot ?
% Ndx=[];
% Ndx(:,1)=(Rt_chs-dur)';
% Ndx(:,2)=(Rt_chs+dur)';
% eegS=[];
% eegF=[];
% for ii=1:length(Rt_chs)
% eegS=eeg(Ndx(ii,1):Ndx(ii,2),ch);
% [B,A]=butter(4,[120 180]/eegsamplingrate*2);
% eegF(:,ii)=filtfilt(B,A,eegS);
% end
% 
% for kk=1:length(Rt_chs)
% subplot(11,5,kk)
% plot(eegF(:,kk))
% hold on 
% plot(dur+1,eegF(dur+1,kk),'r*')
% end
% set(gcf,'color','w')
% print('-depsc',['Re' num2str(k) ])

%%

