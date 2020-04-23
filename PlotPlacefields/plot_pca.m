function [mean_coeaff,std_coeaff]=plot_pca (filename,Celln,sample_win,time_win,recording)
%%  
% 
filename='DJ52_S3C';
Celln=int_celln(1);
sample_win=10;
time_win=2;
 recording='c';
%% chronic

if isempty(strfind(recording,'c'))==0
Elnb=1:6;
samplerate=30000;
nbchperEl=10;
chorder=1:10;
nsamples=60;
load([filename '.mat'])
else
% acute
Elnb=[1:6,8:12];
samplerate=50000;
nbchperEl=10;
chorder=1:10;
nsamples=80;
load([filename '.mat'])
G_C=G;
end

%% make spk.ts
[~,spk,waveform] = SpikeTimeWave3(filename,Elnb,nbchperEl);
spk.ts = spk.i/samplerate;
% make spike time
shankID=fix(G_C(Celln)/100)
subg=mod(G_C(Celln),100)   
load([filename '_clu_',num2str(shankID),'.mat'])
ndx_spk=find(clu.g==subg);
spk_time=selectgroup(spk.ts,spk.g,G_C(Celln));

% load spkw

nchannels = length(chorder);
fp = fopen([filename '.spk.' num2str(shankID)], 'r');
spkW = fread(fp, [nchannels, inf], 'short');
nspike = size(spkW, 2)/nsamples;
spkW = reshape(spkW, [nchannels, nsamples, nspike]);
spkW = spkW(chorder,:,:);

%% find peak channels
use_samples=1:nsamples;

for ch_touse=1:nbchperEl
    
tmp_spkW = reshape(spkW(ch_touse,use_samples,ndx_spk),length(ch_touse)*length(use_samples),length(ndx_spk)); 
[a(ch_touse),b(ch_touse)]=min(mean(tmp_spkW ,2));

end

[~,ch_touse]=sort(a);

%% make pca
% figure('Position',[300 0 800 1000])
% k is channnel peak number
for k=1%:2
    
use_samples=b(ch_touse(k))-sample_win:b(ch_touse(k))+sample_win;
tmp_spkW = reshape(spkW(ch_touse(k),use_samples,ndx_spk),length(use_samples),length(ndx_spk));  
[coeff,score,latent,tsquared,explained,mu]=pca(tmp_spkW);
% coeff=coeff/abs(max(coeff(:,1)));
% make mean and st coeff from pca

h=0:time_win*60:spk_time(end);

Time_ndx=[];
for i=1:length(h)
    [~,t_ndx(i)]=min(abs(h(i)-spk_time(:,1)));
end

[s1,s2]=size(t_ndx);

if isempty(t_ndx)==1
    mean_coeaff=NaN;
    std_coeaff=NaN;
    
else
    
    if  s2 < 3
        Time_ndx=t_ndx;
    else
        
        Time_ndx(:,1)=[t_ndx(1) t_ndx(2:length(t_ndx)-1)+1];
        Time_ndx(:,2)=[t_ndx(2:length(h))];
    end
    
    mean_coeaff=[];
    std_coeaff=[];
    diff_mean=[];
    
    for i=1:length(h)-1
        mean_coeaff=[mean_coeaff mean(coeff(Time_ndx(i,1):Time_ndx(i,2),1))];
        std_coeaff=[std_coeaff std(coeff(Time_ndx(i,1):Time_ndx(i,2),1))];
    end
end
%%

% coeff=coeff/ max(mean_coeaff);
%% plot 
% figure('position', [100 100 500 800] )
subplot(1,2,k)
 plot(coeff(:,1)/(max(mean_coeaff)),spk_time,'k.')
% % ylim([0 max(spk_time)])
% % xlabel(['Ch = ' num2str(ch_touse(k))], 'FontSize', 20)
% 
 T=3.599613200000000e+03;
% t_tr=2.6829e+03;
% 
 ylim([0 T])
% hold on
% plot([-1 2],[t_tr t_tr],'r');
% set(gca,'box','off')
% subplot(1,3,3) ploting the wave amplitude versuse  time
% plot(tmp_spkW(b(ch_touse(1)),:),spk_time,'.')
% xlim([-3000 0])
% ploting the erorbar
 tr_treadmil=22;
 tr_maze=8;
subplot(1,2,2)
errorbar([1:length(Time_ndx)]-0.5,mean_coeaff/(max(mean_coeaff)),std_coeaff/(max(mean_coeaff)),'k.')
% errorbar(1:length(Time_ndx),mean_coeaff),std_coeaff,'k.')
% hold on
% plot(1:length(Time_ndx),mean_coeaff,'k.','markersize',15)
% hold on
%  plot([tr tr],[0.002 0.008],'r');
set(gca,'CameraUpVector',[1,0,0],'YDir','reverse','XAxisLocation','top')
% xlim([0 length(h)])
 xlim([0 tr_treadmil+ tr_maze])
%  xlim([1-0.2 tr+0.2])
% ylim([min(mean_coeaff-std_coeaff)-0.1 max(mean_coeaff+std_coeaff)+0.1])%for normalizing to the max coeff
%  set(gca,'box','off', 'XTick',[],'YTick',[],'TickDir','out')

end



%% ploting the erorbar as shadow between 2 lines
% figure 
% x1=(1 :length(Time_ndx))';
% y1=(mean_coeaff+std_coeaff)';
% y2=(mean_coeaff-std_coeaff(end))';
% 
% x=[x1 ; flipud(x1)];
% y=[y1 ; flipud(y2)];
% fill(x,y,[.9 .9 .9],'linestyle','none')
% hold on
% 
% plot(1:length(Time_ndx),mean_coeaff,'k','Linewidth',2)
% hold on
% 
% plot(x1,y1,'k--',x1,y2,'k--')
% 
% set(gca,'CameraUpVector',[1,0,0],'YDir','reverse','XAxisLocation','top')
% xlim([0 length(h)])















