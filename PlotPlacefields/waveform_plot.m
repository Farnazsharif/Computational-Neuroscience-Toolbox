function waveform_plot(filename,Celln)
% Waveform_plot('DJ57_S2_2Oc',Celln)
% filename='SFA5_S4_TRD1';;
% Celln=231;
%%

load([filename '.mat'])
% G_C=G;% for acute
shankID=fix(G_C(Celln)/100);
subg=mod(G_C(Celln),100);
%%
shankID=2;
subg=9;
%%
if isnan(shankID)==0
    
load([filename '_clu_',num2str(shankID),'.mat']);
ndex=find(clu.g==subg);

%% load spkw
chorder=1:10; 
nsamples=60;% chronic probe
% nsamples=80;% for acute
nchannels = length(chorder);
fp = fopen([filename '.spk.' num2str(shankID)], 'r');
spkW = fread(fp, [nchannels, inf], 'short');
nspike = size(spkW, 2)/nsamples;
spkW = reshape(spkW, [nchannels, nsamples, nspike]);
spkW = spkW(chorder,:,:);
%% plot SPK position
if length(ndex) < 50
    sndx=length(ndex);
else
    sndx=50;
end

chspace=700;% for chronic
% chspace=150;% for acute

for kk=1:10;

 plot(squeeze(spkW(kk,1:60,ndex(1:sndx)))-(kk-1)*chspace,'color',[ 1 0.7 0.4 ],'linewidth',2.5) % chronic probs
% plot(squeeze(spkW(kk,1:60,ndex(1:sndx)))-(kk-1)*chspace,'k') % acute
hold on

end
% hold on
% figure
% for rr=1:10;
% %  
% plot(squeeze(clu.spk_mean(rr,10:50,subg))-(rr-1)*chspace,'color',[ 1 0.7 0.4 ],'linewidth',1.5) 
% hold on
% 
% end

% set(gca,'YTickLabel',[])
% set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[10:-1:1])
% set(gca,'Ytick',[ -9*chspace -8*chspace -7*chspace -6*chspace -5*chspace -4*chspace -3*chspace -2*chspace -1*chspace 0 ])
% ylabel('Ch # ','fontsize',11)
% xlabel('Sample #','fontsize',11)
xlim([0 60]) 
%  ylim([-10*chspace-80 chspace/2])% for acute
ylim([-10*chspace chspace/2])% for chronic
end
