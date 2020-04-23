function plot_laser_rs_PF_ACG (filename,binsize,smooth)


load([filename '.mat'])

eventTime=laserON;
TimeWindow=[-2 3];
%  binsize=0.1;
% Celln=26;
% smooth=10;

for Celln=1:length(G);
figure('Position',[300 0 1000 700])
% subplot(4,2,[1 3 6])
spikeTime=selectgroup(spk.ts,spk.g,G(Celln));
Tperievent(spikeTime,eventTime,TimeWindow,binsize);

subplot(4,2,[2 4 6])
 plotPFSetObj_DJ(filename,Celln,smooth)

 subplot(4,2,8)
 bar(acg.acg1t,acg.acg1(:,Celln))
xlim([-50 50]); %xlim([-timewidth timewidth]);
if length(G_C)==length(G);
title(['Cell ' num2str(Celln)]) 
else
    KK=find(G(Celln)==G_C);
    title(['Cell ' num2str(KK)]) 
end
cd laser_raster
print([filename '_' num2str(G(Celln))], '-depsc' )
cd ..
close all



end
