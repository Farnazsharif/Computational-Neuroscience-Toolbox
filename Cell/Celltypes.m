
load('FM05_1.mat')
load('PrePostG1.mat')
%%
[~,prendx]=ismember(PrePostG1(:,1),G);
[~,postndx]=ismember(PrePostG1(:,2),G);
%%
figure
plot(spkinfo.duration, acg.burstIndex,'x','markersize',10,'LineWidth',2) 
hold on
plot(spkinfo.duration(postndx), acg.burstIndex(postndx),'rd','markersize',16,'LineWidth',2) 
hold on
plot(spkinfo.duration(prendx), acg.burstIndex(prendx),'go','markersize',16,'LineWidth',2) 


xlabel('spk.duration','fontsize',16) 
ylabel('burstIndex','fontsize',16)
%%
figure
plot(acg.refractoryT, acg.burstIndex,'x','markersize',10,'LineWidth',2)
hold on
plot(acg.refractoryT(1,prendx), acg.burstIndex(1,prendx),'go','markersize',16,'LineWidth',2) 
hold on
plot(acg.refractoryT(1,postndx), acg.burstIndex(1,postndx),'rd','markersize',16,'LineWidth',2) 
xlabel('spk.refractoryT','fontsize',16)
ylabel('burstIndex','fontsize',16)
%%
a=acg.refractoryT;
b=acg.burstIndex;
c=spkinfo.duration;
%%
figure
plot3(spkinfo.duration,acg.refractoryT,acg.burstIndex,'x','markersize',10,'LineWidth',1.5) 
grid on
hold on
plot3(spkinfo.duration(prendx),acg.refractoryT(prendx),acg.burstIndex(prendx),'go','markersize',16,'LineWidth',2) 
hold on
plot3(spkinfo.duration(postndx),acg.refractoryT(postndx),acg.burstIndex(postndx),'rd','markersize',16,'LineWidth',2) 
hold on
xlabel('spkinfo.duration','fontsize',16)
ylabel('refractoryT','fontsize',16)
zlabel('burstIndex','fontsize',16)
 set(gcf,'Color','w');

%%
figure
plot3(SpatialC(:,set+1),spkinfo.duration,acg.burstIndex,'x','markersize',10,'LineWidth',1.5) 
grid on
hold on
plot3(SpatialC(prendx,set+1),spkinfo.duration(prendx),acg.burstIndex(prendx),'go','markersize',16,'LineWidth',2) 
hold on
plot3(SpatialC(postndx,set+1),spkinfo.duration(postndx),acg.burstIndex(postndx),'rd','markersize',16,'LineWidth',2) 
hold on
xlabel('SpatialC','fontsize',16)
ylabel('spkinfo.duration','fontsize',16)
zlabel('burstIndex','fontsize',16)


%%
filename3='FM05_1PlaceField';
filename1='FM05_1';
[SpatialC]=spatialcoverage2(filename1,filename3);
%%
 figure;

subplot(4,1,1)
% 
set=2;
plot(SpatialC(:,set+1), acg.burstIndex,'x','markersize',10,'LineWidth',2) 
hold on
plot(SpatialC(postndx,set+1), acg.burstIndex(postndx),'rd','markersize',16,'LineWidth',2) 
hold on
plot(SpatialC(prendx,set+1), acg.burstIndex(prendx),'go','markersize',16,'LineWidth',2) 
% set(gcf,'Color','w');
xlabel('SpatialC','fontsize',16)
ylabel('burstIndex','fontsize',16)
title('Set2','fontsize',16)

subplot(4,1,2)
plot(SpatialC(:,set+1), acg.refractoryT,'x','markersize',10,'LineWidth',2) 
hold on
plot(SpatialC(postndx,set+1), acg.refractoryT(postndx),'rd','markersize',16,'LineWidth',2) 
hold on
plot(SpatialC(prendx,set+1), acg.refractoryT(prendx),'go','markersize',16,'LineWidth',2) 
xlabel('SpatialC','fontsize',16)
ylabel('refractoryT','fontsize',16)

subplot(4,1,3)
plot(SpatialC(:,set+1), spkinfo.duration,'x','markersize',10,'LineWidth',2) 
hold on
plot(SpatialC(postndx,set+1), spkinfo.duration(postndx),'rd','markersize',16,'LineWidth',2) 
hold on
plot(SpatialC(prendx,set+1), spkinfo.duration(prendx),'go','markersize',16,'LineWidth',2) 
xlabel('SpatialC','fontsize',16)
ylabel('spkinfo.duration','fontsize',16)

subplot(4,1,4)
plot(SpatialC(:,set+1), Frate.rest,'x','markersize',10,'LineWidth',2) 
hold on
plot(SpatialC(postndx,set+1), Frate.rest(postndx),'rd','markersize',16,'LineWidth',2) 
hold on
plot(SpatialC(prendx,set+1), Frate.rest(prendx),'go','markersize',16,'LineWidth',2) 
xlabel('SpatialC','fontsize',16)
ylabel('Frate.rest','fontsize',16)

%%






























