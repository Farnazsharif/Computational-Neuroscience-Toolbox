
CellN=8;
clc
SpeedThreshold=0.04;
TX=behav.TXVt;
TX(TX(:,4)==0,:)=[];
TX(isnan(TX(:,2))==1,:)=[];
TX=TX(find(TX(:,3)>SpeedThreshold),:);
TX(:,2)=TX(:,2)./max(TX(:,2));
TX=TX(:,1:2);

if Phase.probeID==3
    if floor(CellN/100)<7
        
        thetaP(:,2)=Phase.thetaP_1;
    else
        thetaP(:,2)=Phase.thetaP_2;
    end
    
else
    
    if floor(CellN/100)<9
        
        thetaP(:,2)=Phase.thetaP_1;
    else
        thetaP(:,2)=Phase.thetaP_2;
    end
    
end
thetaP(:,1)=(0:length(thetaP)-1)/1250;

%% Plot

ses=1;%even
CellG=[925,827,819];
ses=1;%odd
CellG=[305];

ses=4;%even
CellG=[810 402];

ses=8;%even
CellG=[1521];

ses=2;
CellG=[924];
ses=5;
CellG=[931];
ses=7;
CellG=[502];


%%
cd 'Metrics'
% load('tiral_3')
load('tiral_3')
cd ..

%% Plot 
nd=find(G==1521)
CellN=nd;
BinSize=(max(behav.TXVt(:,2))./max(xtvts2(:,1)));
Seg=size(tiral,2)-1;
% close all
figure('position',[500 0 1800 600])

for k=1:Seg;
nBins=max(xtvts2(:,1));

smooth_phase_map=10;
tiral_plot=[];
tiral_plot=tiral{1, k};
subplot(3,Seg,k)
M=[];
M(:,:)=tiral_plot.PF(CellN,:,:);
% imagesc(matnorm(M,2))
imagesc(M)
Linear=tiral_plot.Linear(CellN,:);
title(['SPI=' num2str(round(Linear(2),2)) ' bits/sec'])
xlabel(['Oly1=' num2str(round(Linear(6),3)) '  Oly2=' num2str(round(Linear(7),3))])

subplot(3,Seg,Seg+k)

W=(tiral_plot.Field_edge(CellN,3)-tiral_plot.Field_edge(CellN,1));
rows=tiral_plot.rows(CellN,:); 
a=tiral_plot.Field_edge(CellN,4);
b=tiral_plot.Field_edge(CellN,7);
[~, c]=max(rows);
plot(rows)
hold on
plot([a a],[0 max(rows)],'--')
hold on
plot([b b],[0 max(rows)],'--')
hold on
plot([a b],[max(rows) max(rows)],'r')

hold on
plot([c c],[0 max(rows)],'--')
title(['W=' num2str(round(W))])


subplot(3,Seg,2*Seg+k)
M=[];
M=tiral_plot.In_Spike{1,CellN};
% Slope_range=[5];
% [circ_lin_corr pval slope_deg phi0_deg RR] = cl_corr( M(:,2),radtodeg(M(:,5)),-Slope_range(1),Slope_range(1));
xi=(tiral_plot.intervf(CellN,1));
xe=(tiral_plot.intervf(CellN,2));
B1=xi*BinSize;
B2=xe*BinSize;

plot([M(:,2); M(:,2)],[radtodeg(M(:,5)); radtodeg(M(:,5))+360],'.k','markersize',10)
hold on
plot([B1 B1],[0 720],'r')
hold on
plot([B2 B2],[0 720],'g')

hold on
binrange = [0:2*pi/100:4*pi];
x = binrange;
y1 = -cos(x)/10+B2;
x=radtodeg(binrange);
plot([y1],[x+B2],'r','linewidth',2)
ylim([0 720])

meanO1=radtodeg(tiral_plot.Phase_info(1,1,CellN));
meanO4=radtodeg(tiral_plot.Phase_info(4,1,CellN));
meanOall=radtodeg(tiral_plot.Phase_info(5,1,CellN));
Range=meanO4-meanO1;

title(['\phi1=' num2str(round(meanO1))  '  \phi4=' num2str(round(meanO4))])
xlabel([' \phiall=' num2str(round(meanOall)) '   Range=' num2str(round(Range)) ])
% slope=[];phi0=[];
% t=[0 (max(M(:,2))-min(M(:,2)))]; phi0=phi0_deg; slope=slope_deg;Pval=pval;RR=RR;
% xlabel(['\phi_0 =' num2str(round(phi0)) '  slope =' num2str(round(slope/36)) '  P=' num2str(round(Pval)) '  RR=' num2str(round(RR,2))] )
% 
% hold on
% for ii=1:2
%     far=slope.*t+phi0+360*(ii-1);
%     plot(t+min(M(:,2)),far,'Color','b','linewidth',2)
%     hold on
% end
end

%%

figure('position',[100 -300 1800 200])

meanO_all=[];meanO1=[];meanO4=[];Range=[];Linear=[];W=[];Speed_Spike=[];PeakFR=[];meanFR=[];
for ii=1:(Seg)  
   meanO1=[meanO1; tiral{1,ii}.Phase_info(1,1,CellN)];
   meanO4=[meanO4; tiral{1,ii}.Phase_info(4,1,CellN)]; 
   meanO_all=[meanO_all; tiral{1,ii}.Phase_info(4,1,CellN)]; 
   Linear=[Linear;tiral{1,ii}.Linear(CellN,:)];
   W=[W; (tiral{1,ii}.Field_edge(CellN,3)-tiral{1,ii}.Field_edge(CellN,1))];
   Speed_Spike=[Speed_Spike;nanmean(tiral{1,ii}.SPK_speed)];
   PeakFR=[PeakFR;max(tiral{1,ii}.rows(CellN,:))];
   meanFR=[meanFR;nanmean(tiral{1,ii}.rows(CellN,:))];
end
Range=meanO1-meanO4;

subplot(1,7,1)
plot(radtodeg(meanO1),'.r','markersize',20)
hold on
plot(radtodeg(meanO1),'--','markersize',10)
hold on
binrange = [0:2*pi/100:2*pi];
y1 = -cos(binrange)*2+2;
x=radtodeg(binrange);
plot([y1],[x],'k','linewidth',1)
ylim([0 360])

title('firstQ phase')

subplot(1,7,2)
plot(radtodeg(Range),'.b','markersize',10)
hold on
plot(abs(radtodeg(Range)),'or','linewidth',1)
hold on
plot(abs(radtodeg(Range)),'--','markersize',10)
% hold on
% plot(abs(radtodeg(Range)),'--','markersize',10)
title('Range')

subplot(1,7,3)
plot(Linear(:,2),'.k','markersize',20)
hold on
plot(Linear(:,2),'--','markersize',20)
title('SPI')

subplot(1,7,4)
plot(Linear(:,6),'.k','markersize',20)
hold on
plot(Linear(:,6),'--','markersize',10)
title('Olypher')



subplot(1,7,5)
plot(Linear(:,9),'.k','markersize',20)
hold on
plot(Linear(:,9),'--','markersize',10)
hold on
plot(Speed_Spike,'.r','markersize',20)
hold on
plot(Speed_Spike,'--','markersize',10)
title('Speed')


subplot(1,7,6)
plot(W,'.k','markersize',20)
hold on
plot(W,'--','markersize',10)
title('Width')


subplot(1,7,7)
plot(meanFR,'.k','markersize',20)
hold on
plot(meanFR,'--','markersize',10)
hold on
plot(PeakFR,'.r','markersize',20)
hold on
plot(PeakFR,'--r','markersize',10)

title('FR')
xlabel(['CellN = ' num2str(G(CellN))])


%%
CellN=nd;
BinSize=(max(behav.TXVt(:,2))./max(xtvts2(:,1)));
Seg=size(tiral,2)-1;
% close all
figure('position',[100 -300 300 800])

% for k=1:Seg;
k=Seg+1;

nBins=max(xtvts2(:,1));
smooth_phase_map=10;

tiral_plot=[];
tiral_plot=tiral{1, k};
subplot(6,1,1)
M=[];
M(:,:)=tiral_plot.PF(CellN,:,:);
imagesc(matnorm(M,2))
Linear=tiral_plot.Linear(CellN,:);
title(['SPI=' num2str(round(Linear(2),2)) ' bits/sec', 'Olypher1=' num2str(round(Linear(6),3))])
grid on;
set(gca,'xticklabel',{[]})

subplot(7,1,2)
M=[];
M=tiral_plot.In_Spike{1,CellN};
% Slope_range=[5];
% [circ_lin_corr pval slope_deg phi0_deg RR] = cl_corr( M(:,2),radtodeg(M(:,5)),-Slope_range(1),Slope_range(1));
xi=(tiral_plot.intervf(CellN,1));
xe=(tiral_plot.intervf(CellN,2));
B1=xi*BinSize;
B2=xe*BinSize;

plot([M(:,2); M(:,2)],[radtodeg(M(:,5));radtodeg(M(:,5))+360],'.k','markersize',10)
hold on
plot([B1 B1],[0 720],'r')
hold on
plot([B2 B2],[0 720],'g')

hold on
binrange = [0:2*pi/100:4*pi];
x = binrange;
y1 = -cos(x)/10+B2;
x=radtodeg(binrange);
plot([y1],[x+B2],'r','linewidth',2)
ylim([0 720])

meanO1=radtodeg(tiral_plot.Phase_info(1,1,CellN));
meanO4=radtodeg(tiral_plot.Phase_info(4,1,CellN));
meanOall=radtodeg(tiral_plot.Phase_info(5,1,CellN));
meanO1(meanO1<pi)=meanO1(meanO1<pi)+pi;
Range=meanO1-meanO4;

title(['\phi1=' num2str(round(meanO1))  '  \phi4=' num2str(round(meanO4)),' \phiall=' num2str(round(meanOall)) '   Range=' num2str(round(Range))])
set(gca,'xticklabel',{[]})
% slope=[];phi0=[];
% t=[0 (max(M(:,2))-min(M(:,2)))]; phi0=phi0_deg; slope=slope_deg;Pval=pval;RR=RR;
% xlabel(['\phi_0 =' num2str(round(phi0)) '  slope =' num2str(round(slope/36)) '  P=' num2str(round(Pval)) '  RR=' num2str(round(RR,2))] )
% 
% hold on
% for ii=1:2
%     far=slope.*t+phi0+360*(ii-1);
%     plot(t+min(M(:,2)),far,'Color','b','linewidth',2)
%     hold on
% end

meanO_all=[];meanO1=[];meanO4=[];Range=[];Linear=[];
for ii=1:(Seg)
    
   meanO1=[meanO1; tiral{1,ii}.Phase_info(1,1,CellN)];
   meanO4=[meanO4; tiral{1,ii}.Phase_info(4,1,CellN)]; 
   meanO_all=[meanO_all; tiral{1,ii}.Phase_info(4,1,CellN)]; 
   Linear=[Linear;tiral{1,ii}.Linear(CellN,:)];
   
end
% meanO1(meanO1<pi)= meanO1(meanO1<pi)+2*pi;
% meanO4(meanO4<pi)= meanO4(meanO4<pi)+2*pi;

Range=meanO1-meanO4;

subplot(7,1,3)
plot(radtodeg(meanO1),'.r','markersize',20)
hold on
plot(radtodeg(meanO1),'--','markersize',10)
hold on
% binrange = [pi:2*pi/100:3*pi];
% y1 = sin(binrange)*2+2;
binrange = [0:2*pi/100:2*pi];
y1 = -cos(binrange)*2+2;
x=radtodeg(binrange);
plot([y1],[x],'k','linewidth',1)
% ylim([180 180+360])
ylim([0 360])

title('firstQ phase')


subplot(7,1,4)
plot(radtodeg(meanO1),'.r','markersize',20)
hold on
plot(radtodeg(meanO1),'--','markersize',10)
hold on
plot(radtodeg(meanO4),'.b','markersize',20)
hold on
plot(radtodeg(meanO4),'--','markersize',10)
hold on
plot([y1],[x],'k','linewidth',1)
ylim([0 360])
title('lastQ and firstQ phase')

subplot(7,1,5)
plot(radtodeg(Range),'.b','markersize',10)
hold on
plot(abs(radtodeg(Range)),'or','linewidth',1)
hold on
plot(abs(radtodeg(Range)),'--','markersize',10)
% hold on
% plot(abs(radtodeg(Range)),'--','markersize',10)
title('Range')

subplot(7,1,6)
plot(Linear(:,2),'.k','markersize',20)
hold on
plot(Linear(:,2),'--','markersize',20)
title('SPI')

subplot(7,1,7)
plot(Linear(:,6),'.k','markersize',20)
hold on
plot(Linear(:,6),'--','markersize',10)
title('Olypher')
xlabel(['CellN = ' num2str(G(CellN))])

%%
TR=[1 tr(end)];
figure('position',[500 500 1100 200])

%%###############################################################plot FRmap
subplot(1,5,1)
[matC]=plot_PF_1D_circular(Rate_Matrix,smooth_rate_map,nBins,CellN);
[Linear_Right,Linear_left,Linear_Circ] = linear_maze_PFinfo(sessions{ses},CellN,Rate_Matrix,smooth_rate,TR);
title(['CellN=' num2str(G(CellN)) '  SPI=' num2str(round(Linear_Circ(2),2)) ' bits/sec'])
%%###############################################################plot Speed
M=behav.TXVt;
Speed_int=[];
for i=1:length(int_tr)
    nd=[];
    nd=find(M(:,4)==int_tr(i));
    Speed_int=[Speed_int;nanmean(M(nd,3))*100];
end


Speed_las=[];
for i=1:length(las_tr)
    nd=[];
    nd=find(M(:,4)==las_tr(i));
    Speed_las=[Speed_las;nanmean(M(nd,3))*100];
end

subplot(1,5,2)
plot([ones(length(int_tr),1) 2*ones(length(int_tr),1)],[Speed_int Speed_las(1:4)],'.','markersize',20)
xlim([0 3])
ylim([0 150])
ylabel('Speed')
title('4initial   4last')

%% ###############################################################plot Phase
M=Spike_all;
slope=[];phi0=[];
subplot(1,5,3)
plot([M(:,2); M(:,2)],[radtodeg(M(:,5)); radtodeg(M(:,5))+360],'.k')
t=boundaries_all(1,1):0.01:boundaries_all(1,2); phi0=radtodeg(statsPP_all.intercept); slope=radtodeg(statsPP_all.slope)  ;
Pval=statsPP_all.p;
RR=statsPP_all.r2;

hold on
binrange = [0:2*pi/100:4*pi];
xi = binrange;
y1 = -cos(xi)/40+boundaries_all(2);
xe=radtodeg(binrange);
plot([y1],[xe],'r','linewidth',2)
ylim([0 720])


hold on
for ii=1
    far=slope.*t+phi0+360*(ii-1);
    plot(t,far,'Color','b','linewidth',1.5)
    hold on
end
ylim([0 720])
title(['\phi_0 =' num2str(round(phi0)) '  slope =' num2str(round(slope)) '  P=' num2str(round(Pval)) '  RR=' num2str(round(RR,2))] )
clc
%%
subplot(1,5,4)
[meanR,meanO]=meanphase(M(:,5));

n=14;
h=polarhistogram(M(:,5),n,'linewidth',3,'facecolor',[.6 .6 .6],'edgecolor','k','FaceAlpha',.3,'Normalization','probability')
title(['\Delta\phiall =' num2str(round(radtodeg(meanO))) ])

subplot(1,5,5)

Mphase=Phase_Matrix(:,4+CellN);
matC=reshape(Mphase,nBins,length(Mphase)/nBins)';
% imagesc(matC)
M=matC;
M(isnan(M)==1)=0;
ngrid=nBins;
Fullwin=(-ngrid:ngrid)';
Smoother=exp(-Fullwin.^2/(smooth_phase_map/2)^2);
S=conv2(Smoother, Smoother, M, 'same');
imagesc(S)
% colorbar

[track_info_ph,pos_info_val_ph] = bz_olypherInfo(matC,2,smooth_phase,'phase');
Olypher_1=max(track_info_ph);Olypher_2=max(max(pos_info_val_ph));

title(['Olypher1=' num2str(round(Olypher_1,3)) '  Olypher2=' num2str(round(Olypher_2,3))])



%%
% [curve,stats] = FiringCurve(TX,Spike(:,1),'nBins',nBins,'threshold',0.1,'minSize',14,'minPeak',2,'smooth',3)


% M=tiral.In_Spike_Cric{1,CellN};
% Spike=M;
% Spike(:,2)=(Spike(:,2)-min(Spike(:,2)))/(max(Spike(:,2))-min(Spike(:,2)));
% boundaries(1,1)= min(Spike(:,2));
% boundaries(1,2)= max(Spike(:,2));
% [dataPP,statsPP] = PhasePrecession_F(TX,Spike(:,1),thetaP, 'boundaries',boundaries)
% phi0_deg2=radtodeg(statsPP.intercept); slope_deg2=statsPP.slope ; Pval2=statsPP.p;RR2=statsPP.r2;
% ylabel(['\phi_0 =' num2str(round(phi0_deg2)) '  slope =' num2str(round(slope_deg2*180)/(nBins*pi),2) '  P=' num2str(round(Pval2)) '  RR=' num2str(round(RR2,2))] )



% subplot(1,4,2)
% plot(tiral.rows_Cric(CellN,:,:))
% hold on
% plot([xi-nBins xi-nBins],[0 max(tiral.rows_Cric(CellN,:,:))])
% hold on
% plot([xe-nBins xe-nBins],[0 max(tiral.rows_Cric(CellN,:,:))])
% title([ num2str(Q{h}) '         '])



% subplot(1,4,4)
% Mphase=Phase_Matrix(:,CellN+4);
% Phase_Map=[];
% Phase_Map=reshape(Mphase,nBins,length(Mphase)/nBins)';
% matC=Phase_Map(TRL(1):TRL(2),:);
% M=matC;
% M(isnan(M)==1)=0;
% ngrid=nBins;
% Fullwin=(-ngrid:ngrid)';
% Smoother=exp(-Fullwin.^2/(smooth_phase_map/2)^2);
% S=conv2(Smoother, Smoother, M, 'same');
% imagesc(S)




