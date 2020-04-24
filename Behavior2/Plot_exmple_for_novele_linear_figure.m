%##############################################################################################################
%Plot_exmple_for_novele_linear_figure
%%
trN=5;speed_thr=0.04;
delta=0.3;smooth_rate=10;smooth_phase=10;
Rate_Matrix=xtvts2;
Phase_Matrix=xtvtph;
BinSize=(max(behav.TXVt(:,2))./max(xtvts2(:,1)));
smooth=10;
TR=[];
TR=unique(behav.TXVt(:,4));
TR(TR==0)=[];
r0=[];TRi=[];TR_all=[];TR_even=[];
TR_even=find(mod(TR,2)==0);
TR_all=TR(TR_even);Steps=2;% 1 D linear
r0=[1 trN];
Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
TRi=TR_all(1:Seg*trN);
TRi=reshape(TRi,[trN Seg])'

%%
Matrix=xtvts2;
% tr=TRi(6,:)
% tr=TRi(8,:)
tr=[2:2:20]
trials_port=tr;
nd_int=[];
for i=1:length(trials_port)
    nd_int=[nd_int;find(Matrix(:,4)==trials_port(i))];
end
Matrix=Matrix(nd_int,:);
nPFbin=max(Matrix(:,1));
nd=find(G==1521)
CellN=nd;

%%

M=[];
M(:,1:4)=Matrix(:,1:4);
M(:,5)=Matrix(:,CellN+4);
M=[M; M; M];

size(M)
smoothT=smooth1D(M(:,2),smooth,1);
smoothC=smooth1D(M(:,5),smooth,1);
ncell=length(smoothC(1,:));
rate=smoothC./repmat(smoothT,1,ncell);
xtsr=[M(:,[1 3 4]) rate];
rr=xtsr(:,4);
%matC=reshape(rr,nPFbin,length(rr)/nPFbin)';
matC=reshape(rr,3*nPFbin,length(tr))';
MatC=matC;
figure
%plot(mean(MatC))
imagesc(MatC)

%%

trN=5;speed_thr=0.04;
delta=0.3;smooth_rate=10;smooth_phase=10;
CellN=nd;
[Field_edge,rows,PF] = Feild_width_linearMaze(CellN,smooth_rate,delta,Rate_Matrix,tr);
[Linear] = linear_maze_PFinfo(CellN,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,tr)

%%

figure('position',[500 200 150 300])
M=[];
M=[rows rows rows];
W=(Field_edge(1,3)-Field_edge(1,1));
a=Field_edge(1,1);
b=Field_edge(1,3);
% [~, c]=max(rows);
c=Field_edge(1,2);
plot(mean(MatC))
hold on
plot([a a],[0 max(rows)],'--')
hold on
plot([b b],[0 max(rows)],'--')
hold on
plot([a b],[max(rows) max(rows)],'r')
hold on
plot([c c],[0 max(rows)],'--')
title(['W=' num2str(round(W))])
xlim([108 186])
% set(gca,'fontsize',20)

%% extract Spikes
folder='Cellinfo3';
cd (folder)
load(['TXVtPh_' num2str(CellN) '.mat'])
cd ..

TXVtPh(isnan(TXVtPh(:,2))==1,:)=[];
TXVtPh(find(TXVtPh(:,3)<speed_thr),:)=[];%.04
TXVtPh(TXVtPh(:,4)==(0),:)=[];

trials_port=tr;
nd_int=[];
for i=1:length(trials_port)
    nd_int=[nd_int;find(TXVtPh(:,4)==trials_port(i))];
end
TXVtPh=TXVtPh(nd_int,:);
interv=[Field_edge(1,1) Field_edge(1,3)-1];
[Phase_info_5,In_Spike]=PF_sectionphase(TXVtPh,behav.TXVt,interv,BinSize,'all');

figure('position',[100 200 150 300])
M=[];
M=In_Spike;

Slope_range=[5];
[circ_lin_corr pval slope_deg phi0_deg RR] = cl_corr( M(:,2),radtodeg(M(:,5)),-Slope_range(1),Slope_range(1))
plot([M(:,2); M(:,2)],[radtodeg(M(:,5));radtodeg(M(:,5))+360],'.k','markersize',10)

slope=[];phi0=[];
t=[0 (max(M(:,2))-min(M(:,2)))]; phi0=phi0_deg; slope=slope_deg;Pval=pval;RR=RR;
xlabel(['\phi_0 =' num2str(round(phi0)) '  slope =' num2str(round(slope/36)) '  P=' num2str(round(Pval)) '  RR=' num2str(round(RR,2))] )

hold on
for ii=2
    far=slope.*t+phi0+360*(ii-1);
    plot(t+min(M(:,2)),far,'Color','b','linewidth',2)
    hold on
end
hold on
binrange = [0:2*pi/100:4*pi];
B2=3;
x = binrange;
y1 = -cos(x)/10+B2;
x=radtodeg(binrange);
plot([y1],[x+B2],'r','linewidth',2)

hold on
plot([2.6 3.2],[radtodeg(Phase_info_5(1)) radtodeg(Phase_info_5(1))])
xlim([2.6 3.2])
yticks([0 180 360 540 720])
ylim([-50 720])
% ylim([-50 360])
% set(gca,'fontsize',20)
%%
figure
[meanR,meanO]=meanphase(M(:,5));
n=32;
% h=polarhistogram(M(:,5),n,'linewidth',3,'facecolor',[.6 .6 .6],'edgecolor','k','FaceAlpha',.3,'Normalization','probability')
h=polarhistogram(TXVtPh(:,5),n,'linewidth',3,'facecolor',[.6 .6 .6],'edgecolor','k','FaceAlpha',.3,'Normalization','probability')

%  PlotCircularDistribution(dist,a,s,varargin)


%%
clc
% [dist,binned,stats] = CircularDistribution(M(:,5),'nBins',32,'smooth',2);
[dist,binned,stats] = CircularDistribution([TXVtPh(:,5)],'nBins',32,'smooth',2);
figure
plot([binned binned+2*pi],[dist dist])



