
clearvars;clc;
dirData = 'Y:\OptoMECLEC\';
animal  = {'OML18','OML18','OML18','OML19','OML19'};
sessions  = {'day1','day2','day4','day2','day3'};
ses=3;
dirses = [dirData animal{ses} '\'  sessions{ses}];
cd(dirses);

%%
filename=[animal{ses} sessions{ses}];
load([sessions{ses} '.spikes.cellinfo.mat'])
load([filename '.mat' ])
%%
trN=5;speed_thr=0.03;
delta=0.3;smooth_rate=10;smooth_phase=10;
Rate_Matrix=xtvtls2;
Phase_Matrix=xtvtlph;
nBins=(max(xtvtls2(:,1)))
smooth=10;
TR=[];

%%
mkdir('PFs')
TR=unique(behav.txvtl(:,4));
TR(TR==0)=[];
r0=[];TRi=[];TR_all=[];TR_even=[];
TR_even=find(mod(TR,2)==0);
TR_On=TR(TR_even);

TR=[];
TR=unique(behav.txvtl(:,4));
TR(TR==0)=[];
r0=[];TRi=[];TR_all=[];TR_Odd=[];
TR_Off=find(mod(TR,2)==1);
TR_Off=TR(TR_Off);

tic
for h = 10%:size(spikes.times,2)
    
    
    cellN = h
    TXVtlPh=[];
    TXVtlPh=spikes.maze{1,cellN};
    TXVtlPh(isnan(TXVtlPh(:,2))==1,:)=[];
    TXVtlPh(find(TXVtlPh(:,3)<speed_thr),:)=[];%.04
    TXVtlPh(TXVtlPh(:,4)==(0),:)=[];
    
    
    tr=TR_On;
    trials_port=tr;
    nd_int=[];
    for i=1:length(trials_port)
        nd_int=[nd_int;find(TXVtlPh(:,4)==trials_port(i))];
    end
    TXVtlPh_ON=TXVtlPh(nd_int,:)  ;
    curve=[];stats=[];dataP=[];statsP=[];
    [curve,stats] = FiringCurve(behav.txvtl(:,1:2),TXVtlPh_ON(:,1),'nBins',nBins,'threshold',0.1,'minSize',10,'minPeak',1,'smooth',3);
    Bounds=[];
    if isnan(stats.fieldX)==1
        Bounds=[0 1];
    else  
        Bounds=curve.x(stats.fieldX(1,:));
    end
    [dataP,statsP] = PhasePrecession(behav.txvtl(:,1:2),TXVtlPh_ON(:,1),[Theta.Phase.timestamps    Theta.Phase.phase  ],'boundaries',Bounds);
    if isnan(statsP.slope)~=1
    figure('position',[100 100 200 800])
    PlotPhasePrecession(dataP,statsP)
    title(['ON- ' num2str(h)])
    cd PFs
    print('-djpeg',['Phasecell_ON' num2str(h)]) 
    cd ..
    end
    
    curve_ON{1,h}=curve;stats_ON{1,h}=stats;
    dataP_ON{1,h}=dataP;statsP_ON{1,h}=statsP;
    
    cellN = h
    tr=TR_Off;
    trials_port=tr;
    nd_int=[];
    for i=1:length(trials_port)
        nd_int=[nd_int;find(TXVtlPh(:,4)==trials_port(i))];
    end
    TXVtlPh_Off=TXVtlPh(nd_int,:)  ;
    curve=[];stats=[];dataP=[];statsP=[];
    [curve,stats] = FiringCurve(behav.txvtl(:,1:2),TXVtlPh_Off(:,1),'nBins',nBins,'threshold',0.1,'minSize',10,'minPeak',1,'smooth',3);
    Bounds=[];
    if isnan(stats.fieldX)==1
        Bounds=[0 1];
    else  
        Bounds=curve.x(stats.fieldX(1,:));
    end
    [dataP,statsP] = PhasePrecession(behav.txvtl(:,1:2),TXVtlPh_Off(:,1),[Theta.Phase.timestamps    Theta.Phase.phase  ],'boundaries',Bounds);
    if isnan(statsP.slope)~=1
    figure('position',[100 100 200 800])
    PlotPhasePrecession(dataP,statsP)
    title(['Off- ' num2str(h)])
    curve_Off{1,h}=curve;stats_Off{1,h}=stats;
    dataP_Off{1,h}=dataP;statsP_Off{1,h}=statsP;
    cd PFs
    print('-djpeg',['Phasecell_Off' num2str(h)]) 
    cd ..
    end

end
toc

%%
mkdir('Metrics')
cd Metrics
save(['curve_Off'],'curve_Off')
save(['stats_Off'],'stats_Off')
save(['curve_ON'],'curve_ON')
save(['stats_ON'],'stats_ON')

save(['dataP_Off'],'dataP_Off')
save(['statsP_Off'],'statsP_Off')
save(['dataP_ON'],'dataP_ON')
save(['statsP_ON'],'statsP_ON')
cd ..
%%
figure
plot(curve.rate)
%%
% figure
% plot([TXVtlPh_ON(:,2);TXVtlPh_ON(:,2)],[TXVtlPh_ON(:,6)+2*pi;TXVtlPh_ON(:,6)+4*pi],'.k','markersize',5)
% figure
% plot([TXVtlPh_Off(:,2);TXVtlPh_ON(:,2)],[TXVtlPh_ON(:,6)+2*pi;TXVtlPh_ON(:,6)+4*pi],'.k','markersize',5)
% hold on 
% statsP_ON{1, 1}.slope  
% statsP_ON{1, 1}.intercept  
[dataP,statsP] = PhasePrecession(behav.txvtl(:,1:2),TXVtlPh(:,1),[Theta.Phase.timestamps  Theta.Phase.phase],'boundaries',curve.x(stats.fieldX));

%%
plot([TXVtlPh_ON(:,2);TXVtlPh_ON(:,2)],[TXVtlPh_ON(:,6)+2*pi;TXVtlPh_ON(:,6)+4*pi],'.k','markersize',5)
%%
close all
for j=1:29
Onset_On(j)=radtodeg(statsP_ON{1, j}.intercept);
Onset_Off(j)=radtodeg(statsP_Off{1, j}.intercept);
end
figure
plot(Onset_On,'or')
hold on
plot(Onset_Off,'ob')

%%
figure
PlotPhasePrecession(dataP,statsP_ON)






