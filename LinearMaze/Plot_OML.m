clc
clear

ses=1;
dirData = 'Y:\OptoMECLEC\';
animal  = {'OML18','OML18','OML18','OML18','OML19','OML19'};
sessions  = {'day1','day2','day4','day5','day2','day3'};
folder='Metrics';

dirses = [dirData animal{ses} '\'  sessions{ses}];
cd(dirses);
filename=[animal{ses} sessions{ses}]
load([filename '.mat'])

%%
cd 'Metrics'
load('tiral_2')
% load('tiral_1')
% load('PhasePrec.mat')

LM=nanmean(tiral{1, size(tiral,2)}.In_Spike{1, 5}(:,5))
cd ..
% if LM==0
%     statsPP=statsPPoff;
% else
%     statsPP=statsPPon;
% end

behavMatrix=behav.txvtl;
Ratemap=map.xtvtls2;
BinSize=(max(behavMatrix(:,2))./max(Ratemap(:,1)));
Seg=size(tiral,2)-1;
LM=1;

%%

close all
for CellN=1:spikes.numcells
    CellN
    figure('position',[1400 100 300 800])
    
    % for k=1:Seg;
    k=Seg+1;
    nBins=max(Ratemap(:,1));
    smooth_phase_map=10;
    
%     if isstruct(statsPP{1, CellN})==0
%         Slope=nan;
%         intercept=nan;
%     else
%         Slope=statsPP{1, CellN}.slope;
%         intercept=statsPP{1, CellN}.intercept;
%     end
%     Slope
%     intercept
    
    
    tiral_plot=[];
    tiral_plot=tiral{1, k};
    subplot(4,1,1)
    % CellN=20;
    Rf_ch_shank=Ripple.RefChannelAvg(:,1); 
    clc
    channel_order=[];
    channel_order{1}=[12 11 10 9 8 7 6 5 4 3 2 1];
    channel_order{2}=[16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1];
    CellPosition=[];
    Cellchannel=[spikes.Cellinfo(CellN).shankID spikes.Cellinfo(CellN).peakChOrganized];
    
    if Cellchannel(1,1)==3
        nd=find(spikes.channelXML{Cellchannel(1,1),2}==spikes.Cellinfo(CellN).peakChOrganized);
        CellPosition{2}=[spikes.Cellinfo(CellN).shankID channel_order{2}(nd)];
        CellPosition{1}=[nan nan];
    elseif Cellchannel(1,1)>3
        nd=find(spikes.channelXML{Cellchannel(1,1),2}==spikes.Cellinfo(CellN).peakChOrganized);
        CellPosition{1}=[spikes.Cellinfo(CellN).shankID channel_order{1}(nd)];
        CellPosition{2}=[nan nan];
    else
        nd=find(spikes.channelXML{Cellchannel(1,1),2}==spikes.Cellinfo(CellN).peakChOrganized);
        CellPosition{1}=[spikes.Cellinfo(CellN).shankID channel_order{1}(nd)];
        CellPosition{2}=[nan nan];
    end
    
    probeID=4;
    [CellShankN,CellChannelN,RefCellChannelN,RefCellShankN]=Clayout_AntonioProb(CellPosition,Rf_ch_shank,probeID ,channel_order);
    if RefCellChannelN-CellChannelN >= 0 ;  a='Superficial';else a='Deep'; end ;
%     title([a, ' Sh', num2str(CellShankN) ,' Ch' num2str(CellChannelN), '  Slope=' num2str(Slope) ' \phi0=' num2str(intercept)])
        title([a, ' Sh', num2str(CellShankN) ,' Ch' num2str(CellChannelN)])

    
    subplot(5,1,2)
    M=[];
    M(:,:)=tiral_plot.PF(CellN,:,:);
    ngrid=100;smooth=3;
    Fullwin=(-ngrid:ngrid)';
    Smoother=exp(-Fullwin.^2/(smooth/2)^2);
    M = conv2(Smoother, Smoother, M, 'same');
    imagesc(matnorm(M,2))
    % imagesc(M)
    axis xy
    hold on
    plot((mean(M,1)./max(mean(M,1))).*size(M,1),'w','linewidth',2)
    Linear=tiral_plot.Linear(CellN,:);
    title(['SPI=' num2str(round(Linear(2),2)) 'bits/sec ', ' Olyph=' num2str(round(Linear(6),3))])
    grid on;
    set(gca,'xticklabel',{[]})
    %
    % close all
    
    
    subplot(8,1,4)
    M=[];
    M=tiral_plot.In_Spike{1,CellN};
    Phasecolomn=size(tiral_plot.In_Spike{1,CellN},2);
    LaserMode=nanmean(M(:,5));
    % Slope_range=[5];
    % [circ_lin_corr pval slope_deg phi0_deg RR] = cl_corr( M(:,2),radtodeg(M(:,5)),-Slope_range(1),Slope_range(1));
    xi=(tiral_plot.intervf(CellN,1));
    xe=(tiral_plot.intervf(CellN,2));
    B1=xi*BinSize;
    B2=xe*BinSize;
    
    plot([M(:,2);M(:,2)],[radtodeg(M(:,Phasecolomn));radtodeg(M(:,Phasecolomn))+360],'.k','markersize',2)
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
    % meanO1(meanO1<pi)=meanO1(meanO1<pi)+pi;
    Range_all=meanO1-meanO4;
    
    
    
    title(['\phi1=' num2str(round(meanO1))  '  \phi4=' num2str(round(meanO4)),' \phiall=' num2str(round(meanOall)) '   Range=' num2str(round(Range_all))])
    set(gca,'xticklabel',{[]})
    % slope=[];phi0=[];
    % t=[0 (max(M(:,2))-min(M(:,2)))]; phi0=phi0_deg; slope=slope_deg;Pval=pval;RR=RR;
    % xlabel(['\phi_0 =' num2str(round(phi0)) '  slope =' num2str(round(slope/36)) '  P=' num2str(round(Pval)) '  RR=' num2str(round(RR,2))])
    
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
    
    subplot(8,1,5)
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
    set(gca,'xticklabel',{[]})
    
    subplot(8,1,6)
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
    set(gca,'xticklabel',{[]})
    
    subplot(12,1,10)
    plot(radtodeg(Range),'.b','markersize',10)
    hold on
    plot(abs(radtodeg(Range)),'or','linewidth',1)
    hold on
    plot(abs(radtodeg(Range)),'--','markersize',10)
    % hold on
    % plot(abs(radtodeg(Range)),'--','markersize',10)
    title('Range')
    set(gca,'xticklabel',{[]})
    
    subplot(12,1,11)
    plot(Linear(:,2),'.k','markersize',20)
    hold on
    plot(Linear(:,2),'--','markersize',20)
    title('SPI')
    set(gca,'xticklabel',{[]})
    
    subplot(12,1,12)
    plot(Linear(:,6),'.k','markersize',20)
    hold on
    plot(Linear(:,6),'--','markersize',10)
    title('Olypher')
    xlabel(['CellN = ' num2str(CellN) '   Laser ' num2str(LaserMode)],  'fontsize',15)
    
    cd PFs
    if LM==1
        print('-djpeg',['ON_cell_' num2str(CellN)])
    else
        print('-djpeg',['Off_cell_' num2str(CellN)])
    end
    cd ..
end
close all
%%  Making ppt file for each sesion
% figure('position',[1400 100 300 800])
clc
FigureY=800;
FigureX=300;
CellNumber=spikes.numcells;
FigureName='\ON_cell_'; pptName='_LayerON';
% FigureName='\Off_cell_';pptName='_LayerOff';



Fullpath=[dirData animal{ses} '\'  sessions{ses} '\PFs' ]
% Adjust Number of slides
N_image_slide=4;
a=FigureX;% sapce between the figures
b=-30;%raw onset

N_slide=ceil(CellNumber./N_image_slide);
Shape_position=[0*a,b,FigureX,FigureY;a*1,b,FigureX,FigureY;a*2,b,FigureX,FigureY;a*3,b,FigureX,FigureY]

h = actxserver('PowerPoint.Application')
h.Visible = 1;
Presentation = h.Presentation.Add
blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);

for k=1:N_slide
    k
    Slide(k) = Presentation.Slides.AddSlide(1,blankSlide);
    for i=1:N_image_slide
        A=i+(k-1)*N_image_slide;
        if A<CellNumber+1
            Image = Slide(k).Shapes.AddPicture([Fullpath FigureName num2str(A) '.jpg'], 'msoFalse', 'msoTrue', Shape_position(i,1),Shape_position(i,2),Shape_position(i,3),Shape_position(i,4))
        end
    end
end

Presentation.SaveAs([Fullpath '\' filename pptName '.ppt'])
invoke(Presentation,'Close');

h.Quit;
h.delete;

%%

CellG=[];

%OML18day1
ses=1;
cellgroupOff=[];
cellgroupON=[];
cellgroupOff.deep=[29 28 22 13 27 15]
cellgroupOff.sup=[nan]
cellgroupOff.DG=[16 21]

cellgroupON.deep=[28 15]
cellgroupON.sup=[nan]
cellgroupON.DG=[16]
CellG{1,ses}=cellgroupOff;
CellG{2,ses}=cellgroupON;

%OML18day2
ses=2;
cellgroupOff=[];
cellgroupON=[];

cellgroupOff.deep=[30 23 12 9 8 7 1]
cellgroupOff.sup=[36 28 24 19 14 13 11 10 5]
cellgroupOff.DG=[26 25 16 6 4]

cellgroupON.deep=[33 30 27 23 22 21 17 12 8 7 2 1 ]
cellgroupON.sup=[36 31 28 24 20  19 18 14 13 11 5]
cellgroupON.DG=[32 29 6]

CellG{1,ses}=cellgroupOff;
CellG{2,ses}=cellgroupON;
%OML18day4
ses=3;
cellgroupOff=[];
cellgroupON=[];
cellgroupOff.deep=[36 26 28 21 5 7 3]
cellgroupOff.sup=[33 31 32 27 23 17  20 13  14 1]
cellgroupOff.DG=[30 22 18 6 2]



cellgroupON.deep=[38 34 36 29 28 21 23 5 3]
cellgroupON.sup=[33 31 32 25 27 17 20 14 12 1]
cellgroupON.DG=[30 22 24 18 6 2]

CellG{1,ses}=cellgroupOff;
CellG{2,ses}=cellgroupON;
%OML19day2
ses=4;
cellgroupOff=[];
cellgroupON=[];
cellgroupOff.deep=[29 5 7]
cellgroupOff.sup=[28 21 22 24 19 20 14 15 10 11 6 3 4]
cellgroupOff.DG=[16 12]



cellgroupON.deep=[25 26 7]
cellgroupON.sup=[27 28 21 22 24 19 20 13 14 15 9 10 11 6 8 1 3 4]
cellgroupON.DG=[29 16 12 2]

CellG{1,ses}=cellgroupOff;
CellG{2,ses}=cellgroupON;
%OML19day3
ses=5;
cellgroupOff=[];
cellgroupON=[];
cellgroupOff.deep=[19 14 2]
cellgroupOff.sup=[22 23 17 20 13 15 16 9 11 12 5 7 3 4]
cellgroupOff.DG=[21 18  10 6 8 1]



cellgroupON.deep=[19 2]
cellgroupON.sup=[22 23 17 20 12 5 3]
cellgroupON.DG=[ 21 18 10 6 1]

CellG{1,ses}=cellgroupOff;
CellG{2,ses}=cellgroupON;


































