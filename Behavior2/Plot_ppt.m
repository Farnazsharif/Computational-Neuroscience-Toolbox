%% for linear maze

close all
load([filename '.mat'])
Matrix=xtvts;

smooth=10;
cellN=1;
nPFbin=50;
for i = 1:length(G)
    cellN = i;
    figure('position',[500 500 200 280])
    plot_PF_1D(Matrix,smooth,nPFbin,cellN);
%     plot_PF_1D_circular(Matrix,smooth,nPFbin,cellN);
     
    if ismember(G(i),Cellgroup.Pyr)==1
    title(['Cell# = ' num2str(G(i)) ' Pyr'])
    else
    title(['Cell# = ' num2str(G(i)) ' Int'])  
    end
    cd PFs
    print('-djpeg',['cell_' num2str(i)]) 
    cd ..
    close all
end


%%  PF info
close all
load([filename '.mat'])

cd 'Metrics'
load('tiral_even')
cd ..

% smooth=10;
% cellN=1;
% nPFbin=50;
% for i = 1:length(G)
%     cellN = i;
%     figure('position',[500 500 200 280])
% %     plot_PF_1D(Matrix,smooth,nPFbin,cellN);
%     plot_PF_1D_circular(Matrix,smooth,nPFbin,cellN);
%      
%     if ismember(G(i),Cellgroup.Pyr)==1
%     title(['Cell# = ' num2str(G(i)) ' Pyr'])
%     else
%     title(['Cell# = ' num2str(G(i)) ' Int'])  
%     end
%     cd PFs
%     print('-djpeg',['cell_' num2str(i)]) 
%     cd ..
%     close all
% end

mkdir('PF_trials')

for i = 1:length(G)
    CellN = i;
figure('position',[500 0 300 800])


k=size(tiral,2);

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


    if ismember(G(i),Cellgroup.Pyr)==1
    xlabel(['Cell# = ' num2str(G(i)) ' Pyr'])
    else
    xlabel(['Cell# = ' num2str(G(i)) ' Int'])  
    end
    cd 'PF_trials'
    print('-djpeg',['cell_' num2str(i)]) 
    cd ..
    close all
end

%%  Making ppt file for each sesion

Fullpath='Y:\Data\GrosmarkAD\Achilles\Achilles_10252013\PF_trials';
% Adjust Number of slides
N_image_slide=4;
a=300;% sapce between the figures
b1=-5;%first raw onset
b2=350;%second raw onset
N_slide=ceil(length(G)./N_image_slide);
Shape_position_1=[0*a,b1,300,800;a*1,b1,300,800;a*2,b1,300,800;a*3,b1,300,800]
Shape_position=[Shape_position_1];

h = actxserver('PowerPoint.Application')
h.Visible = 1;
Presentation = h.Presentation.Add
blankSlide = Presentation.SlideMaster.CustomLayouts.Item(7);

for k=1:N_slide
    Slide(k) = Presentation.Slides.AddSlide(1,blankSlide);
for i=1:N_image_slide
    A=i+(k-1)*N_image_slide;
    if A<length(G)+1
    Image = Slide(k).Shapes.AddPicture([Fullpath '\cell_' num2str(A) '.jpg'], 'msoFalse', 'msoTrue', Shape_position(i,1),Shape_position(i,2),Shape_position(i,3),Shape_position(i,4))   
    end
end
end

% for k=1:N_slide
%     Slide(k) = Presentation.Slides.AddSlide(1,blankSlide);
% for i=1:10
%     A=i+(k-1)*10;
%     if A<length(G)+1
%     Image = Slide(k).Shapes.AddPicture([Fullpath '\cell_' num2str(A) '.jpg'], 'msoFalse', 'msoTrue', Shape_position(i,1),Shape_position(i,2),Shape_position(i,3),Shape_position(i,4))   
%     end
% end
% end


Presentation.SaveAs([Fullpath '\' filename '_even' '.ppt'])

% Presentation.SaveAs([Fullpath '\' filename '_odd' '.ppt'])
invoke(Presentation,'Close');

h.Quit;
h.delete;


