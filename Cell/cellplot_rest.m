%cellplot
function cellplot_rest(shankID,CellN,Cellg)
% 
%  shankID=15;
%  CellN=5;
%  Cellg=1505;

filename1='FM05_1';
filename2='FM05_1_clu_';
% filename3='FM05_1PlaceField.mat';
filename3='FM05_1PlaceFieldRest.mat';
filename4=('FM05_1PlaceFieldTotal.mat');
load('PrePostG1.mat')
load('PrePostG4.mat')
% load('SpatialC.mat')
load('SpatialCrest.mat')
load('thetaf_rest.mat')
load('thetaI_rest.mat')
load('thetaf_run.mat')
load('thetaI_run.mat')
load('refractoryT_rest.mat')
load('refractoryT_run.mat')
load('burstIndex_run.mat')
load('burstIndex_rest.mat')
load(filename1)
load(filename3)
load(filename4)
load([filename2,num2str(shankID),'.mat'])
ndex=find(clu.g==CellN);
Cndex=find(G==Cellg);

 %% Plot firing field

plotPFSetObj2(filename1,Cndex,xttscR)


plotPFSetObj_whole(filename1,Cndex,xttscBo)


% title(['shank' num2str(shankID) 'cell' num2str(CellN) ])

%%  Plot mean fields

axes('Position',[.47 .5 .12 .42])
ndxzero=find(xttscR(:,1)==1);
sets=unique(xttscR(:,4));
binN=100;
for i=1:length(sets)
s=find(xttscR(ndxzero,4)==sets(i));
Set(i,1)=s(1);
Set(i,2)=s(end);
end

cc=reshape(xttscR(:,4+Cndex),binN,142);
scale=10;
I=1:2:200;
I2=imresize(I,[1 length(I)*scale],'lanczos3');

for i =1 :length(sets)   
se=cc(:,Set(i,1):Set(i,2));
S(i,:)=imresize(sum(se,2)',[1 binN*scale],'lanczos3');
NS(i,:)=matnorm(imresize(sum(se,2)',[1 binN*scale],'lanczos3'),2);
IS(i,:)=matnorm(sum(se,2),1);
end

plot (I2,NS(3,:),'b')
hold on
plot(I2,NS(2,:)+1,'r')
hold on
plot(I2,NS(1,:)+2,'k')

set(gca,'YTickLabel',[])
set(gca,'YTickLabel',[3 2 1])
set(gca,'Ytick',[0.6 1.5 2.5])
% ylabel('Nornalized spk # per set','fontsize',11)

axes('Position',[0.64, 0.5, 0.12, 0.42])
box on
plot (I2,S(3,:),'b')
hold on
plot(I2,S(2,:)+50,'r')
hold on
plot(I2,S(1,:)+100,'k')
% title('Sum')
% ylabel('spik # ','fontsize',11)

axes('Position',[.47 .93 .12 .07])
box on
imagesc(IS)
set(gca,'XTickLabel',[])

set(gca,'YTickLabel',[])
%% spatial coverage and firin rate
[~,prendx]=ismember(PrePostG1(:,1),G);
[~,postndx]=ismember(PrePostG1(:,2),G);

axes('Position',[.8, 0.5, 0.195, 0.42])

sety=2;
plot(SpatialCrest(:,sety+1), Frate.rest,'x','markersize',5,'LineWidth',1) 
hold on
plot(SpatialCrest(postndx,sety+1), Frate.rest(postndx),'kd','markersize',5,'LineWidth',2) 
hold on
plot(SpatialCrest(prendx,sety+1), Frate.rest(prendx),'go','markersize',5,'LineWidth',2) 
% xlabel('SpatialC')
% ylabel('Frate.rest')
hold on
sety=1;
plot(SpatialCrest(Cndex,sety+1), Frate.rest(Cndex),'ro','markersize',9,'LineWidth',3) 
sety=2;
hold on
plot(SpatialCrest(Cndex,sety+1), Frate.rest(Cndex),'rd','markersize',9,'LineWidth',3) 
sety=3;
hold on
plot(SpatialCrest(Cndex,sety+1), Frate.rest(Cndex),'rs','markersize',9,'LineWidth',3) 
% legend('C','post','pre','S1','S2','S3','Location','northwest')
xlim([0 .3])

%% plot celltypes_rest 
axes('Position',[0.34, 0.03, 0.3, 0.4])


plot3(spkinfo.duration,refractoryT_rest,burstIndex_rest,'*','markersize',5,'LineWidth',1) 
grid on
hold on
plot3(spkinfo.duration(prendx),refractoryT_rest(prendx),burstIndex_rest(prendx),'go','markersize',8,'LineWidth',2) 
hold on
plot3(spkinfo.duration(postndx),refractoryT_rest(postndx),burstIndex_rest(postndx),'kd','markersize',8,'LineWidth',2) 
hold on
xlabel('S.D')
ylabel('RT')
% title('rest')
% zlabel('B.I')
% set(gcf,'Color','w');
hold on
plot3(spkinfo.duration(Cndex),refractoryT_rest(Cndex),burstIndex_rest(Cndex),'rs','markersize',14,'LineWidth',3)
% plot celltypes_run 
axes('Position',[.03 .03 .3 .4])
plot3(spkinfo.duration,refractoryT_run,burstIndex_run,'*','markersize',5,'LineWidth',1) 
grid on
hold on
plot3(spkinfo.duration(prendx),refractoryT_run(prendx),burstIndex_run(prendx),'go','markersize',8,'LineWidth',2) 
hold on
plot3(spkinfo.duration(postndx),refractoryT_run(postndx),burstIndex_run(postndx),'kd','markersize',8,'LineWidth',2) 
hold on
xlabel('S.D')
ylabel('RT')
% zlabel('B.I')
% set(gcf,'Color','w');
hold on
plot3(spkinfo.duration(Cndex),refractoryT_run(Cndex),burstIndex_run(Cndex),'rs','markersize',14,'LineWidth',3)
% title('run')

dim4 = [.03 .34 .27 .07];
str4=('run');
annotation('textbox',dim4,'String',str4,'FitBoxToText','on','LineStyle','none');
%
dim5 = [.32 .34 .27 .07];
str5=('rest');
annotation('textbox',dim5,'String',str5,'FitBoxToText','on','LineStyle','none');

%% Annotations

sp=0.11;
u=[0 sp*1 sp*2 sp*3 sp*4 sp*5 sp*6];
 [a,b]=find(Cellg==PrePostG4) ;  
if b==1
 str1=('Cell=Pre,');
 annotation('textbox',[.44 .94 .15 .07],'String',str1,'FitBoxToText','on','LineStyle','none');
 
    for ii=1:length(a)
    str{ii}=([ 'post' num2str(ii) '=' num2str(PrePostG4(a(ii),2)) ',']);
    annotation('textbox',[.53+u(ii) .94 .15 .07],'String',str{ii},'FitBoxToText','on','LineStyle','none');
    hold on
    end
elseif b==2
 str1=('Cell=Post,');
 annotation('textbox',[.44  .94 .15 .07],'String',str1,'FitBoxToText','on','LineStyle','none');
  
    for ii=1:length(a)
    str{ii}=([ 'post' num2str(ii) '=' num2str(PrePostG4(a(ii),1)) ',']);
    annotation('textbox',[.53+u(ii) .94 .15 .07],'String',str{ii},'FitBoxToText','on','LineStyle','none');
    hold on
    end
 
end

dim2 = [.4 .4 .27 .07];
 t1=thetaI_run(Cndex);
 t2=thetaIrest(Cndex);
 t3=thetaf_run(Cndex);
 t4=thetafrest(Cndex);
str2=(['\thetaIrun=' num2str(t1) ',  ' '\thetaIrest=' num2str(t2) ', ' '\thetafrun=' num2str(t3) ', ' '\thetafrest=' num2str(t4)]);
annotation('textbox',dim2,'String',str2,'FitBoxToText','on','LineStyle','none');

dim3 = [.4 .37 .27 .07];
 t5=burstIndex_run(Cndex);
 t6=burstIndex_rest(Cndex);
 t7=refractoryT_run(Cndex);
 t8=refractoryT_rest(Cndex);
str3=(['B.I.run=' num2str(t5) ',  ' 'B.I.rest' num2str(t6) ', ' 'R.T.run=' num2str(t7) ', ' 'R.T.rest=' num2str(t8) ', ' 'Cell#=' num2str(Cellg)]);
annotation('textbox',dim3,'String',str3,'FitBoxToText','on','LineStyle','none');


%% Plot acg 
axes('Position',[0.663, 0.03, 0.17, 0.23])
%   bar(acg.acg10t,acg.acg10(:,CellN))
%   xlim([-700 700])
bar(acg.acg1t,acg.acg1rest(:,Cndex))
xlim([-50 50])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
%   bar(acg.acg1t,acg.acg1(:,CellN))
%    xlim([-50 50])

axes('Position',[.663 .26 .17 .12])
box on
bar(acg.acg10t,acg.acg10rest(:,Cndex))
xlim([-700 700])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])

%% Acg run
axes('Position',[0.835, 0.03, 0.17, 0.23])
%   bar(acg.acg10t,acg.acg10(:,CellN))
%   xlim([-700 700])

bar(acg.acg1t,acg.acg1run(:,Cndex))
xlim([-50 50])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
%   bar(acg.acg1t,acg.acg1(:,CellN))
%    xlim([-50 50])

axes('Position',[.835 .26 .17 .12])
box on
bar(acg.acg10t,acg.acg10run(:,Cndex))
xlim([-700 700])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])

set(gcf,'color','w');

%%

print('-depsc',['shank' num2str(shankID) 'cell' num2str(CellN) ])
% print('-dtiff',['shank' num2str(shankID) 'cell' num2str(CellN) ])
% close all

% figure
% plot3((mean(acg.acg1(51:91,:),1)),Frate.run,spkinfo.duration,'*','markersize',5,'LineWidth',1) 
% grid on
% hold on
% plot3((mean(acg.acg1(51:91,(prendx)),1)),Frate.run(prendx),spkinfo.duration(prendx),'go','markersize',8,'LineWidth',2)
% hold on
% plot3((mean(acg.acg1(51:91,(postndx)),1)),Frate.run(postndx),spkinfo.duration(postndx),'kd','markersize',8,'LineWidth',2)
