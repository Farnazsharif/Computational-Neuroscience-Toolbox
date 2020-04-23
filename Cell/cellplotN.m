%cellplot
function cellplotN(shankID,CellN,Cellg,filename1,CA)
% This section shoulde be editted that be useable only by mouse name
%  shankID=7;
%  CellN=2;
%  Cellg=702;
% filename1='FM05_1';

filename2=([filename1 '_clu_']);
filename3=([filename1 'PlaceField.mat']);
filename4=([filename1 'PlaceField_Rset.mat']);
filename5=([filename1 'PlaceField_Total.mat']);

load(filename1)
load(filename3)
load(filename4)
load(filename5)
load([filename2,num2str(shankID),'.mat'])

ndex=find(clu.g==CellN);
Cndex=find(G==Cellg);
PrePostG=Cellinfo.PrePostG_T;

%% Plot firing field
gu1=.04; gu2=0.5; gu3=0.23; gu4=0.45;
% plotPFSetObj2(filename1,Cndex,xttsc)
plotPFSetObj_whole(filename1,Cndex,xttsc,gu1, gu2, gu3, gu4)
% title(['shank' num2str(shankID) 'cell' num2str(CellN) ])
%%  Plot mean fields
g_cons_1=0.055; % space betw the figures;
g_cons_2=0.07; % thikness of the imagesc box;
g_1=gu1+gu3+g_cons_1;
g_2=gu2;
g_3=gu3/2;
g_4=(1-gu2-g_cons_2-0.005)/2.35;

axes('Position',[g_1 g_2 g_3 g_4])
ndxzero=find(xttsc(:,1)==1);
sets=unique(xttsc(:,4));
binN=100;
for i=1:length(sets)
    s=find(xttsc(ndxzero,4)==sets(i));
    Set(i,1)=s(1);
    Set(i,2)=s(end);
end

cc=reshape(xttsc(:,4+Cndex),binN,length(xttsc)/binN);
scale=10;
I=1:2:200;
I2=imresize(I,[1 length(I)*scale],'lanczos3');

for i =1 :length(sets)
    se=cc(:,Set(i,1):Set(i,2));
    spk_perset(i)=sum(sum(se));
    S(i,:)=imresize(sum(se,2)',[1 binN*scale],'lanczos3');
    NS(i,:)=matnorm(imresize(sum(se,2)',[1 binN*scale],'lanczos3'),2);
    IS(i,:)=matnorm(sum(se,2),1);
end
spkN_run=sum(spk_perset);
plot (I2,NS(3,:),'b')
hold on
plot(I2,NS(2,:)+1,'r')
hold on
plot(I2,NS(1,:)+2,'k')

set(gca,'YTickLabel',[])
set(gca,'YTickLabel',[fliplr(spk_perset)])
set(gca,'Ytick',[0.5 1.5 2.5])
set(gca,'XTickLabel',[])

axes('Position',[g_1 g_2+g_4 g_3 g_cons_2])
box on
imagesc(IS)
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])

% g_cons_1=0.055; % space betw the figures;
% g_cons_2=0.07; % thikness of the imagesc box;
% g_13=g_9+g_11+g_cons_1;
g_13=g_1;
g_14=g_2+g_4+g_cons_2+0.001;
g_15=gu3/2;
g_16=(1-gu2-g_cons_2-0.005)/2.35;
axes('Position',[g_13 g_14 g_15 g_16])

cc=[];
I=[];
I2=[];
cc=reshape(xttscT(:,4+Cndex),binN,length(xttsc)/binN);
scale=10;
I=1:2:200;
I2=imresize(I,[1 length(I)*scale],'lanczos3');
se=[];
spk_perset=[];
S=[];
NS=[];
IS=[];
for i =1 :length(sets)
    se=cc(:,Set(i,1):Set(i,2));
    spk_perset(i)=sum(sum(se));
    S(i,:)=imresize(sum(se,2)',[1 binN*scale],'lanczos3');
    NS(i,:)=matnorm(imresize(sum(se,2)',[1 binN*scale],'lanczos3'),2);
    IS(i,:)=matnorm(sum(se,2),1);
end
spkN_rest=sum(spk_perset);
plot (I2,NS(3,:),'b')
hold on
plot(I2,NS(2,:)+1,'r')
hold on
plot(I2,NS(1,:)+2,'k')

set(gca,'YTickLabel',[])
set(gca,'YTickLabel',[fliplr(spk_perset)])
set(gca,'Ytick',[0.5 1.5 2.5])
set(gca,'XTickLabel',[])


axes('Position',[g_13 g_14+g_16 g_3 g_cons_2])
box on
imagesc(IS)
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])

%%
% plotPFSetObj2(filename1,Cndex,xttscR)
g_5=g_1+g_3+0.038;
g_6=gu2;
g_7=gu3*(2/3);
g_8=gu4;
plotPFSetObj_whole(filename1,Cndex,xttscR,g_5, g_6, g_7, g_8)
%%
g_9=g_5+g_7+0.038;
g_10=gu2;
g_11=gu3*(2/3);
g_12=gu4;

plotPFSetObj_whole(filename1,Cndex,xttscT,g_9, g_10, g_11, g_12)

% title(['shank' num2str(shankID) 'cell' num2str(CellN) ])

%% Firing rate and refractory period

gdb_1=0.03;
gdb_2=.03;
gdb_3=0.31;
gdb_4=0.2;

gda_1=gdb_1;
gda_2=gdb_2+gdb_4+0.035;
gda_3=gdb_3;
gda_4=gdb_4;

[~,prendx]=ismember(PrePostG(:,1),G);
[~,postndx]=ismember(PrePostG(:,2),G);

% [~,prendx]=ismember(Cellinfo.PrePostG_T(:,1),G);
% [~,postndx]=ismember(Cellinfo.PrePostG_T(:,2),G);



FrRun=Frate.run;
FrRest=Frate.rest;

RefRun=spkinfo.refractoryTrun;
RefRest=spkinfo.refractoryTrest;

axes('Position',[gda_1 gda_2 gda_3 gda_4])
plot(RefRun,FrRun,'*')
hold on
plot(RefRun(prendx),FrRun(prendx),'go','markersize',8,'LineWidth',2)
hold on
plot(RefRun(postndx),FrRun(postndx),'kd','markersize',8,'LineWidth',2)
hold on
plot(RefRun(Cndex),FrRun(Cndex),'r*','markersize',8,'LineWidth',2)

Ref=(RefRun(Cndex)/RefRest(Cndex));
Fr=FrRun(Cndex)/FrRest(Cndex);
str_gda1 = (['Ref=' num2str(Ref) '   Fr=' num2str(Fr) ]);
dim_gda=[(gda_3)*.43 gda_2+gda_4/1.5  gda_4  0.05] ;
annotation('textbox',dim_gda,'String',str_gda1,'FitBoxToText','on','LineStyle','none');

axes('Position',[gdb_1 gdb_2 gdb_3 gdb_4])
plot(RefRest,FrRest,'*')
hold on
plot(RefRest(prendx),FrRest(prendx),'go','markersize',8,'LineWidth',2)
hold on
plot(RefRest(postndx),FrRest(postndx),'kd','markersize',8,'LineWidth',2)
hold on
plot(RefRest(Cndex),FrRest(Cndex),'r*','markersize',8,'LineWidth',2)

%% plot celltypes

gdk_1=gda_1+gda_3+0.05;
gdk_2=0.03;
gdk_3= 0.35;
gdk_4=0.4;

axes('Position',[gdk_1, gdk_2, gdk_3, gdk_4])



plot3(spkinfo.duration,acg.refractoryT,acg.burstIndex,'*','markersize',5,'LineWidth',1)
grid on
hold on
plot3(spkinfo.duration(prendx),acg.refractoryT(prendx),acg.burstIndex(prendx),'go','markersize',8,'LineWidth',2)
hold on
plot3(spkinfo.duration(postndx),acg.refractoryT(postndx),acg.burstIndex(postndx),'kd','markersize',8,'LineWidth',2)
hold on
xlabel('S.D')
ylabel('RT')

hold on
plot3(spkinfo.duration(Cndex),acg.refractoryT(Cndex),acg.burstIndex(Cndex),'rs','markersize',14,'LineWidth',3)

%% plot SPK position
gdc_1=gdk_1+gdk_3+0.005;
gdc_2=.03;
gdc_3=0.08;
gdc_4=0.33;
axes('Position',[gdc_1,gdc_2,gdc_3,gdc_4])

chspace=400;
for rr=1:8;
    plot(squeeze(clu.spk_mean(rr,:,CellN))-(rr-1)*chspace,'color',[0./256  0./256   211./256])
    hold on
end

set(gca,'YTickLabel',[])
xlabel('Sample #','fontsize',11)
xlim([0 41])
ylim([-8*chspace chspace/2])

%% Annotation
g_an1=1-0.17;
g_an2=.93;
g_an3=.2;
g_an4=.07;

dim2 = [g_an1 g_an2 g_an3 g_an4];
str2=(['Cell#= ' num2str(Cellg)]);
annotation('textbox',dim2,'String',str2,'FitBoxToText','on','LineStyle','none');
dim3 = [g_an1 g_an2-0.03 g_an3 g_an4];
str3=([ 'spk#= ' num2str(length(ndex)) ]);
annotation('textbox',dim3,'String',str3,'FitBoxToText','on','LineStyle','none');
dim4 = [g_an1 g_an2-(0.03)*2 g_an3 g_an4];
str4=([ 'spk- rest#= ' num2str(spkN_rest) ]);
annotation('textbox',dim4,'String',str4,'FitBoxToText','on','LineStyle','none');
dim5 = [g_an1 g_an2-(0.03)*3 g_an3 g_an4];
str5=([ 'spk- run#= ' num2str(spkN_run) ]);
annotation('textbox',dim5,'String',str5,'FitBoxToText','on','LineStyle','none');

% Annotations PrePost

[a,b]=find(Cellg==Cellinfo.PrePostG_T) ;
if b==1
    str1=('Cell=Pre,');
    annotation('textbox',[g_an1 g_an2-(0.03)*4 g_an3 g_an4],'String',str1,'FitBoxToText','on','LineStyle','none');
    
    for ii=1:length(a)
        str{ii}=([ 'post' num2str(ii) '=' num2str(Cellinfo.PrePostG_T(a(ii),2)) ',']);
        annotation('textbox',[g_an1 g_an2-(0.03)*(4+ii) g_an3 g_an4],'String',str{ii},'FitBoxToText','on','LineStyle','none');
        hold on
    end
elseif b==2
    str1=('Cell=Post,');
    annotation('textbox',[g_an1 g_an2-(0.03)*4 g_an3 g_an4],'String',str1,'FitBoxToText','on','LineStyle','none');
    
    for ii=1:length(a)
        str{ii}=([ 'Pre' num2str(ii) '=' num2str(Cellinfo.PrePostG_T(a(ii),1)) ',']);
        annotation('textbox',[g_an1 g_an2-(0.03)*(4+ii) g_an3 g_an4],'String',str{ii},'FitBoxToText','on','LineStyle','none');
        hold on
    end
    
end

%% Synaptic connections
gdl_1=gdc_1+gdc_3+0.01;
gdl_2=0.555;
gdl_3=0.16;
gdl_4=0.11;
%
% sp2=gdl_4;
% uu=[0 sp2*1 sp2*2 sp2*3 sp2*4 sp2*5];
q=length(a);
q1=[];
if q>6
    q1=6;
else
    q1=q;
end

if isempty(b)==0
    
    for i=1:q1
        samplerate =(spkinfo.samplerate);
        [ccg,tccg]=CCG(spk.i,spk.g,0.001*samplerate,10,samplerate,[Cellinfo.PrePostG_T(a(i),1) ; Cellinfo.PrePostG_T(a(i),2)],'count');
        axes('Position',[gdl_1, gdl_2-(gdl_4*(i-1)), gdl_3, gdl_4])
        bar(tccg,ccg(:,1,2));
        xlim([tccg(1) tccg(end)]);
        set(gca,'xtick',[],'ytick',[])
        hold on
        plot([0 0],[0 max(ccg(:,1,2))],'r')
        hold off
    end
    
end
%
%% plot cell position

gdj_1=gda_1+gda_3+0.005;
gdj_2=0.12;
gdj_3=gdk_1+gdc_3+0.01;
gdj_4=0.6;

axes('Position',[gdj_1, gdj_2, gdj_3, gdj_4])
% axes('Position',[.55 .12 .44 .6])
box on
%%
[ap,~]=find(fix(Cellinfo.PrePostG_T(:,1)/100)<8);

if CA==1
    
    load('Ripple.mat')
    
    if shankID <9 % this one was wrong
        PrePostG_Prob=PrePostG(ap(1):ap(end),:);
        Rf_ch_shank=Rf_ch_shank(1:8,:)'';
    else
        PrePostG_Prob=PrePostG(ap(end)+1:length(PrePostG),:);
        Rf_ch_shank=Rf_ch_shank(9:16,:)'';
    end
    
    Rf_ch_shank(:,3)=Rf_ch_shank(:,3)/(max(Rf_ch_shank(:,3)));
    
    
else
    
    Rf_ch_shank=0;
    
    if shankID <9
        PrePostG_Prob=PrePostG(ap(1):ap(end),:);
        
    else
        PrePostG_Prob=PrePostG(ap(end)+1:length(PrePostG),:);
        
    end
  
end
%%
Clayout_2(PrePostG_Prob,filename1,Cellg,Rf_ch_shank,CA)
ylim([0 210])

%%
set(gcf,'color','w');

cd Placefield_images
cd (['Shank' num2str(shankID)])
print('-depsc',['shank' num2str(shankID) 'cell' num2str(CellN) ])
cd ..
cd ..







