%cellplot
function cellplot(shankID,CellN,Cellg)
% 
 shankID=9;
 CellN=3;
 Cellg=903;

filename1='FM05_1';
filename2='FM05_1_clu_';
% filename3='FM05_1PlaceField.mat';
filename3='FM05_1PlaceField.mat';
load('PrePostG1.mat')
load('PrePostG4.mat')
% load('SpatialC.mat')
load('SpatialC.mat')
load('thetaf.mat')
load('thetaI.mat')
load(filename1)
load(filename3)
load([filename2,num2str(shankID),'.mat'])
ndex=find(clu.g==CellN);
Cndex=find(G==Cellg);

 %% Plot firing field

plotPFSetObj2(filename1,Cndex,xttsc)
% title(['shank' num2str(shankID) 'cell' num2str(CellN) ])

%%  Plot mean fields

axes('Position',[.28 .5 .15 .42])
ndxzero=find(xttsc(:,1)==1);
sets=unique(xttsc(:,4));
binN=100;
for i=1:length(sets)
s=find(xttsc(ndxzero,4)==sets(i));
Set(i,1)=s(1);
Set(i,2)=s(end);
end

cc=reshape(xttsc(:,4+Cndex),binN,142);
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

plot (I2,NS(3,:),'b')
hold on
plot(I2,NS(2,:)+1,'r')
hold on
plot(I2,NS(1,:)+2,'k')

set(gca,'YTickLabel',[])
set(gca,'YTickLabel',[3 2 1])
set(gca,'Ytick',[0.5 1.5 2.5])
% ylabel('Nornalized spk # per set','fontsize',11)

axes('Position',[0.47, 0.5, 0.17, 0.42])
box on
plot (I2,S(3,:),'b')
hold on
plot(I2,S(2,:)+50,'r')
hold on
plot(I2,S(1,:)+100,'k')
% title('Sum')
% ylabel('spik # ','fontsize',11)

axes('Position',[.28 .93 .15 .07])
box on
imagesc(IS)
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])

%% plot celltypes 
axes('Position',[.68, 0.5, 0.32, 0.44])

[~,prendx]=ismember(PrePostG1(:,1),G);
[~,postndx]=ismember(PrePostG1(:,2),G);

plot3(spkinfo.duration,acg.refractoryT,acg.burstIndex,'*','markersize',5,'LineWidth',1) 
grid on
hold on
plot3(spkinfo.duration(prendx),acg.refractoryT(prendx),acg.burstIndex(prendx),'go','markersize',8,'LineWidth',2) 
hold on
plot3(spkinfo.duration(postndx),acg.refractoryT(postndx),acg.burstIndex(postndx),'kd','markersize',8,'LineWidth',2) 
hold on
xlabel('S.D')
ylabel('RT')
% zlabel('B.I')
% set(gcf,'Color','w');
hold on
plot3(spkinfo.duration(Cndex),acg.refractoryT(Cndex),acg.burstIndex(Cndex),'rs','markersize',14,'LineWidth',3)

%% Annotations

%*** New Annotation= 'spk#= 'num2str(length(ndex));

sp=0.11;
u=[0 sp*1 sp*2 sp*3 sp*4 sp*5 sp*6];
 [a,b]=find(Cellg==PrePostG4) ;  
if b==1
 str1=('Cell=Pre,');
 annotation('textbox',[.44 .9 .15 .07],'String',str1,'FitBoxToText','on','LineStyle','none');
 
    for ii=1:length(a)
    str{ii}=([ 'post' num2str(ii) '=' num2str(PrePostG4(a(ii),2)) ',']);
    annotation('textbox',[.53+u(ii) .9 .15 .07],'String',str{ii},'FitBoxToText','on','LineStyle','none');
    hold on
    end
elseif b==2
 str1=('Cell=Post,');
 annotation('textbox',[.44  .9 .15 .07],'String',str1,'FitBoxToText','on','LineStyle','none');
  
    for ii=1:length(a)
    str{ii}=([ 'post' num2str(ii) '=' num2str(PrePostG4(a(ii),1)) ',']);
    annotation('textbox',[.53+u(ii) .9 .15 .07],'String',str{ii},'FitBoxToText','on','LineStyle','none');
    hold on
    end
 
end

dim2 = [.44 .93 .27 .07];
 t1=acg.refractoryT(Cndex);
 t2=acg.burstIndex(Cndex);
 t3=thetaf(Cndex);
 t4=thetaI(Cndex);
str2=(['Burst=' num2str(t2) ',  ' 'Rf.T=' num2str(t1) ', ' '\thetaI=' num2str(t4) ', ' '\thetaf=' num2str(t3) ', ' 'Cell#=' num2str(Cellg)]);
annotation('textbox',dim2,'String',str2,'FitBoxToText','on','LineStyle','none');

 %%  Plot MI
cd MI 
axes('Position',[.03 .03 .19 .42])
% axes('Position',[.28, 0.05, 0.3, 0.4])

load(['MI_' num2str(Cndex) '.mat'])
PhaseFreqVector=1:250;
PhaseFreq_BandWidth=4;
plot(PhaseFreqVector+PhaseFreq_BandWidth/2,MI)
% xlabel('Phase')
% ylabel('MI')
xlim([1 254])
cd ../
%% spatial coverage and firin rate
axes('Position',[0.25, 0.03, 0.285, 0.42])

sety=2;
plot(SpatialC(:,sety+1), Frate.run,'x','markersize',5,'LineWidth',1) 
hold on
plot(SpatialC(postndx,sety+1), Frate.run(postndx),'kd','markersize',5,'LineWidth',2) 
hold on
plot(SpatialC(prendx,sety+1), Frate.run(prendx),'go','markersize',5,'LineWidth',2) 
% xlabel('SpatialC')
% ylabel('Frate.rest')
hold on
sety=1;
plot(SpatialC(Cndex,sety+1), Frate.run(Cndex),'ro','markersize',9,'LineWidth',3) 
sety=2;
hold on
plot(SpatialC(Cndex,sety+1), Frate.run(Cndex),'rd','markersize',9,'LineWidth',3) 
sety=3;
hold on
plot(SpatialC(Cndex,sety+1), Frate.run(Cndex),'rs','markersize',9,'LineWidth',3) 
legend('C','post','pre','S1','S2','S3','Location','northwest')
xlim([0 .9])

%% load spkw
% chorder=[1 8 2 7 3 6 4 5];
% nsamples=40;
% nchannels = length(chorder);
% fp = fopen([filename1 '.spk.' num2str(shankID)], 'r');
% spkW = fread(fp, [nchannels, inf], 'short');
% nspike = size(spkW, 2)/nsamples;
% spkW = reshape(spkW, [nchannels, nsamples, nspike]);
% spkW = spkW(chorder,:,:);
%% plot SPK position

axes('Position',[0.54, 0.03, 0.12, 0.35])
chspace=400;
% for kk=1:8;
% 
% plot(squeeze(spkW(kk,:,ndex))-(kk-1)*chspace,'color',[211./256  211./256   211./256]) 
% hold on
% 
% end
% hold on
for rr=1:8;

plot(squeeze(clu.spk_mean(rr,:,CellN))-(rr-1)*chspace,'color',[0./256  0./256   211./256]) 
hold on

end

set(gca,'YTickLabel',[])
% set(gca,'YTickLabel',[8 7 6 5 4 3 2 1])
% set(gca,'Ytick',[ -7*chspace -6*chspace -5*chspace -4*chspace -3*chspace -2*chspace -1*chspace 0])
% ylabel('Ch # ','fontsize',11)
xlabel('Sample #','fontsize',11)
xlim([0 41])

ylim([-8*chspace chspace/2])


%% Plot acg 
axes('Position',[0.663, 0.03, 0.17, 0.23])
%   bar(acg.acg10t,acg.acg10(:,CellN))
%   xlim([-700 700])
bar(acg.acg1t,acg.acg1(:,Cndex))
xlim([-50 50])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
%   bar(acg.acg1t,acg.acg1(:,CellN))
%    xlim([-50 50])

axes('Position',[.663 .26 .17 .12])
box on
bar(acg.acg10t,acg.acg10(:,Cndex))
xlim([-700 700])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
%% Synaptic connections

sp2=0.125;
uu=[0 sp2*1 sp2*2 sp2*3];
q=length(a);
q1=[];
if q>3
    q1=3;
else
    q1=q;
end

if isempty(b)==0
   
    for i=1:q1
    samplerate =(spkinfo.samplerate);
    [ccg,tccg]=CCG(spk.i,spk.g,0.001*samplerate,10,samplerate,[PrePostG4(a(i),1) ; PrePostG4(a(i),2)],'count');
    axes('Position',[0.835, 0.26-uu(i), 0.18, 0.12])
    bar(tccg,ccg(:,1,2));
    xlim([tccg(1) tccg(end)]);
    set(gca,'xtick',[],'ytick',[])
    hold on
    plot([0 0],[0 max(ccg(:,1,2))],'r')
    hold off
    end
    
end
%% plot cell position
axes('Position',[.55 .12 .44 .6])
box on
Clayout_2('PrePostG1.mat',Cellg)
ylim([0 200])
set(gcf,'color','w');

%%
% set(subplot(2,2,1), 'Position', [0.1, 0.5, 0.4, 0.4])
% set(subplot(2,2,2), 'Position', [0.59, 0.5, 0.4, 0.4])
% set(subplot(2,2,3), 'Position', [0.1, 0.05, 0.4, 0.4])
% set(subplot(2,2,4), 'Position', [0.59, 0.05, 0.4, 0.4])
print('-depsc',['shank' num2str(shankID) 'cell' num2str(CellN) ])
% print('-dtiff',['shank' num2str(shankID) 'cell' num2str(CellN) ])
% close all

% s11=sum(set1,2);
% s22=sum(set2,2);
% s33=sum(set3,2);
% T1=sum(s11(Rew(1,1):Rew(1,2),:),1)
% T2=sum(s22(Rew(2,1):Rew(2,2),:),1)
% T3=sum(s33(Rew(3,1):Rew(3,2),:),1)
% allset2=cat(2,s11./T1,s22./T2,s33./T3);
