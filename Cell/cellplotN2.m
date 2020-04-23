function cellplotN2(filename,h,CellID,EEgsampling_rate,numchannel)

shankID=CellID(h,1);
CellN=CellID(h,2);
load([filename '.mat']);
%%
% h=504;
%  shankID=9;
%  CellN=3;
%  Cellg=903;
%%

cd MI_Probs
load(['MI_' num2str(h) '.mat'])
cd ../
[a,b]=max(MI(1:10));
MI_precent=a*100;
phaseprocession(filename,b,Cellinfo.Cellchn_T(h,1),Cellinfo.Cellchn_T(h,2),EEgsampling_rate,numchannel);

%%  Plot MI
 
cd MI_probs

ga_1=.66;
ga_2=.5;
ga_3=.33;
ga_4=.45;
axes('Position',[ga_1 ga_2 ga_3 ga_4])

PhaseFreqVector=1:250;
PhaseFreq_BandWidth=4;
load(['MI_' num2str(h) '.mat'])
plot(PhaseFreqVector+PhaseFreq_BandWidth/2,MI)
hold on
load(['MI_rest_' num2str(h) '.mat'])
plot(PhaseFreqVector+PhaseFreq_BandWidth/2,MI_rest,'k')
hold on
load(['MI_run_' num2str(h) '.mat'])
plot(PhaseFreqVector+PhaseFreq_BandWidth/2,MI_run,'r')

xlim([1 254])
legend('Overal','rest','run','Location','best')
legend('boxoff')
cd ../


%% Plot acg Total
gb_1=.61;
gb_2=.03;
gb_3=0.13;
gb_4=0.25;

axes('Position',[gb_1,gb_2, gb_3, gb_4])
%   bar(acg.acg10t,acg.acg10(:,CellN))
%   xlim([-700 700])
bar(acg.acg1t,acg.acg1(:,h),'FaceColor',[0 0.3 0.9],'EdgeColor',[0 0.3 0.9])
xlim([-50 50])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
%   bar(acg.acg1t,acg.acg1(:,CellN))
%    xlim([-50 50])

gc_1=gb_1;
gc_2=gb_2+gb_4;
gc_3=gb_3;
gc_4=.19;
axes('Position',[gc_1 gc_2 gc_3 gc_4])

box on
bar(acg.acg10t,acg.acg10(:,h),'FaceColor',[0 0.5 0.1],'EdgeColor',[0 0.5 0.1])
xlim([-700 700])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
legend('Overal','Location','northoutside')
% Plot acg rest
gd_1=gb_1+gb_3;
gd_2=.03;
gd_3=gb_3;
gd_4=gb_4;

axes('Position',[gd_1, gd_2, gd_3, gd_4])

%   bar(acg.acg10t,acg.acg10(:,CellN))
%   xlim([-700 700])
bar(acg.acg1t,acg.acg1rest(:,h),'FaceColor',[0 0.3 0.9],'EdgeColor',[0 0.3 0.9])
xlim([-50 50])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
%   bar(acg.acg1t,acg.acg1(:,CellN))
%    xlim([-50 50])

ge_1=gd_1;
ge_2=gd_2+gd_4;
ge_3=gd_3;
ge_4=gc_4;

axes('Position',[ge_1, ge_2, ge_3, ge_4])

box on
bar(acg.acg10t,acg.acg10rest(:,h),'FaceColor',[0 0.5 0.1],'EdgeColor',[0 0.5 0.1])
xlim([-700 700])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
legend('rest','Location','northoutside')
% Acg run

gf_1=gd_1+gd_3;
gf_2=.03;
gf_3=gd_3;
gf_4=gd_4;

axes('Position',[gf_1, gf_2, gf_3, gf_4])

%   bar(acg.acg10t,acg.acg10(:,CellN))
%   xlim([-700 700])

bar(acg.acg1t,acg.acg1run(:,h),'FaceColor',[0 0.3 0.9],'EdgeColor',[0 0.3 0.9])
xlim([-50 50])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
%   bar(acg.acg1t,acg.acg1(:,CellN))
%    xlim([-50 50])
gdi_1=gf_1;
gdi_2=gf_2+gf_4;
gdi_3=gf_3;
gdi_4=gc_4;

axes('Position',[gdi_1,gdi_2,gdi_3,gdi_4])
box on
bar(acg.acg10t,acg.acg10run(:,h),'FaceColor',[0 0.5 0.1],'EdgeColor',[0 0.5 0.1])
xlim([-700 700])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])

legend('run','Location','northoutside')


%% Annotation
g_an1=1-0.17;
g_an2=.93;
g_an3=.2;
g_an4=.07;

dim2 = [g_an1 g_an2 g_an3 g_an4];
str2=(['Cell#= ' num2str(CellID(h,3))]);
annotation('textbox',dim2,'String',str2,'FitBoxToText','on','LineStyle','none');

% cd Path_zone
cd Cue_zone
print('-depsc',['cell' num2str(h) ])
cd ..


% cd Phase_images
% cd (['Shank' num2str(shankID)])
% print('-depsc',['shank' num2str(shankID) 'cell' num2str(CellN) ])
% cd ..
% cd ..