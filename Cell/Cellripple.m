%%
function Cellripple(filename,Winrange,CellID,h,EEGsamplerate,recording,ChN)
% ch_Ref=105;
% Winrange=100;
% filename='FM05_1';
% load('Rt4.mat')
% Rt=Rt4;
% h=546;

shankID=CellID(h,1);
CellN=CellID(h,2);
Cellg=CellID(h,3);

%%
load([filename '.mat']);
ch=Cellinfo.Cellchn_T(h,2);
EEGsamplerate=(spkinfo.samplerate./30);
[spikeT,~] = selectgroup(spk.i,spk.g,Cellg);
spikeTime=spikeT./spkinfo.samplerate;
TimeWindow= [-Winrange/EEGsamplerate Winrange/EEGsamplerate];
binsize=1/EEGsamplerate;
%%
load('Ripple.mat')
if isempty(strfind(recording,'c'))==1
R_ndex=find(ismember(Ripple,RippleT_Single)==1); % if Acute
Rt=Rt(R_ndex,:);
ch_Ref=Rf_ch_shank(shankID,1);
end 
ch_Ref=Rf_ch_shank_pow(shankID,6);

%%
eventTime=Rt(:,ch)./EEGsamplerate;
[~,avgPSTH,~,trialPSTH,Ndxspk]=Tperievent(spikeTime,eventTime,TimeWindow,binsize);

eventTime=Rt(:,ch_Ref)./EEGsamplerate;
[~,avgPSTH_Ref,~,trialPSTH_Ref,Ndxspk_Ref]=Tperievent(spikeTime,eventTime,TimeWindow,binsize);

%% GUI
% spikes raster plot(1)
ga1=0.04; ga2=0.68; ga3=0.3; ga4=0.3;
gc1=ga1+ga3-ga3/3; gc2=ga2+ga4-ga4/3+0.02; gc3=ga3/3; gc4=ga4/3;

% Rpple plots
gb1=ga1; gb2=ga2-ga4-0.005; gb3=ga3; gb4=ga4;

% spikes raster plot(2)
gd1=ga1+ga3+0.01; gd2=ga2; gd3=ga3; gd4=ga4;
ge1=gd1+gd3-gd3/3; ge2=gd2+gd4-gd4/3+0.02; ge3=gc3; ge4=gc4;

% Ripple PLot(2)
gf1=gd1; gf2=gb2; gf3=ga3; gf4=ga4;

% Phase and Spkx
gu1_1=ga1; gu1_2=.03; gu1_3=ga3; gu1_4=ga4-0.025;
gu2_1=ga1+ga3*2+0.1; gu2_2=.5 ;  gu2_3=0.25  ;  gu2_4=.41 ;
gu3_1=ga1+ga3+0.01 ;gu3_2=gu1_2; gu3_3=gu1_3 ;gu3_4=gu1_4;
gu4_1=gu2_1 ;gu4_2=.03 ;gu4_3= gu2_3 ;gu4_4=gu2_4 ;

% MI Bar 
gu5_1=ga1+0.05; gu5_2=.31; gu5_3=0.03; gu5_4=0.05;
gu6_1=gu3_1+0.05; gu6_2=gu5_2; gu6_3=gu5_3; gu6_4=gu5_4;

% Annotation
g_an1=1-0.17; g_an2=.91; g_an3=.2; g_an4=.07;
g_ab1=ga1; g_ab2=ga2+ga4+0.02; g_ab3=0.1; g_ab4=0.0105;
g_ac1=gd1; g_ac2=gd2+gd3+0.02; g_ac3=g_ab3; g_ac4=g_ab4;
%%
eeg= readmulti([filename '.lfp'], ChN);
%%
Gu1= [gu1_1 gu1_2 gu1_3 gu1_4]; Gu2=[gu2_1 gu2_2 gu2_3 gu2_4 ];
[~,MI,Avraged_R]=Ripplephase(filename,ch,Ndxspk,Cellg,spikeTime,Rt,Winrange,EEGsamplerate,eeg,Gu1,Gu2);

Gu3=[gu3_1 gu3_2 gu3_3 gu3_4]; Gu4=[gu4_1 gu4_2 gu4_3 gu4_4];
[~,MI_Ref,Avraged_R_Ref]=Ripplephase(filename,ch_Ref,Ndxspk_Ref,Cellg,spikeTime,Rt,Winrange,EEGsamplerate,eeg,Gu3,Gu4);

%% Plot MI bar 

axes('Position',[gu5_1 gu5_2 gu5_3 gu5_4]) ;% for ripple

bar(MI,'b','r')
hold on
xL = get(gca,'XLim');
line([0.5 1.5],[0.01 0.01] ,xL,'Color','k');
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])

if isnan(MI)==1||isnan(MI_Ref)==1
YL = get(gca,'YLim');
YL = YL(end)+0.01 ;
else
YL = 2*(max([MI MI_Ref]));  
end

ylim([0 YL])

axes('Position',[gu6_1 gu6_2 gu6_3 gu6_4]) ;% for ripple

bar(MI_Ref,'b','r')
hold on
xL = get(gca,'XLim');
line([0.5 1.5],[0.01 0.01] ,xL,'Color','k');
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
ylim([0 YL])

%% plot ripples
ax1=axes('Position',[gb1 gb2 gb3 gb4]) ;

imagesc((1:201),(-.5:.5),Avraged_R)
colormap jet
freezeColors(ax1)

hold on
plot(Avraged_R,'k','LineWidth',2)
hold off

set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'Xtick',[1 101 201])
set(gca,'XTickLabel',{'-100 ms','0','100 ms'})
axis xy
% % plot ripples (2)
ax2=axes('Position',[gf1 gf2 gf3 gf4]) ;% for ripple

imagesc((1:201),(-.5:.5),Avraged_R_Ref)
colormap jet
hold on
plot(Avraged_R_Ref,'k','LineWidth',2)

hold off
freezeColors(ax2)

set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
axis xy

%% spikes raster plot(1)

ax3=axes('Position',[ga1 ga2 ga3 ga4]); 

map_raster = [0,.5,.5 ;1,1,1];
% colormap(ax3,map_raster)
imagesc((1:201),(1:171),trialPSTH)
colormap summer
freezeColors(ax3)

hold on

yL = get(gca,'YLim');
line([101 101],yL,'Color','w');

set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
axis xy


axes('Position',[gc1 gc2 gc3 gc4]) ;

plot((1:201),avgPSTH/binsize,'-','Color',[0 0.5 0.7])
hold on
yL = get(gca,'YLim');
line([101 101],yL,'Color','r');

set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
xlim ([1 201])

% spikes raster plot(2)

ax4=axes('Position',[gd1 gd2 gd3 gd4]); % for ripple

% colormap(ax4,map_raster)
imagesc((1:201),(1:171),trialPSTH_Ref)
colormap summer

freezeColors(ax4)

hold on
yL = get(gca,'YLim');
line([101 101],yL,'Color','w');

set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
axis xy


axes('Position',[ge1 ge2 ge3 ge4]) ;

plot((1:201),avgPSTH_Ref/binsize,'-','Color',[0 0.5 0.7])

hold on
yL = get(gca,'YLim');
line([101 101],yL,'Color','r');
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
xlim ([1 201])

%%

dim2 = [g_an1 g_an2 g_an3 g_an4];
str2=(['Cell#= ' num2str(Cellg)]);
annotation('textbox',dim2,'String',str2,'FitBoxToText','on','LineStyle','none');

dim3 = [g_ab1 g_ab2 g_ab3 g_ab4];
str3=(['Ch = ' num2str(ch)]);
annotation('textbox',dim3,'String',str3,'FitBoxToText','on','LineStyle','none','FontWeight','bold');

dim4 = [g_ac1 g_ac2 g_ac3 g_ac4];
str4=(['Ch -  Ref = ' num2str(ch_Ref)]);
annotation('textbox',dim4,'String',str4,'FitBoxToText','on','LineStyle','none','FontWeight','bold');
%%
cd Ripples_images
cd (['Shank' num2str(shankID)])
print('-djpeg',['shank' num2str(shankID) 'cell' num2str(CellN) ])
% savefig(['shank' num2str(shankID) 'cell' num2str(CellN)])
cd ..
cd ..