function phaseprocession(filename,PhaseFreqVector_a,cellN,ch,EEgsampling_rate,numchannel)
% phaseprocession('FM05_1',5,1501,113);
%%
% phaseprocession('FM05_1',b,Cellinfo_Prob_1.Cellchn(h,1),1,Cellinfo_Prob_1.Cellchn(h,2));

% PhaseFreqVector_a=b;
% Cellinfo.Cellchn_T(h,1),Cellinfo.Cellchn_T(h,2);
%%
PhaseFreqVector_b=PhaseFreqVector_a;
load([filename '.mat']);
% ss=find(G==cellN);
% [Cellchn]=Cellsite(ProbN);

% ch=Cellchn(ss,2);

%%
PhaseFreq_BandWidth=4;
Pf1 = PhaseFreqVector_a-(PhaseFreq_BandWidth)./2;
Pf2 = Pf1 + (PhaseFreq_BandWidth)./2;
eeg= readmulti([filename '.lfp'],numchannel,ch);
tic
%%
%[phase,g]=spkphase(spk.i,spk.g,cellN,spkinfo.samplerate,eeg,[Pf1 Pf2],spkinfo.samplerate./25);
% [MI,n,phase]=cellph(cellN,ch,filename,PhaseFreqVector_a,PhaseFreqVector_b);
[MI,MI_rest, MI_run,n,n_rest,n_run,phase,phase_rest,phase_run]=cellph_Overal(cellN,filename,PhaseFreqVector_a,PhaseFreqVector_b,EEgsampling_rate,eeg);
toc
                                           
%% Phase Procession
[spikeT,spikeG] = selectgroup(spk.i,spk.g,cellN);
for i=1:length(spikeT)
[~,ndx(i)]=min(abs(behav.TXDTS(:,1)-(spikeT(i)./spkinfo.samplerate)));
end
spkx=behav.TXDTS(ndx,2);
%%
% figure
gua_1=.03;
gua_2=.69;
gua_3=.38;
gua_4=.27;
axes('Position',[gua_1 gua_2 gua_3 gua_4])

plot(spkx,phase,'.b',spkx,phase+2*pi,'.b','markersize',5,'Color',[0 0.5 0.8])
set(gca,'YTickLabel',[])
set(gca,'Ytick',[0 pi 2*pi  3*pi 4*pi])
set(gca,'YTickLabel',{'0','pi','2pi','3pi','4pi'})
ylabel('Phase')
% xlabel('X(cm)')
title(['Overal   Frequency' num2str(Pf1 + (PhaseFreq_BandWidth./2))  'Hz lfp'])
% xlim([0 200])
ylim([0 4*pi])
set(gcf,'color','w')
set(gca,'XTickLabel',[])
%%
[tt_rest]=selectt(spikeT,spikeG,behav.restT,spkinfo.samplerate);
TF_re=isempty(tt_rest);

if TF_re==1
    spkx_rest=[];
else
   
for i=1:length(tt_rest)
[~,ndx_rest(i)]=min(abs(behav.TXDTS(:,1)-(tt_rest(i)./spkinfo.samplerate)));
end
spkx_rest=behav.TXDTS(ndx_rest,2);

end

%%
gub_1=gua_1;
gub_2=1-gua_2+0.05;
gub_3=gua_3;
gub_4=gua_4;
axes('Position',[gub_1 gub_2 gub_3 gub_4])

plot(spkx_rest,phase_rest,'.b',spkx_rest,phase_rest+2*pi,'.b','markersize',5,'Color',[0 0.5 0.8])
set(gca,'YTickLabel',[])
set(gca,'Ytick',[0 pi 2*pi  3*pi 4*pi])
set(gca,'YTickLabel',{'0','pi','2pi','3pi','4pi'})
ylabel('Phase')
% xlabel('X(cm)')
title(['Rest    Frequency' num2str(Pf1 + (PhaseFreq_BandWidth./2))  'Hz lfp'])
% xlim([0 200])
ylim([0 4*pi])
set(gcf,'color','w')
set(gca,'XTickLabel',[])
%%
[tt_run]=selectt(spikeT,spikeG,behav.runT,spkinfo.samplerate);
TF_r=isempty(tt_run);
if TF_r==1
    spkx_run=[];
else


for i=1:length(tt_run)
[~,ndx_run(i)]=min(abs(behav.TXDTS(:,1)-(tt_run(i)./spkinfo.samplerate)));
end
spkx_run=behav.TXDTS(ndx_run,2);

end
%%

guc_1=gua_1;
guc_2=.03;
guc_3=gua_3;
guc_4=gua_4;
axes('Position',[guc_1 guc_2 guc_3 guc_4])

plot(spkx_run,phase_run,'.',spkx_run,phase_run+2*pi,'.b','markersize',5,'Color',[0 0.5 0.8])
set(gca,'YTickLabel',[])
set(gca,'Ytick',[0 pi 2*pi  3*pi 4*pi])
set(gca,'YTickLabel',{'0','pi','2pi','3pi','4pi'})
ylabel('Phase')
% xlabel('X(cm)')
title(['Run    Frequency' num2str(Pf1 + (PhaseFreq_BandWidth./2))  'Hz lfp'])
% xlim([0 200])
ylim([0 4*pi])
set(gcf,'color','w')

%% Phase Histogram
[a,b]=max(n);
binrange=10:20:720;

gud_1=gua_1+gua_3+0.02;
gud_2=gua_2;
gud_3=0.17;
gud_4=gua_4-0.05;
axes('Position',[gud_1 gud_2 gud_3 gud_4])

[a1,~]=size(n); 
if a1==1
    n=n';
end
bar(binrange,[n/a;n/a],'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9])

ax1 = gca;
xlim([0 720])
set(ax1,'XTickLabel',[])
set(ax1,'XTickLabel',[0 180 360  540 720])
set(ax1,'Xtick',[0 180 360  540 720])
set(gca,'YTickLabel',[])

x1 = [0:0.1:4*pi];
y1 = 2*cos(x1)+17;
ax1_pos = get(ax1,'Position'); 
ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');

line(x1,y1,'Parent',ax2,'Color',[1 .5 0],'LineWidth',2)
ylim([0 20])
xlim([0 4*pi])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'YTickLabel','')
axis off
% xlabel('Phase')
% ylabel('Spik#')
% title(['Overal    Frequency' num2str(Pf1 + (PhaseFreq_BandWidth./2)) 'Hz lfp'])

dim = [gud_1+0.04 gud_2+gud_4+0.02 gud_3 .02];
str = ([ num2str(MI*100) '   ph=' num2str(binrange(b))]  );
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none');
set(gcf,'color','w')
% polt bar
gubar1_1=gud_1+0.01; gubar1_2=gud_2+gud_4+0.005; gubar1_3=0.03; gubar1_4=0.05;
axes('Position',[gubar1_1 gubar1_2 gubar1_3 gubar1_4]) ;% for ripple

bar(MI,'b','r')
hold on
xL = get(gca,'XLim');
line([0.5 1.5],[0.01 0.01] ,xL,'Color','k');
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
% [MImax,~]=max([MI MI_rest MI_run]);

if isnan(MI)==1||isnan(MI_rest)==1||isnan( MI_run)==1
YL = get(gca,'YLim');
YL = YL(end)+0.01 ;
else
YL = 2*(max([MI MI_rest MI_run]));  
end
ylim([0 YL])
%%
gue_1=gud_1;
gue_2=gub_2;
gue_3=gud_3;
gue_4=gud_4;
axes('Position',[gue_1 gue_2 gue_3 gue_4])

[a1,~]=size(n_rest); 
if a1==1
    n_rest=n_rest';
end


[a_rest,b_rest]=max(n_rest);
bar(binrange,[n_rest/a_rest;n_rest/a_rest],'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9])

ax1 = gca;
xlim([0 720])
set(ax1,'XTickLabel',[])
set(ax1,'XTickLabel',[0 180 360  540 720])
set(ax1,'Xtick',[0 180 360  540 720])
set(gca,'YTickLabel',[])

ax1_pos = get(ax1,'Position'); 
ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');

line(x1,y1,'Parent',ax2,'Color',[1 .5 0],'LineWidth',2)
ylim([0 20])
xlim([0 4*pi])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'YTickLabel','')
axis off


dim = [gud_1+0.04 gue_2+gue_4+0.02 gud_3 .02];
str = ([ num2str(MI_rest*100) '  ph=' num2str(binrange(b_rest))]  );
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none');
set(gcf,'color','w')

% polt bar
gubar2_1=gud_1+0.01; gubar2_2=gue_2+gue_4+0.005; gubar2_3=0.03; gubar2_4=0.05;
axes('Position',[gubar2_1 gubar2_2 gubar2_3 gubar2_4]) ;% for ripple

bar(MI_rest,'b','r')
hold on
xL = get(gca,'XLim');
line([0.5 1.5],[0.01 0.01] ,xL,'Color','k');
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
ylim([0 YL])

%%
guf_1=gud_1;
guf_2=0.03;
guf_3=gud_3;
guf_4=gud_4;
axes('Position',[guf_1 guf_2 guf_3 guf_4])

[a1,~]=size(n_run); 
if a1==1
    n_run=n_run';
end

[a_run,b_run]=max(n_run);
bar(binrange,[n_run/a_run;n_run/a_run],'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9])

ax1 = gca;
xlim([0 720])
set(ax1,'XTickLabel',[])
set(ax1,'XTickLabel',[0 180 360  540 720])
set(ax1,'Xtick',[0 180 360  540 720])
set(gca,'YTickLabel',[])

ax1_pos = get(ax1,'Position'); 
ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');

line(x1,y1,'Parent',ax2,'Color',[1 .5 0],'LineWidth',2)
ylim([0 20])
xlim([0 4*pi])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'YTickLabel','')
axis off

dim = [gud_1+0.04 guf_2+guf_4+0.02 gud_3 .02];
str = ([num2str(MI_run*100) '   ph=' num2str(binrange(b_run))]  );
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none');
set(gcf,'color','w')

gubar2_1=gud_1+0.01; gubar2_2=guf_2+guf_4+0.005; gubar2_3=0.03; gubar2_4=0.05;
axes('Position',[gubar2_1 gubar2_2 gubar2_3 gubar2_4]) ;% for ripple

bar(MI_run,'b','r')
hold on
xL = get(gca,'XLim');
line([0.5 1.5],[0.01 0.01] ,xL,'Color','k');
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
ylim([0 YL])


