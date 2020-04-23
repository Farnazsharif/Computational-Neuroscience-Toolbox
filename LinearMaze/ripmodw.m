
% This will display some useful characteristic to classify cells
% USAGE:
%     profile2classify(filename,cluster_group)
%
% INPUTS:
%     filename: base name of the .kwik and .kwx file to be read
%     cluster_group: vector with the "name" of the different
%     cluster_group, it can either start with 0 or 1,
%     %example: filename = 'AB1'
%               cluster_group = 1:8
% 
% Aza 2015

%NOTE: This program needs to be feed with ripples evt calculated with
%Eran's method, for different location as triggered ripples we'll calculate
%modulation of different location cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AB4 01032013
filebase = ['D:\aza\analysis\data' '\AFR2' '\210214' ];
ratname='\AFR2_210214';
cluster_group=1:9;
% mkdir([filebase '\general\'],'profiles');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load cells organize per shank (cluster group)
[spiket, spikeind, numclus, iEleClu] = ReadCluRes([filebase '\clures' ratname],cluster_group);
%Reading neurons classification and spikes
whowhat = readwhowhat([filebase '\general' ratname '.xlsx']);
clc

%PC analysis for... ca1
ca1 = whowhat.id(find(whowhat.region==1));
ca1pyrlayer = ca1(find(whowhat.layer(ca1)==1));
ca1pyrlayerpyrcell = ca1pyrlayer(find(whowhat.type(ca1pyrlayer)==1));
ca1pyrlayerintcell = ca1pyrlayer(find(whowhat.type(ca1pyrlayer)==2));

%PC analysis for... ca2
ca2 = whowhat.id(find(whowhat.region==2));
ca2pyrlayer = ca2;
ca2pyrlayerpyrcell = ca2pyrlayer(find(whowhat.type(ca2pyrlayer)==1));
ca2pyrlayerintcell = ca2pyrlayer(find(whowhat.type(ca2pyrlayer)==2));

%PC analysis for... ca3
ca3 = whowhat.id(find(whowhat.region==3));
ca3pyrlayer = ca3; %todo int and pyr
ca3pyrlayerpyrcell = ca3pyrlayer(find(whowhat.type(ca3pyrlayer)==1));
ca3pyrlayerintcell = ca3pyrlayer(find(whowhat.type(ca3pyrlayer)==2));

%PC analysis for... DG
DG = whowhat.id(find(whowhat.region==4));
dg = DG; %todo int and pyr
dgexc = dg(find(whowhat.type(dg)==1));
dginh = dg(find(whowhat.type(dg)==2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load ripples evt detected with Eran's method, Antonio did it!
%this matrix will include ripples detected in all shanks, just need to use

% evt = load(['Z:\AYA\data\AB4\jdl512_1_m\01032013\ripples3SD' '\AB4_030113ripsALL_3SD.mat']);
evtca1 = load(['Z:\AYA\data\AFR2\210214\ca1' '\AFR2_210214ripsALL_3SD.mat']);
evtca2 = load(['Z:\AYA\data\AFR2\210214\ca2' '\AFR2_210214ripsALL_3SD.mat']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Separate spikes for task and homecage
SR = 20000; Fs = 1250; %1th and 5th sleep, 2-3th task
t1 = 19*60+44; %first -> NEUROSCOPE
t2 = t1+33*60+19; %second
% kk1 = spiket(spiket<t1*SR); kk1ind = spikeind(spiket<t1*SR);
% kk1b = spiket(spiket>t4*SR); kk1indb = spikeind(spiket>t4*SR);
% kk2 = spiket((spiket>SR*t1)&(spiket<SR*t4));kk2ind = spikeind((spiket>SR*t1)&(spiket<SR*t4));
% spiketsleep = cat(1,kk1,kk1b); spikeindsleep = cat(1,kk1ind,kk1indb);
% spiketrun = kk2; spikeindrun = kk2ind;
% clear kk1 kk1b kk1ind kk1indb kk2 kk2ind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Separate neurons groups, FR and spikes during task and home-cage
    %I'll take both edges plus peak
    ripsCA1edges = evtca1.ripsCA1.edges; ripsCA1peak = evtca1.ripsCA1.trigs;
    ripsCA2edges = evtca2.ripsCA2.edges; ripsCA2peak = evtca2.ripsCA2.trigs;
    
    %separate for just sleep
    ripsCA1edgesSLEEP = ripsCA1edges((ripsCA1edges(:,1)<t1*Fs),:);
    ripsCA1peakSLEEP = ripsCA1peak((ripsCA1peak(:,1)<t1*Fs),:);
    ripsCA2edgesSLEEP = ripsCA2edges((ripsCA2edges(:,1)<t1*Fs),:);
    ripsCA2peakSLEEP = ripsCA2peak((ripsCA2peak(:,1)<t1*Fs),:);
    
    
    %separate for just task
%   ripsCA1task = ripsCA1((ripsCA1(:,1)>t1*Fs & ripsCA1(:,1)<t4*Fs),:);
%   ripsCA2task = ripsCA2((ripsCA2(:,1)>t1*Fs & ripsCA2(:,1)<t4*Fs),:);
%   ripsCA3task = ripsCA3((ripsCA3(:,1)>t1*Fs & ripsCA3(:,1)<t4*Fs),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% trigs in peak
% groupcells = ca3pyrlayerpyrcell;
% grouprips = ripsCA2sleep;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Rip Modulation for trig in CA1
% ripmat = []; 
% count = 0;
% 
% %param
% wide = 0.2*Fs; %time for hist in time samp
% stepbin = floor(5e-3*Fs); %time for each bin in time samp
% xlim = [-250:stepbin:250+stepbin]; %5ms
% xlimtic = [-200,-100,0,100,200];
% 
% tic
% for r = 1:length(grouprips)
%     t = grouprips(r);
%     binrange = [(t-wide):stepbin:+(t+wide+stepbin)];
%     for i = 1:length(groupcells) 
%         j = groupcells(i);
%         spikes = floor((spiketsleep(find(spikeindsleep==j))/SR)*Fs);
%         ripspk = spikes(find(spikes>(t-wide) & spikes<(t+wide)));
%         if ~isempty(ripspk)
% %             pause
%             count = count+1;
% %             ripmat2save{1,count} = ripspk;
% %             ripmat2save{2,count} = t;
%             ripmat(count,:) = histc(ripspk,binrange)./sum(histc(ripspk,binrange));
%         end
%         clear ripspk
%     end
%    clear binrange t
% end
% clear count ripspk binrange t
% toc
% 
% 
% 
% % 
% % ripCA2cellsCA2 = ripmat2save;
% % clear ripmat2save
% %% statistic hist
% 
% 
% %averageging
% for l = 1:size(ripmat,2)
% %     ripmat(l,:) = histc(ripCA1cellsCA1{1,l}(:),binrange)./sum(histc(ripCA1cellsCA1{1,l}(:),binrange));
%      ripmattotal(l) = sum(ripmat(:,l)) / size(ripmat,1);
% end
% clear l
% 
% ripca2sleeppyrca3 = ripmattotal;
% 
% 
% %plot
% figure
% bar(xlim,ripmattotal*Fs/stepbin)
% xlim([1 length(xlim)]);
% set(gca,'XTick',[xlim(1) xlim(22) xlim(42) xlim(64) xlim(end-1)]);
% set(gca,'XTickLabel',[xlimtic]);
% xlabel('Time (ms)');
% ylabel('Firing Rate (Hz)');
% x1=-stepbin;
% y1=get(gca,'ylim');
% hold on
% plot([x1 x1],y1,'r','LineWidth',2)
% title(['CA3 pyr modulation with rip trig-CA2 SLEEP']);
% 
% clear ripmattotal ripmat
% 
% %% save
% 
% save('D:\aza\analysis\data\AB4\01032013\general\rip\toshow\pyr\ripca1pyrca1.mat','ripca1pyrca1');
% save('D:\aza\analysis\data\AB4\01032013\general\rip\toshow\pyr\ripca1pyrca2.mat','ripca1pyrca2');
% save('D:\aza\analysis\data\AB4\01032013\general\rip\toshow\pyr\ripca1pyrca3.mat','ripca1pyrca3');
% 
% save('D:\aza\analysis\data\AB4\01032013\general\rip\toshow\pyr\ripca2pyrca1.mat','ripca2pyrca1');
% save('D:\aza\analysis\data\AB4\01032013\general\rip\toshow\pyr\ripca2pyrca2.mat','ripca2pyrca2');
% save('D:\aza\analysis\data\AB4\01032013\general\rip\toshow\pyr\ripca2pyrca3.mat','ripca2pyrca3');

%% trigs in window time
%groupcells = ca3pyrlayerintcell;
count=0;
totalbins = 84; %+bin 0
xlim = [-2.5:5/totalbins:2.5]; %5ms
stepbin = length(xlim)/5;
xlimtic = [-3,-2,-1,1,2,3];
Fs=1250;

grouprips = ripsCA2edgesSLEEP;
tic
for r = 1:size(grouprips,1)
    disp(['rip ' num2str(r)]);
    t1 = grouprips(r,1); t2 = grouprips(r,2); tiempo = t2-t1;
    binini=t1-2*tiempo;binfin=t2+2*tiempo;scaletiempo(r) = ((binfin-binini)/totalbins);
    binrange = [binini:((binfin-binini)/totalbins):binfin]; %around 4ms
    pre_window_t = t1-tiempo; post_window_t = t2 + tiempo;
    
    %ca1
    groupcells = ca1pyrlayerintcell;
    for i = 1:length(groupcells) 
        j = groupcells(i);
        spikes = floor((spiket(find(spikeind==j))/SR)*Fs);
        ripspk = spikes(find(spikes>binini & spikes<binfin));
        ripspkprewin = spikes(find(spikes>pre_window_t & spikes<t1));
        ripspkpostwin = spikes(find(spikes>t2 & spikes<post_window_t));
        ripspkinside = spikes(find(spikes>t1 & spikes<t2));
        %prob of spikes in diff windows..
        if ~isempty(ripspkinside) %inside rip window
           prob{1,1}{1,2}(i,r) = 1;
        else
           prob{1,1}{1,2}(i,r) = 0;
        end
        if ~isempty(ripspkprewin) %with 1 window plus
            prob{1,1}{1,1}(i,r) = 1;
        else
            prob{1,1}{1,1}(i,r) = 0;
        end
        if ~isempty(ripspkpostwin) %with 1 window plus
            prob{1,1}{1,3}(i,r) = 1;
        else
            prob{1,1}{1,3}(i,r) = 0;
        end
        %firing rate in a temporal wind around ripp event
        ripmatca1{1,1}(i,1) = j;
        if ~isempty(ripspk) %with 2 window plus
            prob{1,1}{1,4}(i,r) = 1;
            ripmatca1{1,2}{1,r}(i,:) = histc(ripspk,binrange)/sum(histc(ripspk,binrange));
            ripmatca1{1,3}{1,r}(i,:) = histc(ripspk,binrange);
        else
            ripmatca1{1,2}{1,r}(i,:) = zeros(1,length(xlim));
            ripmatca1{1,3}{1,r}(i,:) = zeros(1,length(xlim));
            prob{1,1}{1,4}(i,r) = 0;
        end
        clear ripspkprewin ripspkpostwin ripspkinside ripspk
    end
   clear i groupcells

    %ca2
    groupcells = ca2pyrlayerintcell;
    for i = 1:length(groupcells) 
        j = groupcells(i);
        spikes = floor((spiket(find(spikeind==j))/SR)*Fs);
        ripspk = spikes(find(spikes>binini & spikes<binfin));
        ripspkprewin = spikes(find(spikes>pre_window_t & spikes<t1));
        ripspkpostwin = spikes(find(spikes>t2 & spikes<post_window_t));
        ripspkinside = spikes(find(spikes>t1 & spikes<t2));
        %prob of spikes in diff windows..
        if ~isempty(ripspkinside) %inside rip window
           prob{1,2}{1,2}(i,r) = 1;
        else
           prob{1,2}{1,2}(i,r) = 0;
        end
        if ~isempty(ripspkprewin) %with 1 window plus
            prob{1,2}{1,1}(i,r) = 1;
        else
            prob{1,2}{1,1}(i,r) = 0;
        end
        if ~isempty(ripspkpostwin) %with 1 window plus
            prob{1,2}{1,3}(i,r) = 1;
        else
            prob{1,2}{1,3}(i,r) = 0;
        end
        %firing rate in a temporal wind around ripp event
        ripmatca2{1,1}(i,1) = j;
        if ~isempty(ripspk) %with 2 window plus
%             count = count+1;
            prob{1,2}{1,4}(i,r) = 1;
            ripmatca2{1,2}{1,r}(i,:) = histc(ripspk,binrange)./sum(histc(ripspk,binrange));
            ripmatca2{1,3}{1,r}(i,:) = histc(ripspk,binrange);
        else
            ripmatca2{1,2}{1,r}(i,:) = zeros(1,length(xlim));
            ripmatca2{1,3}{1,r}(i,:) = zeros(1,length(xlim));
            prob{1,2}{1,4}(i,r) = 0;
        end
        clear ripspkprewin ripspkpostwin ripspkinside ripspk
    end
    clear i groupcells
    
    %ca3
    groupcells = ca3pyrlayerintcell;
    for i = 1:length(groupcells) 
        j = groupcells(i);
        spikes = floor((spiket(find(spikeind==j))/SR)*Fs);
        ripspk = spikes(find(spikes>binini & spikes<binfin));
        ripspkprewin = spikes(find(spikes>pre_window_t & spikes<t1));
        ripspkpostwin = spikes(find(spikes>t2 & spikes<post_window_t));
        ripspkinside = spikes(find(spikes>t1 & spikes<t2));
        %prob of spikes in diff windows..
        if ~isempty(ripspkinside) %inside rip window
           prob{1,3}{1,2}(i,r) = 1;
        else
           prob{1,3}{1,2}(i,r) = 0;
        end
        if ~isempty(ripspkprewin) %with 1 window plus
            prob{1,3}{1,1}(i,r) = 1;
        else
            prob{1,3}{1,1}(i,r) = 0;
        end
        if ~isempty(ripspkpostwin) %with 1 window plus
            prob{1,3}{1,3}(i,r) = 1;
        else
            prob{1,3}{1,3}(i,r) = 0;
        end
        %firing rate in a temporal wind around ripp event
        ripmatca3{1,1}(i,1) = j;
        if ~isempty(ripspk) %with 2 window plus
%             count = count+1;
            prob{1,3}{1,4}(i,r) = 1;
%             contadorca1(i,r) = 1; % esto nunca puede ser mayor que el numero de cells select            
            ripmatca3{1,2}{1,r}(i,:) = histc(ripspk,binrange)./sum(histc(ripspk,binrange));
            ripmatca3{1,3}{1,r}(i,:) = histc(ripspk,binrange);
        else
            ripmatca3{1,2}{1,r}(i,:) = zeros(1,length(xlim));
            ripmatca3{1,3}{1,r}(i,:) = zeros(1,length(xlim));
            prob{1,3}{1,4}(i,r) = 0;
        end
        clear ripspkprewin ripspkpostwin ripspkinside ripspk
    end
   clear i groupcells    
   
    %dg
    groupcells = dginh;
    for i = 1:length(groupcells) 
        j = groupcells(i);
        spikes = floor((spiket(find(spikeind==j))/SR)*Fs);
        ripspk = spikes(find(spikes>binini & spikes<binfin));
        ripspkprewin = spikes(find(spikes>pre_window_t & spikes<t1));
        ripspkpostwin = spikes(find(spikes>t2 & spikes<post_window_t));
        ripspkinside = spikes(find(spikes>t1 & spikes<t2));
        %prob of spikes in diff windows..
        if ~isempty(ripspkinside) %inside rip window
           prob{1,4}{1,2}(i,r) = 1;
        else
           prob{1,4}{1,2}(i,r) = 0;
        end
        if ~isempty(ripspkprewin) %with 1 window plus
            prob{1,4}{1,1}(i,r) = 1;
        else
            prob{1,4}{1,1}(i,r) = 0;
        end
        if ~isempty(ripspkpostwin) %with 1 window plus
            prob{1,4}{1,3}(i,r) = 1;
        else
            prob{1,4}{1,3}(i,r) = 0;
        end
        %firing rate in a temporal wind around ripp event
        ripmatdg{1,1}(i,1) = j;
        if ~isempty(ripspk) %with 2 window plus
%             count = count+1;
            prob{1,4}{1,4}(i,r) = 1;
            ripmatdg{1,2}{1,r}(i,:) = histc(ripspk,binrange)./sum(histc(ripspk,binrange));
            ripmatdg{1,3}{1,r}(i,:) = histc(ripspk,binrange);
        else
            ripmatdg{1,2}{1,r}(i,:) = zeros(1,length(xlim));
            ripmatdg{1,3}{1,r}(i,:) = zeros(1,length(xlim));
            prob{1,4}{1,4}(i,r) = 0;
        end
        clear ripspkprewin ripspkpostwin ripspkinside ripspk
    end
    clear i groupcells
end
clear count r
toc

%save in vars
mat{1,1} = ripmatca1;
mat{1,2} = ripmatca2;
mat{1,3} = ripmatca3;
mat{1,4} = ripmatdg;
clear ripmatca1 ripmatca2 ripmatca3 ripmatdg

for r = 1:size(prob,2)
    for j = 1:size(prob{1,r},2)
        for i = 1:size(prob{1,r}{j},1) %neurons!
            prob2plot{1,r}{1,j}(i,1) = sum(prob{1,r}{j}(i,:))/size(prob{1,r}{j},2);
        end
        clear i
    end
    clear j
end
clear r

matca2 = mat; clear mat;
probca2 = prob; clear prob;
prob2plotca2 = prob2plot; clear prob2plot;

% %calculo probabilidad

%% plot...

%smooth
Smooth_Factor = 20;
Sigma = 1/Smooth_Factor; 
r = (-length(xlim):length(xlim))/length(xlim);
Smoother = exp(-r.^2/Sigma^2/2);
Smoother = Smoother/(sum(Smoother));

%averageging
toav = matca2;
for quien = 1:size(toav,2)
    ripmattotal{quien} = zeros(length(xlim),1);
    for que = 3 %2=normalizado 3=sin normalizar
        for r = 1:size(toav{1,quien}{1,que},2)
            for l = 1:size(toav{1,quien}{1,que}{1,r},2)%hist normalizado
                ripmattotal{quien}(l,1) = ripmattotal{quien}(l,1) + sum(toav{1,quien}{1,que}{1,r}(:,l)) / size(toav{1,quien}{1,que}{1,r},1);
            end
            clear l
        end
        ripmattotal{quien} = ripmattotal{quien}./size(toav{1,quien}{1,que},2);
    end
end
clear quien que r toav
ca2ripfiring=ripmattotal;clear ripmattotal

toanalyze=ca2ripfiring;
s1=conv2(Smoother, Smoother, toanalyze{1}, 'same');
s2=conv2(Smoother, Smoother, toanalyze{2}, 'same');
s3=conv2(Smoother, Smoother, toanalyze{3}, 'same');
s4=conv2(Smoother, Smoother, toanalyze{4}, 'same');
for i = 1:2
    FR(i) = max(toanalyze{i}*Fs);
end
clear i toanalyze
% % 
save(['D:\aza\analysis\data\AFR2\210214\general\rip\probability\wintrig\int\matca2.mat'],'matca2');
save(['D:\aza\analysis\data\AFR2\210214\general\rip\probability\wintrig\int\probca2.mat'],'probca2');
save(['D:\aza\analysis\data\AFR2\210214\general\rip\probability\wintrig\int\prob2plotca2.mat'],'prob2plotca2');

save(['D:\aza\analysis\data\AFR2\210214\general\rip\probability\wintrig\int\ca2ripca1firing' num2str(Smooth_Factor) '.mat'],'s1');
save(['D:\aza\analysis\data\AFR2\210214\general\rip\probability\wintrig\int\ca2ripca2firing' num2str(Smooth_Factor) '.mat'],'s2');

ticca1=max(s1);yticca1=find(s1==ticca1);
ticca2=max(s2);yticca2=find(s2==ticca2);
ticca3=max(s3);yticca3=find(s3==ticca3);
ticca4=max(s4);yticca4=find(s4==ticca4);
%nota al plotear el eje count sacarlo con el histograma!
figure
plot(xlim,s1,'color',[1,0.75,0.3],'LineWidth',3); hold on
plot(xlim,s2,'color',[0.4,0.6,0.9],'LineWidth',3); hold on
plot(xlim,s3,'color',[0.7,0,0.7],'LineWidth',3); hold on
plot(xlim,s4,'r','LineWidth',3); hold on
quienmax = max([ticca1,ticca2,ticca3,ticca4]);
xlim([1 length(xlim)]); %ylim([0 quienmax+quienmax*0.4]) %y axis missing after convolution
set(gca,'XTick',[xlim(1) xlim(stepbin*1) xlim(stepbin*2) xlim(stepbin*3) xlim(stepbin*4) xlim(end)]);
set(gca,'XTickLabel',[xlimtic]);
xlabel('ripple window'); ylabel('Firing Rate (Hz)');
x1=-xlim(stepbin*2);x2=+xlim(stepbin*2);
y1=get(gca,'ylim');y2=y1; hold on;
plot([x1 x1],y1,'--r','LineWidth',2);hold on;
plot([x2 x2],y2,'--r','LineWidth',2);hold on;
% ylim([0 0.01]);
plot([xlim(yticca1) xlim(yticca1)],get(gca,'ylim')/10+quienmax+quienmax/5,'color',[1,0.75,0.3],'LineWidth',3);hold on;
plot([xlim(yticca2) xlim(yticca2)],get(gca,'ylim')/10+quienmax+quienmax/5,'color',[0.4,0.6,0.9],'LineWidth',3);hold on;
plot([xlim(yticca3) xlim(yticca3)],get(gca,'ylim')/10+quienmax+quienmax/5,'color',[0.7,0,0.7],'LineWidth',3);hold on;
plot([xlim(yticca4) xlim(yticca4)],get(gca,'ylim')/10+quienmax+quienmax/5,'r','LineWidth',3);hold on;
hold on
suptitle(['CA3ripp, int cells, FR ' num2str(FR)]);
legend('ca1','ca2','ca3','dg')


%% probability plots

figure(2)
subplot(1,3,1)
xy (:,1) = prob2plotca1{1,1}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,1}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,2}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,3}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,4}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA3 rip');
title('pre-rip window');
legend('ca1','ca2','ca3','dg');
subplot(1,3,2)
xy (:,1) = prob2plotca1{1,1}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,1}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,2}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,3}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,4}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA3 rip');
title('rip window');
legend('ca1','ca2','ca3','dg');
subplot(1,3,3)
clear xy
xy (:,1) = prob2plotca1{1,1}{1,3}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,1}{1,3}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,3}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,2}{1,3}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,3}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,3}{1,3}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,3}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,4}{1,3}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA3 rip');
title('post-rip window');
legend('ca1','ca2','ca3','dg');
suptitle('ca1-ca3')

figure(3)
subplot(1,3,1)
xy (:,1) = prob2plotca1{1,1}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);
x1=0;y1=0;hold on;
plot([x1 x2],[y1 y2],'--r','LineWidth',2);hold on;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA2 rip');
title('pre-rip window');
legend('ca1','ca2','ca3','dg');
subplot(1,3,2)
clear xy
xy (:,1) = prob2plotca1{1,1}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);
x1=0;y1=0;hold on;
plot([x1 x2],[y1 y2],'--r','LineWidth',2);hold on;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA2 rip');
title('rip window');
legend('ca1','ca2','ca3','dg');
subplot(1,3,3)
xy (:,1) = prob2plotca1{1,1}{1,3}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,3}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,3}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,3}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,3}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,3}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,3}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,3}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);
x1=0;y1=0;hold on;
plot([x1 x2],[y1 y2],'--r','LineWidth',2);hold on;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA2 rip');
title('post-rip window');
legend('ca1','ca2','ca3','dg');
suptitle('ca1-ca2')

figure(4)
subplot(1,3,1)
xy (:,1) = prob2plotca3{1,1}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca3{1,2}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca3{1,3}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca3{1,4}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);
x1=0;y1=0;hold on;
plot([x1 x2],[y1 y2],'--r','LineWidth',2);hold on;
xlabel('Firing probability during CA3 rip');
ylabel('Firing probability during CA2 rip');
title('pre-rip window');
legend('ca1','ca2','ca3','dg');
subplot(1,3,2)
clear xy
xy (:,1) = prob2plotca3{1,1}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca3{1,2}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca3{1,3}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca3{1,4}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);
x1=0;y1=0;hold on;
plot([x1 x2],[y1 y2],'--r','LineWidth',2);hold on;
xlabel('Firing probability during CA3 rip');
ylabel('Firing probability during CA2 rip');
title('rip window');
legend('ca1','ca2','ca3','dg');
subplot(1,3,3)
clear xy
xy (:,1) = prob2plotca3{1,1}{1,3}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,3}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca3{1,2}{1,3}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,3}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca3{1,3}{1,3}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,3}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca3{1,4}{1,3}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,3}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);
x1=0;y1=0;hold on;
plot([x1 x2],[y1 y2],'--r','LineWidth',2);hold on;
xlabel('Firing probability during CA3 rip');
ylabel('Firing probability during CA2 rip');
title('post-rip window');
legend('ca1','ca2','ca3','dg');
suptitle('ca3-ca2')

%aqui pinto toda la ventana que considero during ripp time plus and minus
%si puedo separarla en 5 trozos y pintar los dos primeros el siguiente y
%los dos ultimos tengo lo anterior pero alargado en tiempo, quiza...
figure(5)
subplot(1,3,1)
xy (:,1) = prob2plotca1{1,1}{1,4}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,1}{1,4}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,4}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,2}{1,4}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,4}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,3}{1,4}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,4}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,4}{1,4}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);
x1=0;y1=0;hold on;
plot([x1 x2],[y1 y2],'--r','LineWidth',2);hold on;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA3 rip');
title('ca1-ca3');
legend('ca1','ca2','ca3','dg');
subplot(1,3,2)
xy (:,1) = prob2plotca1{1,1}{1,4}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,4}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,4}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,4}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,4}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,4}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,4}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,4}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);
x1=0;y1=0;hold on;
plot([x1 x2],[y1 y2],'--r','LineWidth',2);hold on;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA2 rip');
title('ca1-ca2');
legend('ca1','ca2','ca3','dg');
subplot(1,3,3)
xy (:,1) = prob2plotca3{1,1}{1,4}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,4}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca3{1,2}{1,4}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,4}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca3{1,3}{1,4}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,4}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca3{1,4}{1,4}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,4}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);
x1=0;y1=0;hold on;
plot([x1 x2],[y1 y2],'--r','LineWidth',2);hold on;
xlabel('Firing probability during CA3 rip');
ylabel('Firing probability during CA2 rip');
title('ca3-ca2');
legend('ca1','ca2','ca3','dg');


%%

%% trigs in ripp peak

%groupcells = ca3pyrlayerintcell;
count=0;
Fs=1250;
wide = 0.2*Fs; %time for hist in time samp
stepbin = floor(4e-3*Fs); %time for each bin in time samp
xlim = [-250:stepbin:250]; %5ms
xlimtic = [-200,-100,0,100,200];
tiempo1 = 60*1e-3*Fs; tiempo2 = 120*1e-3*Fs;

grouprips = ripsCA2peakSLEEP;
tic
for r = 1:size(grouprips,1)
    disp(['rip ' num2str(r)]);
    t = grouprips(r);
    binrange = [(t-wide):stepbin:+(t+wide)];
    pre_window_t1 = t-tiempo1; post_window_t1 = t + tiempo1;
    pre_window_t2 = t-tiempo2; post_window_t2 = t + tiempo2;
    
    %ca1
    groupcells = ca1pyrlayerintcell;
    for i = 1:length(groupcells) 
        j = groupcells(i);
        spikes = floor((spiket(find(spikeind==j))/SR)*Fs);
        ripspk = spikes(find(spikes>(t-wide) & spikes<(t+wide)));
        ripspk60 = spikes(find(spikes>(t-60*Fs*1e-3) & spikes<(t+60*Fs*1e-3)));                        
        ripspkprewin1 = spikes(find(spikes>pre_window_t1 & spikes<t));
        ripspkpostwin1 = spikes(find(spikes>t & spikes<post_window_t1));
        ripspkprewin2 = spikes(find(spikes>pre_window_t2 & spikes<t));
        ripspkpostwin2 = spikes(find(spikes>t & spikes<post_window_t2));
        
        %prob of spikes in diff windows..
        if ~isempty(ripspkprewin2) 
            prob{1,1}{1,1}(i,r) = 1;
        else
            prob{1,1}{1,1}(i,r) = 0;
        end
        if ~isempty(ripspkprewin1) %with 1 window plus
            prob{1,1}{1,2}(i,r) = 1;
        else
            prob{1,1}{1,2}(i,r) = 0;
        end
        if ~isempty(ripspkpostwin1) %with 1 window plus
            prob{1,1}{1,3}(i,r) = 1;
        else
            prob{1,1}{1,3}(i,r) = 0;
        end
        if ~isempty(ripspkpostwin2) %with 1 window plus
            prob{1,1}{1,4}(i,r) = 1;
        else
            prob{1,1}{1,4}(i,r) = 0;
        end
        if ~isempty(ripspk60) %with 1 window plus
            prob{1,1}{1,5}(i,r) = 1;
        else
            prob{1,1}{1,5}(i,r) = 0;
        end
        %firing rate in a temporal wind around ripp event
        ripmatca1{1,1}(i,1) = j;
        if ~isempty(ripspk) %with 2 window plus
            prob{1,1}{1,6}(i,r) = 1;
            ripmatca1{1,2}{1,r}(i,:) = histc(ripspk,binrange)/sum(histc(ripspk,binrange));
            ripmatca1{1,3}{1,r}(i,:) = histc(ripspk,binrange);
        else
            ripmatca1{1,2}{1,r}(i,:) = zeros(1,length(xlim));
            ripmatca1{1,3}{1,r}(i,:) = zeros(1,length(xlim));
            prob{1,1}{1,6}(i,r) = 0;
        end
        clear ripspkprewin ripspkpostwin ripspkinside ripspk
    end
   clear i groupcells

    %ca2
    groupcells = ca2pyrlayerintcell;
    for i = 1:length(groupcells) 
        j = groupcells(i);
        spikes = floor((spiket(find(spikeind==j))/SR)*Fs);
        ripspk = spikes(find(spikes>(t-wide) & spikes<(t+wide)));
        ripspk60 = spikes(find(spikes>(t-60*Fs*1e-3) & spikes<(t+60*Fs*1e-3)));                
        ripspkprewin1 = spikes(find(spikes>pre_window_t1 & spikes<t));
        ripspkpostwin1 = spikes(find(spikes>t & spikes<post_window_t1));
        ripspkprewin2 = spikes(find(spikes>pre_window_t2 & spikes<t));
        ripspkpostwin2 = spikes(find(spikes>t & spikes<post_window_t2));
        
        %prob of spikes in diff windows..
        if ~isempty(ripspkprewin2) %with 1 window plus
            prob{1,2}{1,1}(i,r) = 1;
        else
            prob{1,2}{1,1}(i,r) = 0;
        end        
        if ~isempty(ripspkprewin1) %with 1 window plus
            prob{1,2}{1,2}(i,r) = 1;
        else
            prob{1,2}{1,2}(i,r) = 0;
        end
        if ~isempty(ripspkpostwin1) %with 1 window plus
            prob{1,2}{1,3}(i,r) = 1;
        else
            prob{1,2}{1,3}(i,r) = 0;
        end
        if ~isempty(ripspkpostwin2) %with 1 window plus
            prob{1,2}{1,4}(i,r) = 1;
        else
            prob{1,2}{1,4}(i,r) = 0;
        end
        if ~isempty(ripspkpostwin2) %with 1 window plus
            prob{1,2}{1,5}(i,r) = 1;
        else
            prob{1,2}{1,5}(i,r) = 0;
        %firing rate in a temporal wind around ripp event
        end
        ripmatca2{1,1}(i,1) = j;
        if ~isempty(ripspk) %with 2 window plus
            prob{1,2}{1,6}(i,r) = 1;
            ripmatca2{1,2}{1,r}(i,:) = histc(ripspk,binrange)./sum(histc(ripspk,binrange));
            ripmatca2{1,3}{1,r}(i,:) = histc(ripspk,binrange);
        else
            ripmatca2{1,2}{1,r}(i,:) = zeros(1,length(xlim));
            ripmatca2{1,3}{1,r}(i,:) = zeros(1,length(xlim));
            prob{1,2}{1,6}(i,r) = 0;
        end
        clear ripspkprewin ripspkpostwin ripspkinside ripspk
    end
    clear i groupcells
    
    %ca3
    groupcells = ca3pyrlayerintcell;
    for i = 1:length(groupcells) 
        j = groupcells(i);
       spikes = floor((spiket(find(spikeind==j))/SR)*Fs);
        ripspk = spikes(find(spikes>(t-wide) & spikes<(t+wide)));
        ripspk60 = spikes(find(spikes>(t-60*Fs*1e-3) & spikes<(t+60*Fs*1e-3)));        
        ripspkprewin1 = spikes(find(spikes>pre_window_t1 & spikes<t));
        ripspkpostwin1 = spikes(find(spikes>t & spikes<post_window_t1));
        ripspkprewin2 = spikes(find(spikes>pre_window_t2 & spikes<t));
        ripspkpostwin2 = spikes(find(spikes>t & spikes<post_window_t2));
        
        %prob of spikes in diff windows..
        if ~isempty(ripspkprewin2) %with 1 window plus
            prob{1,3}{1,1}(i,r) = 1;
        else
            prob{1,3}{1,1}(i,r) = 0;
        end        
        if ~isempty(ripspkprewin1) %with 1 window plus
            prob{1,3}{1,2}(i,r) = 1;
        else
            prob{1,3}{1,2}(i,r) = 0;
        end
        if ~isempty(ripspkpostwin1) %with 1 window plus
            prob{1,3}{1,3}(i,r) = 1;
        else
            prob{1,3}{1,3}(i,r) = 0;
        end
        if ~isempty(ripspkpostwin2) %with 1 window plus
            prob{1,3}{1,4}(i,r) = 1;
        else
            prob{1,3}{1,4}(i,r) = 0;
        end
        if ~isempty(ripspk60) %with 1 window plus
            prob{1,3}{1,5}(i,r) = 1;
        else
            prob{1,3}{1,5}(i,r) = 0;
        end
        %firing rate in a temporal wind around ripp event
        ripmatca3{1,1}(i,1) = j;
        if ~isempty(ripspk) %with 2 window plus
            prob{1,3}{1,6}(i,r) = 1;
            ripmatca3{1,2}{1,r}(i,:) = histc(ripspk,binrange)./sum(histc(ripspk,binrange));
            ripmatca3{1,3}{1,r}(i,:) = histc(ripspk,binrange);
        else
            ripmatca3{1,2}{1,r}(i,:) = zeros(1,length(xlim));
            ripmatca3{1,3}{1,r}(i,:) = zeros(1,length(xlim));
            prob{1,3}{1,6}(i,r) = 0;
        end
        clear ripspkprewin ripspkpostwin ripspkinside ripspk
    end
   clear i groupcells    
   
    %dg
    groupcells = dginh;
    for i = 1:length(groupcells) 
        j = groupcells(i);
       spikes = floor((spiket(find(spikeind==j))/SR)*Fs);
        ripspk = spikes(find(spikes>(t-wide) & spikes<(t+wide)));
        ripspk60 = spikes(find(spikes>(t-60*Fs*1e-3) & spikes<(t+60*Fs*1e-3)));
        ripspkprewin1 = spikes(find(spikes>pre_window_t1 & spikes<t));
        ripspkpostwin1 = spikes(find(spikes>t & spikes<post_window_t1));
        ripspkprewin2 = spikes(find(spikes>pre_window_t2 & spikes<t));
        ripspkpostwin2 = spikes(find(spikes>t & spikes<post_window_t2));
        
        %prob of spikes in diff windows..
        if ~isempty(ripspkprewin2) %with 1 window plus
            prob{1,4}{1,1}(i,r) = 1;
        else
            prob{1,4}{1,1}(i,r) = 0;
        end        
        if ~isempty(ripspkprewin1) %with 1 window plus
            prob{1,4}{1,2}(i,r) = 1;
        else
            prob{1,4}{1,2}(i,r) = 0;
        end
        if ~isempty(ripspkpostwin1) %with 1 window plus
            prob{1,4}{1,3}(i,r) = 1;
        else
            prob{1,4}{1,3}(i,r) = 0;
        end
        if ~isempty(ripspkpostwin2) %with 1 window plus
            prob{1,4}{1,4}(i,r) = 1;
        else
            prob{1,4}{1,4}(i,r) = 0;
        end
        if ~isempty(ripspk60) %with 1 window plus
            prob{1,4}{1,5}(i,r) = 1;
        else
            prob{1,4}{1,5}(i,r) = 0;
        end
        %firing rate in a temporal wind around ripp event
        ripmatdg{1,1}(i,1) = j;
        if ~isempty(ripspk) %with 2 window plus
            prob{1,4}{1,6}(i,r) = 1;
            ripmatdg{1,2}{1,r}(i,:) = histc(ripspk,binrange)./sum(histc(ripspk,binrange));
            ripmatdg{1,3}{1,r}(i,:) = histc(ripspk,binrange);
        else
            ripmatdg{1,2}{1,r}(i,:) = zeros(1,length(xlim));
            ripmatdg{1,3}{1,r}(i,:) = zeros(1,length(xlim));
            prob{1,4}{1,6}(i,r) = 0;
        end
        clear ripspkprewin ripspkpostwin ripspkinside ripspk
    end
    clear i groupcells
end
clear count r
toc

%save in vars
mat{1,1} = ripmatca1;
mat{1,2} = ripmatca2;
mat{1,3} = ripmatca3;
mat{1,4} = ripmatdg;
clear ripmatca1 ripmatca2 ripmatca3 ripmatdg

for r = 1:size(prob,2)
    for j = 1:size(prob{1,r},2)
        for i = 1:size(prob{1,r}{j},1) %neurons!
            prob2plot{1,r}{1,j}(i,1) = sum(prob{1,r}{j}(i,:))/size(prob{1,r}{j},2);
        end
        clear i
    end
    clear j
end
clear r

matca2 = mat; clear mat;
probca2 = prob; clear prob;
prob2plotca2 = prob2plot; clear prob2plot;

%% plot...

%smooth
Smooth_Factor = 20;
Sigma = 1/Smooth_Factor; 
r = (-length(xlim):length(xlim))/length(xlim);
Smoother = exp(-r.^2/Sigma^2/2);
Smoother = Smoother/(sum(Smoother));

%averageging
toav = matca2;
for quien = 1:size(toav,2)
    ripmattotal{quien} = zeros(length(xlim),1);
    for que = 3 %2=normalizado 3=sin normalizar
        for r = 1:size(toav{1,quien}{1,que},2)
            for l = 1:size(toav{1,quien}{1,que}{1,r},2)%hist normalizado
                ripmattotal{quien}(l,1) = ripmattotal{quien}(l,1) + sum(toav{1,quien}{1,que}{1,r}(:,l)) / size(toav{1,quien}{1,que}{1,r},1);
            end
            clear l
        end
        ripmattotal{quien} = ripmattotal{quien}./size(toav{1,quien}{1,que},2);
    end
end
clear quien que r toav
ca2ripfiring=ripmattotal;clear ripmattotal

toanalyze=ca2ripfiring;
s1=conv2(Smoother, Smoother, toanalyze{1}, 'same');
s2=conv2(Smoother, Smoother, toanalyze{2}, 'same');
s3=conv2(Smoother, Smoother, toanalyze{3}, 'same');
s4=conv2(Smoother, Smoother, toanalyze{4}, 'same');
for i = 1:2
    FR(i) = max(toanalyze{i}*Fs);
end
clear i toanalyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(['D:\aza\analysis\data\AFR2\210214\general\rip\probability\peaktrig\int\matca2.mat'],'matca2');
save(['D:\aza\analysis\data\AFR2\210214\general\rip\probability\peaktrig\int\probca2.mat'],'probca2');
save(['D:\aza\analysis\data\AFR2\210214\general\rip\probability\peaktrig\int\prob2plotca2.mat'],'prob2plotca2');

save(['D:\aza\analysis\data\AFR2\210214\general\rip\probability\peaktrig\int\ca2ripca1firing' num2str(Smooth_Factor) '.mat'],'s1');
save(['D:\aza\analysis\data\AFR2\210214\general\rip\probability\peaktrig\int\ca2ripca2firing' num2str(Smooth_Factor) '.mat'],'s2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ticca1=max(s1);yticca1=find(s1==ticca1);
ticca2=max(s2);yticca2=find(s2==ticca2);
ticca3=max(s3);yticca3=find(s3==ticca3);
ticca4=max(s4);yticca4=find(s4==ticca4);
%nota al plotear el eje count sacarlo con el histograma!
figure
plot(xlim,s1,'color',[1,0.75,0.3],'LineWidth',3); hold on
plot(xlim,s2,'color',[0.4,0.6,0.9],'LineWidth',3); hold on
plot(xlim,s3,'color',[0.7,0,0.7],'LineWidth',3); hold on
plot(xlim,s4,'r','LineWidth',3); hold on
quienmax = max([ticca1,ticca2,ticca3,ticca4]);
xlim([1 length(xlim)]); %ylim([0 quienmax+quienmax*0.4]) %y axis missing after convolution
set(gca,'XTick',[xlim(1) xlim(26) xlim(51) xlim(76) xlim(end)]);
set(gca,'XTickLabel',[xlimtic]);
xlabel('ripple window'); ylabel('Firing Rate (Hz)');
x1=0;y1=get(gca,'ylim');hold on;
plot([x1 x1],y1,'--r','LineWidth',2);hold on;
% ylim([0 0.01]);
plot([xlim(yticca1) xlim(yticca1)],get(gca,'ylim')/10+quienmax+quienmax/5,'color',[1,0.75,0.3],'LineWidth',3);hold on;
plot([xlim(yticca2) xlim(yticca2)],get(gca,'ylim')/10+quienmax+quienmax/5,'color',[0.4,0.6,0.9],'LineWidth',3);hold on;
plot([xlim(yticca3) xlim(yticca3)],get(gca,'ylim')/10+quienmax+quienmax/5,'color',[0.7,0,0.7],'LineWidth',3);hold on;
plot([xlim(yticca4) xlim(yticca4)],get(gca,'ylim')/10+quienmax+quienmax/5,'r','LineWidth',3);hold on;
hold on
suptitle(['CA3 ripp, int cells, FR ' num2str(FR)]);
legend('ca1','ca2','ca3','dg')

%% prueba
clear xy
figure(2)
subplot(1,5,1)
xy (:,1) = prob2plotca1{1,1}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,1}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,2}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,3}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,4}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA3 rip');
title('pre-rip window 120');
legend('ca1','ca2','ca3','dg');
subplot(1,5,2)
xy (:,1) = prob2plotca1{1,1}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,1}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,2}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,3}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,4}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA3 rip');
title('rip window pre 60');
legend('ca1','ca2','ca3','dg');
subplot(1,5,3)
clear xy
xy (:,1) = prob2plotca1{1,1}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,1}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,2}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,3}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,4}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA3 rip');
title('rip window pre 60');
legend('ca1','ca2','ca3','dg');
subplot(1,5,4)
xy (:,1) = prob2plotca1{1,1}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,1}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,2}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,3}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,4}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA3 rip');
title('post-rip window 60ms');
legend('ca1','ca2','ca3','dg');
subplot(1,5,5)
xy (:,1) = prob2plotca1{1,1}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,1}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,2}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,3}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca3{1,4}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA3 rip');
title('two sides win 60ms');
legend('ca1','ca2','ca3','dg');
suptitle('Firing Probability')
clear xy

clear xy
figure(2)
subplot(1,5,1)
xy (:,1) = prob2plotca1{1,1}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA2 rip');
title('pre-rip window 120');
legend('ca1','ca2','ca3','dg');
subplot(1,5,2)
xy (:,1) = prob2plotca1{1,1}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA2 rip');
title('rip window pre 60');
legend('ca1','ca2','ca3','dg');
subplot(1,5,3)
xy (:,1) = prob2plotca1{1,1}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA2 rip');
title('rip window pre 60');
legend('ca1','ca2','ca3','dg');
subplot(1,5,4)
xy (:,1) = prob2plotca1{1,1}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA2 rip');
title('post-rip window 60ms');
legend('ca1','ca2','ca3','dg');
subplot(1,5,5)
xy (:,1) = prob2plotca1{1,1}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca1{1,2}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca1{1,3}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca1{1,4}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA1 rip');
ylabel('Firing probability during CA2 rip');
title('two sides win 60ms');
legend('ca1','ca2','ca3','dg');
suptitle('Firing Probability')
clear xy


clear xy
figure(4)
subplot(1,5,1)
xy (:,1) = prob2plotca3{1,1}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca3{1,2}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca3{1,3}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca3{1,4}{1,1}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,1}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA3 rip');
ylabel('Firing probability during CA2 rip');
title('pre-rip window 120');
legend('ca1','ca2','ca3','dg');
subplot(1,5,2)
xy (:,1) = prob2plotca3{1,1}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca3{1,2}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca3{1,3}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca3{1,4}{1,2}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,2}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA3 rip');
ylabel('Firing probability during CA2 rip');
title('rip window pre 60');
legend('ca1','ca2','ca3','dg');
subplot(1,5,3)
xy (:,1) = prob2plotca3{1,1}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca3{1,2}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca3{1,3}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca3{1,4}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA3 rip');
ylabel('Firing probability during CA2 rip');
title('rip window pre 60');
legend('ca1','ca2','ca3','dg');
subplot(1,5,4)
xy (:,1) = prob2plotca3{1,1}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca3{1,2}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca3{1,3}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca3{1,4}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA3 rip');
ylabel('Firing probability during CA2 rip');
title('post-rip window 60ms');
legend('ca1','ca2','ca3','dg');
subplot(1,5,5)
xy (:,1) = prob2plotca3{1,1}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,1}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'y','fill');hold on;
kk1 = max(xy);x2a=kk1(1);y2a=kk1(2);clear kk1 xy
xy (:,1) = prob2plotca3{1,2}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,2}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'b','fill');hold on;
kk2 = max(xy);x2b=kk2(1);y2b=kk2(2);clear kk2 xy
xy (:,1) = prob2plotca3{1,3}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,3}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'m','fill');hold on;
kk3 = max(xy);x2c=kk3(1);y2c=kk3(2);clear kk3 xy
xy (:,1) = prob2plotca3{1,4}{1,5}(:,1); %prob ca1 in ca1 rip time window
xy (:,2) = prob2plotca2{1,4}{1,5}(:,1); %prob ca1 in ca2 rip
scatter(xy(:,1),xy(:,2),'r','fill');
kk4 = max(xy);x2d=kk4(1);y2d=kk4(2);clear kk4 xy
x2=max([x2a,x2b,x2c,x2d]);y2=max([y2a,y2b,y2c,y2d]);bisec=max([x2,y2]);
x1=0;y1=0;hold on;
plot([x1 bisec],[y1 bisec],'--r','LineWidth',2);hold on;clear bisec;
xlabel('Firing probability during CA3 rip');
ylabel('Firing probability during CA2 rip');
title('two sides win 60ms');
legend('ca1','ca2','ca3','dg');
suptitle('Firing Probability')
clear xy




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%