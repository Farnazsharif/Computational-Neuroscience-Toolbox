
function ripplemap(filename,cellN,Rt)
% ripplemap('FM05_1',903,Rt3)
% filename='FM05_1';
% cellN=903;
% Rt=Rt3;

load([filename '.mat']);
EEGsamplerate=(spkinfo.samplerate./25);
Srate=(spkinfo.samplerate./25); 
TimeWindow= [-0.1 0.1];
binsize=0.02;

[spikeT,spikeG] = selectgroup(spk.i,spk.g,cellN);
spikeTime=spikeT./spkinfo.samplerate;

%%

chorder=[1 8 2 7 3 6 4 5];
nsite=8;

reorder=[];
for ii = 1:16
    reorder = [reorder chorder+(ii-1)*nsite];
end
M=reorder(65:128);
W=reshape(M,8,8);

%% 
tic
for kk=1:8

for ll=1:8
ch=W(kk,ll);  
    

eventTime=Rt(:,ch-64)./EEGsamplerate;
subplot(8,8,ll+(kk-1)*8)
[t,avgPSTH,spiketrial,trialPSTH,Ndxspk]=Tperievent(spikeTime,eventTime,TimeWindow,binsize);

    plot([spiketrial(1,:);spiketrial(1,:)],length(eventTime)+1-[spiketrial(2,:);spiketrial(2,:)+0.9],'b')
    xlim(TimeWindow)
    ylim([0 length(eventTime)])
%     axis off
    
    hold on
    
    plot(t,avgPSTH/binsize,'k')
    ylabel('Hz')
    xlabel('s')
    xlim(TimeWindow)

end   
end

%%
set(gcf,'color','w')
dim = [.4 .7 .3 .3];
str = (['Cell#= ' num2str(cellN) ] );
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none');
toc