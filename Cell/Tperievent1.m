%function [t,avgPSTH,spiketrial,trialPSTH]=Tperievent(spikeTime,eventTime,TimeWindow,binsize)
%
%t: time vector of peri-event histogram
%avgPSTH: average peri-event histogram (PSTH)
%spiketrial: 2 row matrix, row1 the spike times, row2 the trial# (for raster plots) 
%trialPSTH: A matrix where rows are individual trials PSTH
%
%spiteTime: spike time in s
%eventTime: time (in s) of stimulus events
%TimeWindow: [-2 3] for example set the time window from -2 to 3 s around the events
%binsize: size of bin (in s)
%
%Plot if no output specified

function [t,avgPSTH,spiketrial,trialPSTH,Ndxspk]=Tperievent1(spikeTime,eventTime,TimeWindow,binsize)
%%
% Srate=(spkinfo.samplerate./25);
% eventTime=Rt(:,79-64)./Srate;
% spikeTime=spikeT./spkinfo.samplerate;
% TimeWindow= [-0.1 0.1];
% binsize=0.02;
%%
t=TimeWindow(1):binsize:TimeWindow(2);

spiketrial=[];

trialPSTH=[];

Ndxspk=[];

for ii =1:length(eventTime)
    
    Tlimit = TimeWindow + eventTime(ii);
    Ndx=find(spikeTime>=Tlimit(1)-binsize/2&spikeTime<=Tlimit(2)+binsize/2);
    spkt=spikeTime(Ndx) - eventTime(ii);
    
    Ndxspk=[Ndxspk Ndx'];
    
    spiketrial=[spiketrial [spkt';ii*ones(size(spkt'))]];
	
	trialPSTH=[trialPSTH;hist(spkt,t)];
    
end

avgPSTH = mean(trialPSTH,1);

end