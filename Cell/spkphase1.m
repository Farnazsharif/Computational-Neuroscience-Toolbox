%function [phase,g]=spkphase(T,G,gsubset,smplrate,eeg,freqband,eegsmplrate)
%
%T: spk.i
%G
%gsubset: ex: G(1:10)
%smplrate: spkinfo.samplingrate
%eeg: if only 1 column, use that one for all spike phase calculation. If multiple, use eeg(:,1) for G==100:199, eeg(:,2) for G==200:299, etc...
%freqband: ex: [6 9] for theta
%eegsmplrate: ex: spkinfo.samplingrate/25 

function [phase]=spkphase1(T,G,gsubset,smplrate,eegsmplrate,ph,tt)


%%
% PhaseFreq=eegfilt(eeg',eegsmplrate,6,9); % this is just filtering  
% ph=angle(hilbert(PhaseFreq))'; % this is getting the phase time series
%% 
%  	ph=(phaseCont(eeg,freqband,eegsmplrate));
    
    ph = ph + cumsum([0;diff(ph)<0])*2*pi;%%

% 	[tt,g]=selectgroup(T,G,gsubset);

	ph = interp1((1:length(ph))/eegsmplrate,ph,tt/smplrate); % seb version, when you take the time from selectgroup, you need to divide to spkinfo.srate
%     ph = interp1((1:length(ph))/eegsmplrate,ph,tt); % in my version, it is already devided to srate
    
    phase = rem(ph,2*pi);


end