function [phase]=cellph_Overal(filename,PhaseFreqVector,PhaseFreq_BandWidth,srate,lfp,spikeT)     
%% Define the Amplitude- and Phase- Frequencies
load([filename '.mat']);
% PhaseFreqVector=PhaseFreqVector_a:PhaseFreqVector_b;
% PhaseFreq_BandWidth=4;
% data_length=length(lfp);


%%
% [spikeT,spikeG]=selectgroup(spk.i,spk.g,celln);
% period1=behav.restT; 
% [tt_rest]=selectt(spikeT,spikeG,period1,spkinfo.samplerate);
% period2=behav.runT; 
% [tt_run]=selectt(spikeT,spikeG,period2,spkinfo.samplerate);


%% Do filtering and Hilbert transform on CPU

'CPU filtering'

 
for jj=1:length(PhaseFreqVector)
    ph=[];
    Pf1 = PhaseFreqVector(jj);
    Pf2 = Pf1 + PhaseFreq_BandWidth;
    PhaseFreq=eegfilt(lfp,srate,Pf1,Pf2); % this is just filtering 
    ph = angle(hilbert(PhaseFreq))';
        
%    [phase]=spkphase1(spk.i,spk.g,celln,spkinfo.samplerate,srate,ph,spikeT);
    ph = ph + cumsum([0;diff(ph)<0])*2*pi;
    ph = interp1((1:length(ph))/srate,ph,spikeT);
    ph= rem(ph,2*pi);
   phase(:,jj)=ph;
end


end
