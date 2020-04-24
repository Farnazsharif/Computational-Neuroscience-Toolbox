function [ph]=LFP_Phase_linear(lfp,PhaseFreqVector,PhaseFreq_BandWidth,srate)
% [phase]=cellph_Overal(filename,PhaseFreqVector,PhaseFreq_BandWidth,srate,lfp,spikeT) 
% PhaseFreqVector=[7];PhaseFreq_BandWidth=2;EEGsamplerate=1000;numchannel=64;
    Pf1 = PhaseFreqVector;
    Pf2 = Pf1 + PhaseFreq_BandWidth;
    PhaseFreq=eegfilt(lfp,srate,Pf1,Pf2); % this is just filtering 
    ph = angle(hilbert(PhaseFreq))';     
    ph = ph + cumsum([0;diff(ph)<0])*2*pi;
    ph= rem(ph,2*pi);
end

