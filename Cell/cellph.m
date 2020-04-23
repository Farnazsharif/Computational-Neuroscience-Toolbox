function [MI,n,phase]=cellph(celln,ch,filename,PhaseFreqVector_a,PhaseFreqVector_b)
%% Define the Amplitude- and Phase- Frequencies
% ch=10;
%  filename='FM05_1';
load([filename '.mat']);

% PhaseFreqVector=1:250;
PhaseFreqVector=PhaseFreqVector_a:PhaseFreqVector_b;
PhaseFreq_BandWidth=4;

eeg= readmulti([filename '.lfp'], 128,ch);

data_length=length(eeg);
lfp=eeg';
srate=spkinfo.samplerate./25;

%% Do filtering and Hilbert transform on CPU

'CPU filtering'
tic
% Comodulogram=single(zeros(length(PhaseFreqVector),length(AmpFreqVeAmpFreqTransformed = zeros(length(AmpFreqVector), data_length);
PhaseFreqTransformed = zeros(length(PhaseFreqVector), data_length);

for jj=1:length(PhaseFreqVector)
    Pf1 = PhaseFreqVector(jj);
    Pf2 = Pf1 + PhaseFreq_BandWidth;
    PhaseFreq=eegfilt(lfp,srate,Pf1,Pf2); % this is just filtering 
    PhaseFreqTransformed(jj, :) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
end

%%

for ii=1:length(PhaseFreqVector)
    
ph=PhaseFreqTransformed(ii,:)';
[phase,g]=spkphase1(spk.i,spk.g,celln,spkinfo.samplerate,spkinfo.samplerate./25,ph);
n=[];
n=histc(phase,0:2*pi/18:2*pi);
n(end)=[];
nbin=length(n);
MeanAmp=n;
MI(ii)=(log(nbin)-(-sum((MeanAmp/sum(MeanAmp)).*log((MeanAmp/sum(MeanAmp))))))/log(nbin);
end
toc

end






