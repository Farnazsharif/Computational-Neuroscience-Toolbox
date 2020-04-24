
filename='Achilles_10252013';

Spk = LoadSpk(filename, 10, 32, 1);

'step 1 done'
%%
NumberProbs=2;
probeID=3;
Winrange=100;
ShankN=6;%16;
eeg_smplrate=1250;
smplrate=20000;
Elnb=6;nbchperEl=10;
%%
spk=[];
spkinfo.samplerate=20000;
spk.g=sessInfo.Spikes.SpikeIDs;
spk.ts=sessInfo.Spikes.SpikeTimes;
spk.i=round(sessInfo.Spikes.SpikeTimes*20000);


