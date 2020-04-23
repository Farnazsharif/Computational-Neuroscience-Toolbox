

function CellMI(filename,Cellchn,Cell_Start,numchannel)
%%
% ProbN=2;
% cellN=1;
% filename='FM05_1';
load([filename '.mat']);
%%
% [Cellchn]=Cellsite(ProbN);
PhaseFreqVector_b=250;
PhaseFreqVector_a=1;
eeg= readmulti([filename '.lfp'],numchannel);
%   [T,Ndx]=binT(filename,EEGsamplerate,numchannel);
for i=Cell_Start:length(Cellchn)
    tic
    i
    celln=Cellchn(i,1);
    ch=Cellchn(i,2);
%   [MI]=cellph(celln,ch,filename);
if  isnan(sum(sum (spkinfo.waveform.average(:,:,i))))==1
%     down= [13 20 39 53 78  83 91 115 127 153 155 163 164 165 174 177  235 237 246 ];
    MI=nan(PhaseFreqVector_a,PhaseFreqVector_b);
    MI_run=nan(PhaseFreqVector_a,PhaseFreqVector_b);
    MI_rest=nan(PhaseFreqVector_a,PhaseFreqVector_b);

else
    [MI,MI_rest, MI_run]=cellph_Overal(celln,filename,PhaseFreqVector_a,PhaseFreqVector_b,spkinfo.samplerate/30,eeg(:,ch));
end

    cd MI_Probs
    save  (['MI_', num2str(i)],'MI')    
    save  (['MI_rest_', num2str(i)],'MI_rest')
    save  (['MI_run_', num2str(i)],'MI_run')
    cd ..
    toc
end

%%
end