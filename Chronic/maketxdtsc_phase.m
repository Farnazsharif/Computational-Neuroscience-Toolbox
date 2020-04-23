%function txdtsc=maketxdtsc(TXDTS,spkt,spkg)
%
%

function [txdtsc_phase,txdtsc_Spike]=maketxdtsc_PhS(TXDTS,gsub,filename,speed_lim)
%%
% filename='SFA4_S4_TRD1';
% TXDTS=behav.TXDTS;
load([filename '.mat'])

numchannel=64;
eeg= readmulti([filename '.lfp'],numchannel);
EEGsamplerate=1000;
PhaseFreqVector=[7];PhaseFreq_BandWidth=2;


%%
[~,b]=size(TXDTS);
txdtsc=TXDTS;
if txdtsc(1,1)<0
    txdtsc=txdtsc(find(txdtsc(:,1)>0),:);
end

[~,runT]=findrunrestT(behav.TXDTS(:,1),behav.TXDTS(:,3),speed_lim,2,1);

txdtsc_phase=txdtsc;
txdtsc_Spike=txdtsc;
%%
for ii = 1:length(gsub)
    
    [tt,gg]=selectgroup(spk.i,spk.g,gsub(ii));
    [tt]=selectt(tt,gg,runT,spkinfo.samplerate);
    tt=tt./spkinfo.samplerate;
    tt=tt(find(tt>=txdtsc(1,1)&tt<=txdtsc(end,1)));
    
    % ********************************** % Select local channel or refrence
    ch=Cellinfo.Cellchn_T(ii,2);
    lfp=eeg(:,ch)';
    [phase]=cellph_Overal(filename,PhaseFreqVector,PhaseFreq_BandWidth,EEGsamplerate,lfp,tt);
    
    
    % ********************************** % Bin data
    [bincounts,ind] = histc(tt,txdtsc(:,1));
    
    % ********************************** % Spike_count
    txdtsc_Spike(:,b+ii)=bincounts;
    
    % ********************************** % Phase_count
    ndx_1_spk=[];ndx_multi_spk=[];Ph=[];
    
    ndx_1_spk=find(bincounts>0 & bincounts<2);
    Ph=phase(ismember(ind,ndx_1_spk));
    
    ndx_multi_spk=find(bincounts>1);
    for j=1:length(ndx_multi_spk)
        [~,ndx_multi_spk(j,2)]= meanphase(phase(ismember(ind,ndx_multi_spk(j))) );
    end
    
    txdtsc_phase(ndx_1_spk,b+ii)=Ph;
    txdtsc_phase(ndx_multi_spk(:,1),b+ii)=ndx_multi_spk(:,2);
 
end

%%





