
function [phase]=spkphase_period(T,G,gsubset,smplrate,eegsmplrate,ph,tt)

%     [spikeT,spikeG]=selectgroup(T,G,gsubset);
   
%     [tt,g]=selectt(spikeT,spikeG,period,smplrate);
    
    ph = ph + cumsum([0;diff(ph)<0])*2*pi;

	ph = interp1((1:length(ph))/eegsmplrate,ph,tt/smplrate);
    
    phase = rem(ph,2*pi);