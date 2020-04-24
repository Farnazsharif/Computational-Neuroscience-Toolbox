function [Mat_spike,Mat_phase] = makeTXYtc_PhS_Antoni(Matrix,spiket,spikeind,phase,EEG_srate,Spike_sarte)
%%
% spiket=spike_time;
% spikeind=spike_ind;
% Matrix=TXYt;
%%
SR=Spike_sarte;
Fs=EEG_srate;
thetaP=phase;
Mat_spike=Matrix;
Mat_phase=Matrix;


%%
[~,b]=size(Matrix);
Cellgroup=unique(spikeind);

for ii =1:length(Cellgroup)
    
    SPh=[];ind=[];bincounts=[]; 
    SPh(:,1)=spiket(spikeind==Cellgroup(ii));
    
    if max(SPh(:,1))./SR > size(thetaP,1)./Fs  
        bad=find(SPh(:,1)./SR > size(thetaP,1)./Fs);
        SPh(bad)=[];
    end
    
    SPh(:,2)=mod(thetaP(round(SPh(:,1)*Fs/SR)),2*pi);
    SPh(:,1)=SPh(:,1)./SR;
    SPh=SPh(find(SPh(:,1)>=Matrix(1,1) & SPh(:,1)<=Matrix(end,1)),:);
    % ********************************** % Bin data
    [bincounts,ind] = histc(SPh(:,1),Matrix(:,1));
    % ********************************** % Spike_count
    Mat_spike(:,b+ii)=[ bincounts];  
    AS=[];AS=unique(ind);
    for kk=1:length(AS)
        [~,AS(kk,2)]=meanphase(SPh(find(ind==AS(kk)),2));
    end
    Mat_phase(AS(:,1),b+ii)=AS(:,2);
end

end

