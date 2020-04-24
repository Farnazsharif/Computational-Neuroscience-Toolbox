function [Mat_spike,Mat_phase] = makeTXYt_PhS_Antoni(filename,Matrix,spiket,spikeind,phase,EEG_srate,Spike_sarte)
%%
% spiket=spk.ts;
% spikeind=spk.g;
% Matrix=behav.TXVt;
% Spike_sarte=20000;
% EEG_srate=1250;
% phase=thetaP;
%%
SR=Spike_sarte;
Fs=EEG_srate;
[~,b]=size(Matrix);
CellGroup=unique(spikeind);
%%
Mat_phase=[];Mat_spike=[];
Mat_phase=nan(size(Matrix,1),b+length(CellGroup));
Mat_spike=zeros(size(Matrix,1),b+length(CellGroup));
Mat_spike(:,1:4)=Matrix;
Mat_phase(:,1:4)=Matrix;
%%
for ii =1:length(CellGroup)
    SPh=[];ind=[];bincounts=[]; 
    SPh(:,1)=spiket(spikeind==CellGroup(ii));
    if size(phase,2)>1
        
        CellN=CellGroup(ii);
        load([filename '.mat'])
        
        if Phase.probeID==3
            if floor(CellN/100)<7
                
                thetaP=Phase.thetaP_1;
            else
                thetaP=Phase.thetaP_2;
            end
            
        else
            
            if floor(CellN/100)<9
                
                thetaP=Phase.thetaP_1;
            else
                thetaP=Phase.thetaP_2;
            end
            
        end
        
    else
        thetaP=phase;
    end
    
    
    
%     SPh(:,2)=mod(thetaP(round(SPh(:,1)*Fs/SR)),2*pi);
    SPh(:,2)=mod(thetaP(round(SPh(:,1)*Fs)),2*pi);
%     SPh(:,1)=SPh(:,1)./SR;
    SPh(:,1)=SPh(:,1);
    SPh=SPh(find(SPh(:,1)>=Matrix(1,1) & SPh(:,1)<=Matrix(end,1)),:);
    % ********************************** % Bin data
    [bincounts,ind] = histc(SPh(:,1),Matrix(:,1));
    % ********************************** % Spike_count
    Mat_spike(:,b+ii)=[bincounts];  
    AS=[];AS=unique(ind);
    for kk=1:length(AS)
        [~,AS(kk,2)]=meanphase(SPh(find(ind==AS(kk)),2));
    end
    if isempty(AS)==1
        Mat_phase(:,b+ii)=nan(size(bincounts,1),1);
        
    else
    Mat_phase(AS(:,1),b+ii)=AS(:,2);
    end
end

end

