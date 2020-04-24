
mkdir('Cellinfo1')
filename=sessions{ses};
Fs=1250;
for k=1:length(G)
    thetaP=[];
    CellN=G(k)
    k
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
    
    TXVtPh=Cell_TXVtPh(filename,CellN,Fs,thetaP);
    cd Cellinfo
    save(['TXVtPh_', num2str(k)]); % Wrong
    cd ..
end

%%
    Matrix=behav.TXVt;
    spiket=spk.ts;
    spikeind=spk.g;
    phase=[Phase.thetaP_1 Phase.thetaP_1];
    EEG_srate=1250;
    Spike_sarte=20000;
    [Mat_spike,Mat_phase] = makeTXYt_PhS_Antoni(filename,Matrix,spiket,spikeind,phase,EEG_srate,Spike_sarte);
    behav.TXVt_Phas=Mat_phase;
    behav.TXVt_rate=Mat_spike;
    save([filename '.mat'],'-append','behav')
    
    SpeedThreshold=.04;
    xbinNumber=round(max(behav.TXVt(:,2))*100,1)/2;
    [xtvts1,xtvts2,xtvtph]=make_rate_matrix_1D(filename,xbinNumber,SpeedThreshold);
    save([filename '.mat'],'-append','xtvts1','xtvts2','xtvtph')







