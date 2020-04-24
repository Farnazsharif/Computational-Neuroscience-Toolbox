%% Novel linear track for Farnaz

clearvars;clc;
dirData = 'Y:\Data\GrosmarkAD\';
% dirData = 'A:\Data\GrosmarkAD\';
sessions  = {'Achilles_10252013','Achilles_11012013','Buddy_06272013','Cicero_09012014','Cicero_09102014',...
    'Cicero_09172014','Gatsby_08282013','Gatsby_08022013'};
animal = {'Achilles','Achilles','Buddy','Cicero','Cicero','Cicero','Gatsby','Gatsby'};


%'Metrics_xtvts',   trials, xtvts,  delta=0.1; smooth_rate=1; smooth_phase=1; Binsize=2

%'Metrics_smooth_1',  trials, xtvts1, delta=0.1; smooth_rate=1; smooth_phase=1; Binsize=2

% xtvts2_3 binsize=3 has no metrics

%'CFC3' low and mid gamma ndx nad CFC for 10 trials

%'CFC2' CFC for 5 trials
 
% Cellinfo3 TXVtPh
% Cellinfo2 trash
% Cellinfo trash


for ses =1%:length(sessions)%[2 5 7],[1 3 4 6 8]%6:8
    
    ses
    dirses = [dirData animal{ses} '\' sessions{ses}];
    cd(dirses);
    load([sessions{ses} '_sessInfo.mat'])
    load([sessions{ses} '.mat'])
    filename=sessions{ses};
    
        %%
%     Fs=1250;
%     speed_thr=0.04;
%     TR=[];
%     TR=unique(behav.TXVt(:,4));
%     TR(TR==0)=[];
%     
%     if ses==2 | ses==5 | ses==7
%         trN=10;
%         
%         r0=[];TRi=[];TR_all=[];
%         TR_all=TR;Steps=1;
%         r0=[1 trN];
%         Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
%         TRi=TR_all(1:Seg*trN);
%         TRi=reshape(TRi,[trN Seg])';
%         tr=TRi;
%         
%         cd CFC3
%         load('fase_1')
%         load('amp_1')
%         cd ..
%         Com=[];
%         [Com]=CFC_linear(filename,Amp,Phase,tr,Fs);
%         
%         cd CFC3
%         save(['Com_1.mat'],'Com')
%         cd ..
%         
%         cd CFC3
%         load('fase_2')
%         load('amp_2')
%         cd .. 
%         
%         Com=[];
%         [Com]=CFC_linear(filename,Amp,Phase,tr,Fs);
%         
%         cd CFC3
%         save(['Com_2.mat'],'Com')
%         cd ..
%     else
%         
%         trN=20;
%         
%         r0=[];TRi=[];TR_all=[];
%         TR_all=TR;Steps=1;
%         r0=[1 trN];
%         Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
%         TRi=TR_all(1:Seg*trN);
%         TRi=reshape(TRi,[trN Seg])';
%         tr=TRi;
%         
%         cd CFC3
%         load('fase_1')
%         load('amp_1')
%         cd ..
%         
%         Com=[];
%         [Com]=CFC_linear(filename,Amp,Phase,tr,Fs);
%         
%         cd CFC3
%         save(['Com_1.mat'],'Com')
%         cd ..
%         
%         cd CFC3
%         load('fase_2')
%         load('amp_2')
%         cd ..
%         
%         Com=[];
%         [Com]=CFC_linear(filename,Amp,Phase,tr,Fs);
%         
%         cd CFC3
%         save(['Com_2.mat'],'Com')
%         cd ..
%     end
%     
%     
%     toc
    
    
%     Fs=1250;
%     Frq_low=[25:50];
%     Frq_mid=[65:90];
%     
%     ND=round(sessInfo.Epochs.MazeEpoch*Fs);
%     load([filename '.mat'])
%     [lfp] = bz_GetLFP('all');
%     LFPtheta_1=double(lfp.data(ND(1):ND(2),Phase.RefCh_Prob_1));
%     LFPtheta_2=double(lfp.data(ND(1):ND(2),Phase.RefCh_Prob_2));
%     
%     TR=[];
%     TR=unique(behav.TXVt(:,4));
%     TR(TR==0)=[];
%     
%     if ses==2 | ses==5 | ses==7
%         trN=5;
%         
%         r0=[];TRi=[];TR_all=[];
%         TR_all=TR;Steps=1;
%         r0=[1 trN];
%         Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
%         TRi=TR_all(1:Seg*trN);
%         TRi=reshape(TRi,[trN Seg])';
%         tr=TRi;
%         
%         lfp=[];
%         lfp=LFPtheta_1;
%         [Ndx_lowg,Ndx_midg]=gamma_finder_linear(filename,lfp,Frq_low,Frq_mid,Phase.thetaP_1,Fs,tr);
%         
%         cd CFC3
%         save(['Ndx_lowg_10.mat'],'Ndx_lowg')
%         save(['Ndx_midg_10.mat'],'Ndx_midg')
%         cd ..
%         
%         lfp=[];
%         lfp=LFPtheta_2;
%         [Ndx_lowg,Ndx_midg]=gamma_finder_linear(filename,lfp,Frq_low,Frq_mid,Phase.thetaP_2,Fs,tr);
%         
%         cd CFC3
%         save(['Ndx_lowg_20.mat'],'Ndx_lowg')
%         save(['Ndx_midg_20.mat'],'Ndx_midg')
%         cd ..
%         
%     else
%         
%         trN=10;
%         r0=[];TRi=[];TR_all=[];
%         TR_all=TR;Steps=1;
%         r0=[1 trN];
%         Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
%         TRi=TR_all(1:Seg*trN);
%         TRi=reshape(TRi,[trN Seg])';
%         tr=TRi;
%         lfp=[];
%         lfp=LFPtheta_1;
%         [Ndx_lowg,Ndx_midg]=gamma_finder_linear(filename,lfp,Frq_low,Frq_mid,Phase.thetaP_1,Fs,tr);
%         
%         cd CFC3
%         save(['Ndx_lowg_10.mat'],'Ndx_lowg')
%         save(['Ndx_midg_10.mat'],'Ndx_midg')
%         cd ..
%         
%         lfp=[];
%         lfp=LFPtheta_2;
%         [Ndx_lowg,Ndx_midg]=gamma_finder_linear(filename,lfp,Frq_low,Frq_mid,Phase.thetaP_2,Fs,tr);
%         
%         cd CFC3
%         save(['Ndx_lowg_20.mat'],'Ndx_lowg')
%         save(['Ndx_midg_20.mat'],'Ndx_midg')
%         cd ..
%         
%     end
    
    toc
end

%% Make cells and rates matrixes

%     filename=sessions{ses};

%     mkdir('Cellinfo3')
%     filename=sessions{ses};
%     Fs=1250;
%     for k=1:length(G)
%         thetaP=[];
%         CellN=G(k)
%         k
%         if Phase.probeID==3
%             if floor(CellN/100)<7
%
%                 thetaP=Phase.thetaP_1;
%             else
%                 thetaP=Phase.thetaP_2;
%             end
%
%         else
%
%             if floor(CellN/100)<9
%
%                 thetaP=Phase.thetaP_1;
%             else
%                 thetaP=Phase.thetaP_2;
%             end
%
%         end
%         TXVtPh=Cell_TXVtPh(filename,CellN,Fs,thetaP);
%         cd Cellinfo3
%         save(['TXVtPh_' num2str(k) '.mat'],'TXVtPh')
%         cd ..
%     end
%
%     clc
%     ses
%
%     display('cells are made, xtvts1 IN PROSSES')
%
%     Matrix=behav.TXVt;
%     spiket=spk.ts;
%     spikeind=spk.g;
%     phase=[Phase.thetaP_1 Phase.thetaP_1];
%     EEG_srate=1250;
%     Spike_sarte=20000;
%     [Mat_spike,Mat_phase] = makeTXYt_PhS_Antoni(filename,Matrix,spiket,spikeind,phase,EEG_srate,Spike_sarte);
%     behav.TXVt_Phas=Mat_phase;
%     behav.TXVt_rate=Mat_spike;
%     save([filename '.mat'],'-append','behav')

%     folder= 'Cellinfo3';
%     SpeedThreshold=.04;
%     xbinNumber=round(max(behav.TXVt(:,2))*100,1)/2;
%     [xtvts1,xtvts2,xtvtph]=make_rate_matrix_1D(filename,xbinNumber,SpeedThreshold,folder);
%     save([filename '.mat'],'-append','xtvts1','xtvts2','xtvtph')
%     ses
%     display('step 2 is done')

