

    folder='CFC3';
    mkdir(folder)
    frecs1=2:0.1:10;%2:0.5:15
    frecs2=25:2:200;%20:2:150
    Fs=1250;
    load([filename '.mat'])
    [lfp] = bz_GetLFP('all');
    ND=round(sessInfo.Epochs.MazeEpoch*Fs);
    
    LFPtheta_1=double(lfp.data(ND(1):ND(2),Phase.RefCh_Prob_1));
    LFPtheta_2=double(lfp.data(ND(1):ND(2),Phase.RefCh_Prob_2));
    signals_1=double(lfp.data(ND(1):ND(2),Phase.Gam_Ch_Prob_1));
    signals_2=double(lfp.data(ND(1):ND(2),Phase.Gam_Ch_Prob_2));
    
    lfp=[];lfp=LFPtheta_1;
    LFPtheta_Filt=eegfilt(lfp',Fs,4,15);
    lfp=[];lfp=signals_1;
    signals_Filt=eegfilt(lfp',Fs,25,200);
    
    Phase=[];Amp=[];
    [Phase,Amp]=Amp_Ph_CFC(LFPtheta_Filt',signals_Filt',frecs1,frecs2,Fs);
    cd (folder)
    save(['fase_1.mat'],'Phase','-v7.3')
    save(['amp_1.mat'],'Amp','-v7.3')
    cd ..
    
    
    lfp=[];lfp=LFPtheta_2;
    LFPtheta_Filt=eegfilt(lfp',Fs,4,15);
    lfp=[];lfp=signals_2;
    signals_Filt=eegfilt(lfp',Fs,25,200);
    
    Phase=[];Amp=[];
    [Phase,Amp]=Amp_Ph_CFC(LFPtheta_Filt',signals_Filt',frecs1,frecs2,Fs);
    cd (folder)
    save(['fase_2.mat'],'Phase','-v7.3')
    save(['amp_2.mat'],'Amp','-v7.3')
    cd .. 
    toc
    
    %%
    Fs=1250;
    speed_thr=0.04;
    TR=[];
    TR=unique(behav.TXVt(:,4));
    TR(TR==0)=[];
    
    if ses==2 | ses==5 | ses==7
        trN=5;
        
        r0=[];TRi=[];TR_all=[];
        TR_all=TR;Steps=1;
        r0=[1 trN];
        Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
        TRi=TR_all(1:Seg*trN);
        TRi=reshape(TRi,[trN Seg])';
        tr=TRi;
        
        cd CFC2
        load('fase_1')
        load('amp_1')
        cd ..
        Com=[];
        [Com]=CFC_linear(filename,Amp,Phase,tr,Fs);
        
        cd CFC2
        save(['Com_1.mat'],'Com')
        cd ..
        
        cd CFC2
        load('fase_2')
        load('amp_2')
        cd .. 
        
        Com=[];
        [Com]=CFC_linear(filename,Amp,Phase,tr,Fs);
        
        cd CFC2
        save(['Com_2.mat'],'Com')
        cd ..
    else
        
        trN=10;
        
        r0=[];TRi=[];TR_all=[];
        TR_all=TR;Steps=1;
        r0=[1 trN];
        Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
        TRi=TR_all(1:Seg*trN);
        TRi=reshape(TRi,[trN Seg])';
        tr=TRi;
        
        cd CFC2
        load('fase_1')
        load('amp_1')
        cd ..
        
        Com=[];
        [Com]=CFC_linear(filename,Amp,Phase,tr,Fs);
        
        cd CFC2
        save(['Com_1.mat'],'Com')
        cd ..
        
        cd CFC2
        load('fase_2')
        load('amp_2')
        cd ..
        
        Com=[];
        [Com]=CFC_linear(filename,Amp,Phase,tr,Fs);
        
        cd CFC2
        save(['Com_2.mat'],'Com')
        cd ..
    end
    
 %%   
    toc
clearvars;clc;
dirData = 'Y:\Data\GrosmarkAD\';
% dirData = 'A:\Data\GrosmarkAD\';

sessions  = {'Achilles_10252013','Achilles_11012013','Buddy_06272013','Cicero_09012014','Cicero_09102014',...
    'Cicero_09172014','Gatsby_08282013','Gatsby_08022013'};
animal = {'Achilles','Achilles','Buddy','Cicero','Cicero','Cicero','Gatsby','Gatsby'};
j=0;
% Trial_all=[];
Com_all=[];
for ses =1:length(sessions)
    ses
    dirses = [dirData animal{ses} '\' sessions{ses}];
    cd(dirses);
    load([sessions{ses} '_sessInfo.mat'])
    load([sessions{ses} '.mat'])
    filename=sessions{ses};
  
    cd CFC3
        Com=[];
        load('Com_1.mat') 
        Com_all{1,ses}=Com; 
        Com=[];
        load('Com_2.mat')
        Com_all{2,ses}=Com;

    end
    
%% CFC statistics
a=2;
M=[];
N_l=nan(1,20);
N_m=nan(1,20);
for j=1:8
    
    size(Com_all{a, j},2)
    for i=1:size(Com_all{1, j},2)
        
%         M=[M; Com_all{a, j}{1, i}];
        M=[];
        M = Com_all{a, j}{1, i};
        lowg = mean(mean(M(13:17,12:18)));
        midg = mean(mean(M(14:17,51:58)));
        N_l(j,i)=lowg;
        N_m(j,i)=midg;
    end
%     figure('position',[100 100 1200 300])
%     imagesc(M')
%     axis xy
end

%% CFC plots
close all
a=2;
M=[];
for j=8%:8
    
    size(Com_all{a, j},2)
        for i=1:4%size(Com_all{1, j},2)
        M=[M; Com_all{a, j}{1, i}];
    end
    figure('position',[100 100 1500 300])
    imagesc(M')
    axis xy

end

%%
M=[];
M(:,:)=N_m;
Xpos=1:20;
size(M)
figure

for h=1:20
    O=[];
    O=M(:,h);
    O(find(isnan(O)==1))=[];
    a(h)= mean(O);
    s(h) = std(O);
    b1=errorbar(Xpos(h),a(h),(s(h))/sqrt(3*length(O)));
%       b1=errorbar(Xpos(h),a(h),(s(h))/sqrt(5));
    set(b1,'linewidth',2,'color','k')
    hold on
end


hold on
plot(a,'.k','markersize',30)
hold on
plot(a,'--b','linewidth',3)
% xlim([1 10])    
%% Gamma Script

Fs=1250;
Frq_low=[25:45];
Frq_mid=[60:90];

ND=round(sessInfo.Epochs.MazeEpoch*Fs);
load([filename '.mat'])
[lfp] = bz_GetLFP('all');
LFPtheta_1=double(lfp.data(ND(1):ND(2),Phase.RefCh_Prob_1));
LFPtheta_2=double(lfp.data(ND(1):ND(2),Phase.RefCh_Prob_2));

Fs=1250;
speed_thr=0.04;
TR=[];
TR=unique(behav.TXVt(:,4));
TR(TR==0)=[];

if ses==2 | ses==5 | ses==7
    trN=5;
    
    r0=[];TRi=[];TR_all=[];
    TR_all=TR;Steps=1;
    r0=[1 trN];
    Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
    TRi=TR_all(1:Seg*trN);
    TRi=reshape(TRi,[trN Seg])';
    tr=TRi;
    
    lfp=[];
    lfp=LFPtheta_1;
    [Ndx_lowg,Ndx_midg]=gamma_finder_linear(filename,lfp,Frq_low,Frq_mid,Phase.thetaP_1,Fs,tr);
    
    cd CFC2
    save(['Ndx_lowg_1.mat'],'Ndx_lowg')
    save(['Ndx_midg_1.mat'],'Ndx_midg')
    cd ..
    
    lfp=[];
    lfp=LFPtheta_1;
    [Ndx_lowg,Ndx_midg]=gamma_finder_linear(filename,lfp,Frq_low,Frq_mid,Phase.thetaP_2,Fs,tr);
    
    cd CFC2
    save(['Ndx_lowg_2.mat'],'Ndx_lowg')
    save(['Ndx_midg_2.mat'],'Ndx_midg')
    cd ..
 
else
    
    trN=10;
    
    r0=[];TRi=[];TR_all=[];
    TR_all=TR;Steps=1;
    r0=[1 trN];
    Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
    TRi=TR_all(1:Seg*trN);
    TRi=reshape(TRi,[trN Seg])';
    tr=TRi;
    lfp=[];
    lfp=LFPtheta_1;
    [Ndx_lowg,Ndx_midg]=gamma_finder_linear(filename,lfp,Frq_low,Frq_mid,Phase.thetaP_1,Fs,tr);
    
    cd CFC2
    save(['Ndx_lowg_1.mat'],'Ndx_lowg')
    save(['Ndx_midg_1.mat'],'Ndx_midg')
    cd ..
    
    lfp=[];
    lfp=LFPtheta_1;
    [Ndx_lowg,Ndx_midg]=gamma_finder_linear(filename,lfp,Frq_low,Frq_mid,Phase.thetaP_2,Fs,tr);
    
    cd CFC2
    save(['Ndx_lowg_2.mat'],'Ndx_lowg')
    save(['Ndx_midg_2.mat'],'Ndx_midg')
    cd ..
    
end

toc




    
    