
clearvars;clc;
% dirData = 'Y:\Data\GrosmarkAD\';
dirData = 'A:\Data\GrosmarkAD\';
sessions  = {'Achilles_10252013','Achilles_11012013','Buddy_06272013','Cicero_09012014','Cicero_09102014',...
    'Cicero_09172014','Gatsby_08282013','Gatsby_08022013'};
animal = {'Achilles','Achilles','Buddy','Cicero','Cicero','Cicero','Gatsby','Gatsby'};


for ses =2:length(sessions)%[2 5 7],[1 3 4 6 8]%6:8
    
    ses
    dirses = [dirData animal{ses} '\' sessions{ses}];
    cd(dirses);
    load([sessions{ses} '_sessInfo.mat'])
    load([sessions{ses} '.mat'])
    filename=sessions{ses};
    

    folder='CFC3';
    mkdir(folder)
    Frq=[4:15];
    
    Fs=1250;
    load([filename '.mat'])
    [lfp] = bz_GetLFP('all');
    ND=round(sessInfo.Epochs.MazeEpoch*Fs);
    
    LFPtheta_1=double(lfp.data(ND(1):ND(2),Phase.RefCh_Prob_1));
    LFPtheta_2=double(lfp.data(ND(1):ND(2),Phase.RefCh_Prob_2));
    
    S=[];Theta_P1=[];
    scale=frq2scal(Frq,Fs);
    S=cwt(LFPtheta_1,scale,'morl');
    Theta_P1= (envelop(S.*S))';
    
    S=[];Theta_P2=[];
    S=cwt(LFPtheta_2,scale,'morl');
    Theta_P2= (envelop(S.*S))';
    
    cd (folder)
    
    save(['Theta_P1.mat'],'Theta_P1','-v7.3')
    save(['Theta_P2.mat'],'Theta_P2','-v7.3')
    
    cd .. 

   toc
end
%%


clearvars;clc;
dirData = 'Y:\Data\GrosmarkAD\';
% dirData = 'A:\Data\GrosmarkAD\';
sessions  = {'Achilles_10252013','Achilles_11012013','Buddy_06272013','Cicero_09012014','Cicero_09102014',...
    'Cicero_09172014','Gatsby_08282013','Gatsby_08022013'};
animal = {'Achilles','Achilles','Buddy','Cicero','Cicero','Cicero','Gatsby','Gatsby'};


for ses =1:length(sessions)%[2 5 7],[1 3 4 6 8]%6:8
    
    ses
    dirses = [dirData animal{ses} '\' sessions{ses}];
    cd(dirses);
    load([sessions{ses} '_sessInfo.mat'])
    load([sessions{ses} '.mat'])
    filename=sessions{ses};
    
Fs=1250;
speed_thr=0.04;
TR=[];
TR=unique(behav.TXVt(:,4));
TR(TR==0)=[];

NDX_sync=[];
ND=round(sessInfo.Epochs.MazeEpoch*Fs);
NDX_sync=ND(1)-1;

for x=1:size(TR,1)
    nd=find(behav.TXVt(:,4)==TR(x));
    Trial_ND(x,1)=round(behav.TXVt(nd(1),1)*Fs)-NDX_sync;
    Trial_ND(x,2)=round(behav.TXVt(nd(end),1)*Fs)-NDX_sync;
end



if ses==2 | ses==5 | ses==7
    trN=5;
    
    r0=[];TRi=[];TR_all=[];
    TR_all=TR;Steps=1;
    r0=[1 trN];
    Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
    TRi=TR_all(1:Seg*trN);
    TRi=reshape(TRi,[trN Seg])';
    tr=TRi;
    
    cd CFC3
    load('Theta_P1')
    load('Theta_P2')
    cd ..
    
    Pow=Theta_P1;
    P_t=[];
    for x=1:size(tr,1)
        
        for y=1:size(tr,2)
            P=[];
            P=(Pow(Trial_ND(tr(x,y),1):Trial_ND(tr(x,y),2),3:5));
            P_t(x,y)=mean(mean(P));
        end
        
    end
    Theta_Pow{ses,1}=P_t;
    
    Pow=Theta_P2;
    P_t=[];
    for x=1:size(tr,1)
        
        for y=1:size(tr,2)
            P=[];
            P=(Pow(Trial_ND(tr(x,y),1):Trial_ND(tr(x,y),2),3:5));
            P_t(x,y)=mean(mean(P));
        end
        
    end
    Theta_Pow{ses,2}=P_t;
    
else
    
    
    trN=10;
    
    r0=[];TRi=[];TR_all=[];
    TR_all=TR;Steps=1;
    r0=[1 trN];
    Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
    TRi=TR_all(1:Seg*trN);
    TRi=reshape(TRi,[trN Seg])';
    tr=TRi;
    
    cd CFC3
    load('Theta_P1')
    load('Theta_P2')
    cd ..
    
    Pow=Theta_P1;
    P_t=[];
    for x=1:size(tr,1)
        
        for y=1:size(tr,2)
            P=[];
            P=(Pow(Trial_ND(tr(x,y),1):Trial_ND(tr(x,y),2),3:5));
            P_t(x,y)=mean(mean(P));
        end
        
    end
    Theta_Pow{ses,1}=P_t;
    
    Pow=Theta_P2;
    P_t=[];
    for x=1:size(tr,1)
        
        for y=1:size(tr,2)
            P=[];
            P=(Pow(Trial_ND(tr(x,y),1):Trial_ND(tr(x,y),2),3:5));
            P_t(x,y)=mean(mean(P));
        end
        
    end
    Theta_Pow{ses,2}=P_t;
    
    
end

end

%%


save(['Theta_Pow'],'Theta_Pow')












