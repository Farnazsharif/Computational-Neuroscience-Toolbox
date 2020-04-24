    function [Ndx_g_low,Ndx_g_mid]=gamma_finder_linear(filename,lfp,Frq_low,Frq_mid,fase,Fs,tr)

    load([filename '_sessInfo.mat'])
    load([filename '.mat'])
    Speed_thr=0.035;
    Thr=2;
    AmpFreq_BandWidth=4;
    ND=round(sessInfo.Epochs.MazeEpoch*Fs);
    
%% ########################################################################  

    Frq=Frq_low;

    Ndx_g=[];
    Ndx_g=gamma_win(lfp,Fs,Frq,AmpFreq_BandWidth,Thr);
    Ndx_g(:,1)=Ndx_g(:,1)+ND(1)-1;
    Ndx_g(:,3)=fase(Ndx_g(:,1));
    
    nd=[];
    for i=1:size(Ndx_g,1)
        [~,nd(1,i)]=min(abs((Ndx_g(i,1)/Fs)-behav.TXVt(:,1)));
    end
    Ndx_g(:,4:6)=behav.TXVt(nd,2:4);

    Ndx_g((Ndx_g(:,6)==0),:)=[];
    Ndx_g((Ndx_g(:,5)<Speed_thr),:)=[];
    
    g=[];
    for i=1:size(tr,1)
        nd=[];
        for j=1:size(tr,2)
            nd=[nd; find(Ndx_g(:,6)==tr(i,j))];
        end
        g{1,i}=Ndx_g(nd,:);
    end
    g{1,i+1}=Ndx_g;
    Ndx_g_low=g;

    
%% ########################################################################  
    Frq=Frq_mid;
    
    Ndx_g=[];
    Ndx_g=gamma_win(lfp,Fs,Frq,AmpFreq_BandWidth,Thr);
    Ndx_g(:,1)=Ndx_g(:,1)+ND(1)-1;
    Ndx_g(:,3)=fase(Ndx_g(:,1));
    nd=[];
    for i=1:size(Ndx_g,1)
        [~,nd(1,i)]=min(abs((Ndx_g(i,1)/Fs)-behav.TXVt(:,1)));
    end
    
    Ndx_g(:,4:6)=behav.TXVt(nd,2:4);
    Ndx_g((Ndx_g(:,6)==0),:)=[];
    Ndx_g((Ndx_g(:,5)<Speed_thr),:)=[];  
    
    g=[];
    for i=1:size(tr,1)
        nd=[];
        for j=1:size(tr,2)
            nd=[nd; find(Ndx_g(:,6)==tr(i,j))];
        end
        g{1,i}=Ndx_g(nd,:);
    end
    g{1,i+1}=Ndx_g;

    Ndx_g_mid=g;
    
    
    
    
    