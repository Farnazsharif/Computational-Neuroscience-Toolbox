 %#######################################################################################################################################   
%Gamma counting script

    Fs=1250;
    Frq_low=[25:50];
    Frq_mid=[65:90];
    
    ND=round(sessInfo.Epochs.MazeEpoch*Fs);
    load([filename '.mat'])
    [lfp] = bz_GetLFP('all');
    LFPtheta_1=double(lfp.data(ND(1):ND(2),Phase.RefCh_Prob_1));
    LFPtheta_2=double(lfp.data(ND(1):ND(2),Phase.RefCh_Prob_2));
    
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
        
        cd CFC3
        save(['Ndx_lowg_1.mat'],'Ndx_lowg')
        save(['Ndx_midg_1.mat'],'Ndx_midg')
        cd ..
        
        lfp=[];
        lfp=LFPtheta_2;
        [Ndx_lowg,Ndx_midg]=gamma_finder_linear(filename,lfp,Frq_low,Frq_mid,Phase.thetaP_2,Fs,tr);
        
        cd CFC3
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
        
        cd CFC3
        save(['Ndx_lowg_1.mat'],'Ndx_lowg')
        save(['Ndx_midg_1.mat'],'Ndx_midg')
        cd ..
        
        lfp=[];
        lfp=LFPtheta_2;
        [Ndx_lowg,Ndx_midg]=gamma_finder_linear(filename,lfp,Frq_low,Frq_mid,Phase.thetaP_2,Fs,tr);
        
        cd CFC3
        save(['Ndx_lowg_2.mat'],'Ndx_lowg')
        save(['Ndx_midg_2.mat'],'Ndx_midg')
        cd ..
        
    end
    
    toc

%% #######################################################################################################################################    
% Counting gamma

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
        Ndx_lowg=[];
        load('Ndx_lowg_1.mat') 
        Lowg_all{1,ses}=Ndx_lowg; 
        Ndx_midg=[];
        load('Ndx_midg_1.mat')
        Midg_all{1,ses}=Ndx_midg;
        
        Ndx_lowg=[];
        load('Ndx_lowg_2.mat') 
        Lowg_all{2,ses}=Ndx_lowg; 
        Ndx_midg=[];
        load('Ndx_midg_2.mat')
        Midg_all{2,ses}=Ndx_midg;
end
  
%% #######################################################################################################################################   


clearvars;clc;
dirData = 'Y:\Data\GrosmarkAD\';
% dirData = 'A:\Data\GrosmarkAD\';

sessions  = {'Achilles_10252013','Achilles_11012013','Buddy_06272013','Cicero_09012014','Cicero_09102014',...
             'Cicero_09172014','Gatsby_08282013','Gatsby_08022013'};
animal = {'Achilles','Achilles','Buddy','Cicero','Cicero','Cicero','Gatsby','Gatsby'};
clearvars;clc;
dirData = 'Y:\Data\GrosmarkAD\';
% dirData = 'A:\Data\GrosmarkAD\';

sessions  = {'Achilles_10252013','Achilles_11012013','Buddy_06272013','Cicero_09012014','Cicero_09102014',...
             'Cicero_09172014','Gatsby_08282013','Gatsby_08022013'};
animal = {'Achilles','Achilles','Buddy','Cicero','Cicero','Cicero','Gatsby','Gatsby'};

Trials_Seg=[];
for ses =1:length(sessions)
    ses
    dirses = [dirData animal{ses} '\' sessions{ses}];
    cd(dirses);
    load([sessions{ses} '_sessInfo.mat'])
    load([sessions{ses} '.mat'])
    filename=sessions{ses};

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
        
            else
        
        trN=10;
        r0=[];TRi=[];TR_all=[];
        TR_all=TR;Steps=1;
        r0=[1 trN];
        Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
        TRi=TR_all(1:Seg*trN);
        TRi=reshape(TRi,[trN Seg])';
        tr=TRi; 
    end
    
    TRN=[];
    
    for i=1:size(tr,1)
       TRN{1,i} =tr(i,:);
    end
    
    TRN{1,size(tr,1)+1}=behav.TXVt;
    Trials_Seg{1,ses}=TRN;
    
end   

%%
save(['Trials_Seg.mat'],'Trials_Seg')
%% CFC statistics
a=1;r=150;
N_ratio=nan(8,20);
N_Low=nan(8,20);
N_Mid=nan(8,20);

for j=1:8
    Seg=size(Lowg_all{a, j},2)
    Ndx_g_lg=[];Ndx_g_hg=[];
    
    Ndx_g_lg=Lowg_all{a, j}{1, Seg};
    Ndx_g_hg=Midg_all{a, j}{1, Seg};
    [~,C_L]=setdiff(Ndx_g_lg(:,1),Ndx_g_hg(:,1));%delet mid gamma
    [~,C_h]=setdiff(Ndx_g_hg(:,1),Ndx_g_lg(:,1));%delet low gamma
    Ndx_g_lg=Ndx_g_lg(C_L,:);
    Ndx_g_hg=Ndx_g_hg(C_h,:);
    
    d=[];M=[];
    M=Ndx_g_hg;
    for i=1:size(Ndx_g_lg,1)
        idx=rangesearch(M(:,1),Ndx_g_lg(i,1),r);
        d=[d idx{1}];
    end
    M(d,:)=[];
    
            N_Mid_all(j)=size(M,1);
% for Power ###############################################################
%     if isempty(M)==0
%         
%         N_Mid_all(j)=max(M(:,2));
%     else
%         N_Mid_all(j)=nan;
%     end
% #########################################################################    
    
    M=[];d=[];
    M=Ndx_g_lg;
    for i=1:size(Ndx_g_hg,1)
        idx=rangesearch(M(:,1),Ndx_g_hg(i,1),r);
        d=[d idx{1}];
    end
    M(d,:)=[];
            N_Low_all(j)=size(M,1);
            
% for Power ############################################################### 
%     if isempty(M)==0
%         
%         N_Low_all(j)=max(M(:,2));
%     else
%         N_Low_all(j)=nan;
%     end
% #########################################################################    
    
    
    
    
    for k=1:(Seg-1)
        Ndx_g_lg=[];Ndx_g_hg=[];
        
        Ndx_g_lg=Lowg_all{a, j}{1, k};
        Ndx_g_hg=Midg_all{a, j}{1, k};
        [~,C_L]=setdiff(Ndx_g_lg(:,1),Ndx_g_hg(:,1));%delet mid gamma
        [~,C_h]=setdiff(Ndx_g_hg(:,1),Ndx_g_lg(:,1));%delet low gamma
        Ndx_g_lg=Ndx_g_lg(C_L,:);
        Ndx_g_hg=Ndx_g_hg(C_h,:);
        
        d=[];M=[];
        M=Ndx_g_hg;
        for i=1:size(Ndx_g_lg,1)
            idx=rangesearch(M(:,1),Ndx_g_lg(i,1),r);
            d=[d idx{1}];
        end
        M(d,:)=[];
        
                N_Mid(j,k)=size(M,1);
        %         N_Mid(j,k)=size(M,1)./N_Mid_all(j);
% for Power ############################################################### 

%         N_Mid(j,k)=mean(M(:,2)./N_Mid_all(j));
    
        
        M=[];d=[];
        M=Ndx_g_lg;
        for i=1:size(Ndx_g_hg,1)
            idx=rangesearch(M(:,1),Ndx_g_hg(i,1),r);
            d=[d idx{1}];
        end
        M(d,:)=[];
        
                N_Low(j,k)=size(M,1);
        %         N_Low(j,k)=size(M,1)/N_Low_all(j);
% for Power ############################################################### 
%         N_Low(j,k)=mean(M(:,2)./N_Low_all(j));
%########################################################################## 
        
        N_ratio(j,k)=N_Mid(j,k)./N_Low(j,k);
        
% for Power ############################################################### 
%        N_ratio(j,k)=(N_Mid(j,k)/N_Mid_all(j))./(N_Low(j,k)/N_Low_all(j));
%##########################################################################     

    end
end


%%
% close all
hold on
% figure('position',[200 200 250 250])
M=[];

% M(:,:)=N_ratio;
M(:,:)=N_Low;
% M(:,:)=N_Mid;
% M(:,:)=N_Low;M(:,9)=M(:,9)+.1;
% M(2,:)=[];
% M(:,5)=[];

Xpos=1:19;
size(M)

for h=1:19
    O=[];
    O=M(:,h);
    O(find(isnan(O)==1))=[];
    a(h)= mean(O);
    s(h) = std(O);

      b1=errorbar(Xpos(h),a(h),(s(h))/sqrt(40*length(O)));
%     if h==10
%     b1=errorbar(Xpos(h),a(h),(s(h))/sqrt(200*length(O)));   
%     else
%         b1=errorbar(Xpos(h),a(h),(s(h))/sqrt(40*length(O)));
%     end

    set(b1,'linewidth',2,'color','k')
    hold on
end


hold on
plot(a,'.k','markersize',20)
hold on
plot(a,':b','linewidth',2)
xlim([0 10])

%%
clc
a=1;
Ndx_g_lg=[];Ndx_g_hg=[];
Ndx_g_lg=Lowg_all{a, 7}{1, 17};
Ndx_g_hg=Midg_all{a, 7}{1, 17};

%% separat gammas

CA3_Phase=[degtorad(45) degtorad(160)];
EC3_Phase=[degtorad(200) degtorad(360)];


[~,C_L]=setdiff(Ndx_g_lg(:,1),Ndx_g_hg(:,1));%delet mid gamma
[~,C_h]=setdiff(Ndx_g_hg(:,1),Ndx_g_lg(:,1));%delet low gamma
[~,C_bL]=intersect(Ndx_g_hg(:,1),Ndx_g_lg(:,1)); % both gamma 
[~,C_bh]=intersect(Ndx_g_lg(:,1),Ndx_g_hg(:,1)); % both gamma  

size(C_L)
size(C_h)
size(C_bh)
size(C_bL)
Ndx_g_lg=Ndx_g_lg(C_L,:);
Ndx_g_hg=Ndx_g_hg(C_h,:);

%% separat mid gammas
l=2;
clc
d=[];
r=100;
M=[];
M=Ndx_g_hg;
for i=1:size(Ndx_g_lg,1)
        idx=rangesearch(M(:,1),Ndx_g_lg(i,1),r);
        d=[d idx{1}];
end

size(M)
M(d,:)=[];
size(M)

ND_TR=nan(100,2);
TR=unique(M(:,6));
for jj=1:length(TR)
    nd=[];
    nd=find(M(:,6)==TR(jj));
    ND_TR(TR(jj),2)=length(nd);
    ND_TR(TR(jj),1)=TR(jj);
    ND_TR(TR(jj),1)=TR(jj);
    ND_TR(TR(jj),3)=mean(M(nd,2));
end
ND_TR_mid=ND_TR;

% seperat phase
ph_g=EC3_Phase;
ndx_ph=find(M(:,3)>= ph_g(1,1) & M(:,3)<= ph_g(1,2));
M=M(ndx_ph,:);
size(M)
ND_TR=nan(100,2);
TR=unique(M(:,6));
for jj=1:length(TR)
    nd=[];
    nd=find(M(:,6)==TR(jj));
    ND_TR(TR(jj),2)=length(nd);
    ND_TR(TR(jj),1)=TR(jj);
    ND_TR(TR(jj),3)=mean(M(nd,2));
end
ND_mid_phase=ND_TR;


% separat low gammas
l=2;
clc
d=[];

M=[];
M=Ndx_g_lg;
for i=1:size(Ndx_g_hg,1)
        idx=rangesearch(M(:,1),Ndx_g_hg(i,1),r);
        d=[d idx{1}];
end

size(M)
M(d,:)=[];
size(M)

ND_TR=nan(100,2);
TR=unique(M(:,6));
for jj=1:length(TR)
    nd=[];
    nd=find(M(:,6)==TR(jj));
    ND_TR(TR(jj),2)=length(nd);
    ND_TR(TR(jj),1)=TR(jj);
    ND_TR(TR(jj),3)=mean(M(nd,2));
end
ND_TR_low=ND_TR;

% seperat phase
ph_g=EC3_Phase;
ndx_ph=find(M(:,3)>= ph_g(1,1) & M(:,3)<= ph_g(1,2));
M=M(ndx_ph,:);
size(M)
ND_TR=nan(100,2);
TR=unique(M(:,6));
for jj=1:length(TR)
    nd=[];
    nd=find(M(:,6)==TR(jj));
    ND_TR(TR(jj),2)=length(nd);
    ND_TR(TR(jj),1)=TR(jj);
    ND_TR(TR(jj),3)=mean(M(nd,2));
end
ND_low_phase=ND_TR;
%%

%%

close all
figure('position',[100 0 200 800])
event1=ND_TR_mid;
event1_N=event1(:,2)./length(event1);

event2=ND_mid_phase;
event2_N=event2(:,2)./length(event2);
xlim([1 50])

subplot(8,1,1)
plot(event1(:,1),event1(:,2),'.-k','markersize',10,'linewidth',2)
hold on
plot(event2(:,1),event2(:,2),'.-r','markersize',10,'linewidth',2)
title('Midg N')
xlim([1 50])

subplot(8,1,2)
plot(event1(:,1),event1_N,'.-k','markersize',10,'linewidth',2)
hold on
plot(event2(:,1),event2_N,'.-r','markersize',10,'linewidth',2)
title('Normal Midg N')
xlim([1 50])

subplot(8,1,3)
plot(event1(:,1),event1(:,3),'.-k','markersize',10,'linewidth',2)
hold on
plot(event2(:,1),event2(:,3),'.-r','markersize',10,'linewidth',2)
title('Midg Power')
xlim([1 50])

event1=ND_TR_low;
event1_N=event1(:,2)./length(event1);

event2=ND_low_phase;
event2_N=event2(:,2)./length(event2);

subplot(8,1,4)
plot(event1(:,1),event1(:,2),'.-k','markersize',10,'linewidth',2)
hold on
plot(event2(:,1),event2(:,2),'.-r','markersize',10,'linewidth',2)
title('Lowg N')
xlim([1 50])

subplot(8,1,5)
plot(event1(:,1),event1_N,'.-k','markersize',10,'linewidth',2)
hold on
plot(event2(:,1),event2_N,'.-r','markersize',10,'linewidth',2)
title('Normal Lowg N')
xlim([1 50])

subplot(8,1,6)
plot(event1(:,1),event1(:,3),'.-k','markersize',10,'linewidth',2)
hold on
plot(event2(:,1),event2(:,3),'.-r','markersize',10,'linewidth',2)
title('Lowgamma Power')
xlim([1 50])

subplot(8,1,7)

M_mid=ND_TR_mid;
M_low=ND_TR_low;
M_mid(isnan(M_mid(:,1))==1,:)=0;
M_low(isnan(M_low(:,1))==1,:)=0;
plot(M_mid(:,1),M_mid(:,2)./M_low(:,2),'.k','markersize',10,'linewidth',2)
title('Lowg/Midg')
xlim([1 50])
subplot(8,1,8)

M_mid=ND_mid_phase;
M_low=ND_low_phase;
M_mid(isnan(M_mid(:,1))==1,:)=0;
M_low(isnan(M_low(:,1))==1,:)=0;
plot(M_mid(:,1),M_mid(:,2)./M_low(:,2),'.r','markersize',10,'linewidth',2)
title('Lowg/Midg-ph')
xlabel(['r' num2str(r)])
xlim([1 50])

%% Gamma statistics
a=1;r=150;
N_ratio=nan(8,20);
N_Low=nan(8,20);
N_Mid=nan(8,20);

for j=1:8
    Seg=size(Lowg_all{a, j},2)
    Ndx_g_lg=[];Ndx_g_hg=[];
    
    Ndx_g_lg=Lowg_all{a, j}{1, Seg};
    Ndx_g_hg=Midg_all{a, j}{1, Seg};
    [~,C_L]=setdiff(Ndx_g_lg(:,1),Ndx_g_hg(:,1));%delet mid gamma
    [~,C_h]=setdiff(Ndx_g_hg(:,1),Ndx_g_lg(:,1));%delet low gamma
    Ndx_g_lg=Ndx_g_lg(C_L,:);
    Ndx_g_hg=Ndx_g_hg(C_h,:);
    
    d=[];M=[];
    M=Ndx_g_hg;
    for i=1:size(Ndx_g_lg,1)
        idx=rangesearch(M(:,1),Ndx_g_lg(i,1),r);
        d=[d idx{1}];
    end
    M(d,:)=[];
    
    %         N_Mid_all(j)=size(M,1);
    if isempty(M)==0
        
        N_Mid_all(j)=max(M(:,2));
    else
        N_Mid_all(j)=nan;
    end
    
     
    M=[];d=[];
    M=Ndx_g_lg;
    for i=1:size(Ndx_g_hg,1)
        idx=rangesearch(M(:,1),Ndx_g_hg(i,1),r);
        d=[d idx{1}];
    end
    M(d,:)=[];
    %         N_Low_all(j)=size(M,1);
    
    if isempty(M)==0
        
        N_Low_all(j)=max(M(:,2));
    else
        N_Low_all(j)=nan;
    end
    
    
    
    
    for k=1:(Seg-1)
        Ndx_g_lg=[];Ndx_g_hg=[];
        
        Ndx_g_lg=Lowg_all{a, j}{1, k};
        Ndx_g_hg=Midg_all{a, j}{1, k};
        [~,C_L]=setdiff(Ndx_g_lg(:,1),Ndx_g_hg(:,1));%delet mid gamma
        [~,C_h]=setdiff(Ndx_g_hg(:,1),Ndx_g_lg(:,1));%delet low gamma
        Ndx_g_lg=Ndx_g_lg(C_L,:);
        Ndx_g_hg=Ndx_g_hg(C_h,:);
        
        d=[];M=[];
        M=Ndx_g_hg;
        for i=1:size(Ndx_g_lg,1)
            idx=rangesearch(M(:,1),Ndx_g_lg(i,1),r);
            d=[d idx{1}];
        end
        M(d,:)=[];
        
        %         N_Mid(j,k)=size(M,1);
        %         N_Mid(j,k)=size(M,1)./N_Mid_all(j);
        N_Mid(j,k)=mean(M(:,2)./N_Mid_all(j));
        
        
        M=[];d=[];
        M=Ndx_g_lg;
        for i=1:size(Ndx_g_hg,1)
            idx=rangesearch(M(:,1),Ndx_g_hg(i,1),r);
            d=[d idx{1}];
        end
        M(d,:)=[];
        
        %         N_Low(j,k)=size(M,1);
        %         N_Low(j,k)=size(M,1)/N_Low_all(j);
        N_Low(j,k)=mean(M(:,2)./N_Low_all(j));
        
        N_ratio(j,k)=N_Mid(j,k)./N_Low(j,k);
        %         N_ratio(j,k)=(N_Mid(j,k)/N_Mid_all(j))./(N_Low(j,k)/N_Low_all(j));
        
        
    end
end
%%


save(['Midg_all.mat'],'Midg_all')
save(['Lowg_all.mat'],'Lowg_all')






























