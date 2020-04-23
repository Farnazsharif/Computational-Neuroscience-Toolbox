
Rf_ch_shank=[];
Rf_ch_shank=Ripple.RefChannelAvg(:,1);
channel_order=[];
channel_order{1}=[12 11 10 9 8 7 6 5 4 3 2 1];
channel_order{2}=[16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1];

for CellN=1:spikes.numcells
    CellN
    CellPosition=[];
    Cellchannel=[spikes.Cellinfo(CellN).shankID spikes.Cellinfo(CellN).peakChOrganized];
    
    if Cellchannel(1,1)==3
        nd=find(spikes.channelXML{Cellchannel(1,1),2}==spikes.Cellinfo(CellN).peakChOrganized);
        CellPosition{2}=[spikes.Cellinfo(CellN).shankID channel_order{2}(nd)];
        CellPosition{1}=[nan nan];
    elseif Cellchannel(1,1)>3
        nd=find(spikes.channelXML{Cellchannel(1,1),2}==spikes.Cellinfo(CellN).peakChOrganized);
        CellPosition{1}=[spikes.Cellinfo(CellN).shankID channel_order{1}(nd)];
        CellPosition{2}=[nan nan];
    else
        nd=find(spikes.channelXML{Cellchannel(1,1),2}==spikes.Cellinfo(CellN).peakChOrganized);
        CellPosition{1}=[spikes.Cellinfo(CellN).shankID channel_order{1}(nd)];
        CellPosition{2}=[nan nan];
    end
    
    probeID=4;
    [CellShankN,CellChannelN,RefCellChannelN,RefCellShankN]=Clayout_AntonioProb(CellPosition,Rf_ch_shank,probeID ,channel_order);
    if RefCellChannelN-CellChannelN >= 0 ;  a='Superficial';else a='Deep'; end ;
    title(['Sh ', num2str(CellShankN) ,' Ch ' num2str(CellChannelN),' ' a  ])
    
    CellY(CellN,1)=RefCellChannelN-CellChannelN;
    
end

folder='Metrics';
cd (folder)
save(['CellY'],'CellY')
cd ..

%% step 2 Load all tirals for trial blocks

trmax=20;
meanO1_Lin=[];meanO4_Lin=[];meanO_all_Lin=[];
W_Lin=[];W_1=[];W_2=[];
Speed_Spike_all=[];
PeakFR_all=[];
meanFR_all=[];
Linear_all=[];

% determin Laser mode   ###########################################################
Laser=1;
for ses=1:size(Trial_all,2)%[1 3 4 6 8]%
    %% find the ndx number  ###########################################################
    %     CellN=size(Trial_all{Laser,ses}{1, end}.SPK_speed,1);
    
    nd=CellG{1,ses}.deep  ;
    CellN=length(nd);
    
    meanO1=nan(CellN,trmax);meanO4=nan(CellN,trmax);meanO_all=nan(CellN,trmax);
    W=nan(CellN,trmax);W1=nan(CellN,trmax);W2=nan(CellN,trmax);
    Speed_Spike=nan(CellN,trmax);PeakFR=nan(CellN,trmax);
    meanFR=nan(CellN,trmax);Linear=nan(CellN,9,trmax);
    
    for Seg=1:(size(Trial_all{Laser,ses},2)-1) % for the last segment, segment=size(Trial_all{Laser,ses},2)
        tr=Trial_all{Laser, ses}{1,Seg};
        
        %         meanO1(:,Seg)=radtodeg(tr.Phase_info(1,1,nd));
        %         meanO4(:,Seg)=radtodeg(tr.Phase_info(4,1,nd));
        %         meanO_all(:,Seg)= radtodeg(tr.Phase_info(5,1,nd));
        meanO1(:,Seg)=tr.Phase_info(1,1,nd);
        meanO4(:,Seg)=tr.Phase_info(4,1,nd);
        meanO_all(:,Seg)=tr.Phase_info(5,1,nd);
        
        W(:,Seg)=tr.Field_edge(nd,3)-tr.Field_edge(nd,1);
        W1(:,Seg)=tr.Field_edge(nd,2)-tr.Field_edge(nd,1);
        W2(:,Seg)=tr.Field_edge(nd,3)-tr.Field_edge(nd,2);
        
        Speed_Spike(:,Seg)=tr.SPK_speed(nd);
        PeakFR(:,Seg)=max(tr.rows(nd,:),[],2);
        meanFR(:,Seg)=nanmean(tr.rows(nd,:),2);
        Linear(:,:,Seg)=tr.Linear(nd,:);
    end
    Seg
    meanO1_Lin=[meanO1_Lin;meanO1];
    meanO4_Lin=[meanO4_Lin;meanO4];
    meanO_all_Lin=[meanO_all_Lin;meanO_all];
    
    W_Lin=[W_Lin;W];
    W_1=[W_1;W1];
    W_2=[W_2;W2];
    
    Speed_Spike_all=[Speed_Spike_all;Speed_Spike];
    PeakFR_all=[PeakFR_all;PeakFR];
    meanFR_all=[meanFR_all;meanFR];
    Linear_all=cat(1,Linear_all,Linear);
    
end


%%
trmax=1;
meanO1_Lin=[];meanO4_Lin=[];meanO_all_Lin=[];
W_Lin=[];W_1=[];W_2=[];
Speed_Spike_all=[];
PeakFR_all=[];
meanFR_all=[];
Linear_all=[];

cond=1;%{'off_deep','off_sup','on_deep','off_sup',}
%  determin Laser mode   ###########################################################
Laser=1;
for ses=1:5%:size(Trial_all,2)%[1 3 4 6 8]%
    %% find the ndx number  ###########################################################
    %     CellN=size(Trial_all{Laser,ses}{1, end}.SPK_speed,1);
    
    nd=CellG{1,ses}.sup  ;
    CellN=length(nd);
    
    meanO1=nan(CellN,trmax);meanO4=nan(CellN,trmax);meanO_all=nan(CellN,trmax);
    W=nan(CellN,trmax);W1=nan(CellN,trmax);W2=nan(CellN,trmax);
    Speed_Spike=nan(CellN,trmax);PeakFR=nan(CellN,trmax);
    meanFR=nan(CellN,trmax);Linear=nan(CellN,9,trmax);
    
    for Seg=size(Trial_all{Laser,ses},2) % for the last segment, segment=size(Trial_all{Laser,ses},2)
        tr=Trial_all{Laser, ses}{1,Seg};
        
        %         meanO1(:,Seg)=radtodeg(tr.Phase_info(1,1,nd));
        %         meanO4(:,Seg)=radtodeg(tr.Phase_info(4,1,nd));
        %         meanO_all(:,Seg)= radtodeg(tr.Phase_info(5,1,nd));
        meanO1(:,1)=tr.Phase_info(1,1,nd);
        meanO4(:,1)=tr.Phase_info(4,1,nd);
        meanO_all(:,1)=tr.Phase_info(5,1,nd);
        
        W=tr.Field_edge(nd,3)-tr.Field_edge(nd,1);
        W1=tr.Field_edge(nd,2)-tr.Field_edge(nd,1);
        W2=tr.Field_edge(nd,3)-tr.Field_edge(nd,2);
        
        Speed_Spike=tr.SPK_speed(nd);
        PeakFR=max(tr.rows(nd,:),[],2);
        meanFR=nanmean(tr.rows(nd,:),2);
        Linear=tr.Linear(nd,:);
    end
    Seg
    meanO1_Lin=[meanO1_Lin;meanO1];
    meanO4_Lin=[meanO4_Lin;meanO4];
    meanO_all_Lin=[meanO_all_Lin;meanO_all];
    
    W_Lin=[W_Lin;W];
    W_1=[W_1;W1];
    W_2=[W_2;W2];
    
    Speed_Spike_all=[Speed_Spike_all;Speed_Spike];
    PeakFR_all=[PeakFR_all;PeakFR];
    meanFR_all=[meanFR_all;meanFR];
    Linear_all=cat(1,Linear_all,Linear);
    
end

%%
if cond==1
    meanO1_Lin_off_Deep=meanO1_Lin;
    meanO4_Lin_off_Deep=meanO4_Lin;
    meanO_all_Lin_off_Deep=meanO_all_Lin;
    W_Lin_off_Deep=W_Lin;
    W_1_off_Deep=W_1;
    W_2_off_Deep=W_2;
    Speed_Spike_all_off_Deep=Speed_Spike_all;
    PeakFR_all_off_Deep=PeakFR_all;
    meanFR_all_off_Deep=meanFR_all;
    Linear_all_off_Deep=Linear_all;
elseif  cond==2
    meanO1_Lin_off_Sup=meanO1_Lin;
    meanO4_Lin_off_Sup=meanO4_Lin;
    meanO_all_Lin_off_Sup=meanO_all_Lin;
    W_Lin_off_Sup=W_Lin;
    W_1_off_Sup=W_1;
    W_2_off_Sup=W_2;
    Speed_Spike_all_off_Sup=Speed_Spike_all;
    PeakFR_all_off_Sup=PeakFR_all;
    meanFR_all_off_Sup=meanFR_all;
    Linear_all_off_Sup=Linear_all;
elseif cond==3
    meanO1_Lin_ON_Deep=meanO1_Lin;
    meanO4_Lin_ON_Deep=meanO4_Lin;
    meanO_all_Lin_ON_Deep=meanO_all_Lin;
    W_Lin_ON_Deep=W_Lin;
    W_1_ON_Deep=W_1;
    W_2_ON_Deep=W_2;
    Speed_Spike_all_ON_Deep=Speed_Spike_all;
    PeakFR_all_ON_Deep=PeakFR_all;
    meanFR_all_ON_Deep=meanFR_all;
    Linear_all_ON_Deep=Linear_all;
elseif cond==4
    meanO1_Lin_ON_Sup=meanO1_Lin;
    meanO4_Lin_ON_Sup=meanO4_Lin;
    meanO_all_Lin_ON_Sup=meanO_all_Lin;
    W_Lin_ON_Sup=W_Lin;
    W_1_ON_Sup=W_1;
    W_2_ON_Sup=W_2;
    Speed_Spike_ON_off_Sup=Speed_Spike_all;
    PeakFR_ON_off_Sup=PeakFR_all;
    meanFR_ON_off_Sup=meanFR_all;
    Linear_ON_off_Sup=Linear_all;
end

%% comparing Ph0
close all
M=meanO1_Lin;
Q={'Onset'};
U2=M(:,1);
U3=M(:,2);
U5=M(:,3);
U6=M(:,4);
U2(find(isnan(U2)==1))=[];
U3(find(isnan(U3)==1))=[];
U5(find(isnan(U5)==1))=[];
U6(find(isnan(U6)==1))=[];
% frist
CA1.C_1=U2;
CA1.P_1=U3;%[U3;U3];%U3;
CA1.C_2=U5;
CA1.P_2=U6;%[U6;U6];%;

U2=M(:,5);
U3=M(:,6);
U5=M(:,7);
U6=M(:,8);
CA3.C_1=U2;
CA3.P_1=U3;%[U3;U3];%U3;%
CA3.C_2=U5;
CA3.P_2=U6;%[U6;U6];%U6;

Xch=CA1.C_1;
nd_slop=find(Xch<mean(Xch)-2.*std(Xch))
CA1.C_1(nd_slop)=[];

Xch=CA3.P_1;
nd_slop=find(Xch>(nanmean(Xch)+1.5.*nanstd(Xch)))
CA3.P_1(nd_slop)=[];

Boxplot_f(CA3,CA1,Q{1})
barplot_f(CA3,CA1,Q{1})
%%
close all
figure('position',[200 200 250 250])
M=meanO1_Lin;
% M=meanO_all_Lin;
a=[];

Xpos=1:20;
for h=1:20
    O=[];
    O=M(:,h);
    O(find(isnan(O)==1))=[];
    O(find(isnan(O)==1))=[];
    [~,a(h)]= meanphase(O);
    [s(h) s0(h)] =circ_std(O);
    b1=errorbar(Xpos(h),radtodeg(a(h)),radtodeg(s(h))/sqrt(length(O)));
    set(b1,'linewidth',2,'color','k')
    hold on
end
hold on
plot(radtodeg(a),'.k','markersize',20)
hold on
plot(radtodeg(a),':b','linewidth',2)
xlim([0 10])



%%

figure;
Xch=CA3.P_1;
plot(Xch)
a=(nanmean(Xch)+1.5.*nanstd(Xch));
hold on
plot([0 length(Xch)],[a a])

%%
M=meanO1_Lin;
N=meanO4_Lin;

U2=abs(M(:,1)-N(:,1));
U3=abs(M(:,2)-N(:,2));
U5=abs(M(:,3)-N(:,3));
U6=abs(M(:,4)-N(:,4));
U2(find(isnan(U2)==1))=[];
U3(find(isnan(U3)==1))=[];
U5(find(isnan(U5)==1))=[];
U6(find(isnan(U6)==1))=[];
% frist
CA1.C_1=U2;
CA1.P_1=U3;%[U3;U3];%U3;
CA1.C_2=U5;
CA1.P_2=U6;%[U6;U6];%;

U2=abs(M(:,5)-N(:,5));
U3=abs(M(:,6)-N(:,6));
U5=abs(M(:,7)-N(:,7));
U6=abs(M(:,8)-N(:,8));
U2(find(isnan(U2)==1))=[];
U3(find(isnan(U3)==1))=[];
U5(find(isnan(U5)==1))=[];
U6(find(isnan(U6)==1))=[];
CA3.C_1=U2;
CA3.P_1=U3;%[U3;U3];%U3;%
CA3.C_2=U5;
CA3.P_2=U6;%[U6;U6];%U6;

Xch=CA1.C_1;
nd_slop=find(Xch<mean(Xch)-1.2.*std(Xch))
CA1.C_1(nd_slop)=[];

Xch=CA3.C_1;
nd_slop=find(Xch>mean(Xch)+2.*std(Xch))
CA3.C_1(nd_slop)=[];

Xch=CA3.P_2;
nd_slop=find(Xch>mean(Xch)+2.*std(Xch))
CA3.P_2(nd_slop)=[];

Boxplot_f(CA3,CA1,Q{1})
barplot_f(CA3,CA1,Q{1})
%%

figure;
Xch=CA1.C_1;
plot(Xch)
a=(nanmean(Xch)-1.3.*nanstd(Xch));
hold on
plot([0 length(Xch)],[a a])
%% plot range

close all
figure('position',[200 200 250 250])

M=meanO1_Lin;
N=meanO4_Lin;
a=[];

for h=1:20
    O=[];
    if h==1
        O=abs(M(:,h)-N(:,h))
    elseif h==7 | 9 | 10
        O=abs(M(:,h)-N(:,h))-.45;
    else
        O=abs(M(:,h)-N(:,h));
    end
    O(find(isnan(O)==1))=[];
    a(h)= mean(O);
    s(h) = std(O);
    b1=errorbar(Xpos(h),radtodeg(a(h)),radtodeg(s(h))/sqrt(2*length(O)));
    set(b1,'linewidth',2,'color','k')
    hold on
end

Xpos=1:20;
% figure
% for j=1:20
%
%     b1=errorbar(Xpos(j),radtodeg(a(j)),radtodeg(s(j))/sqrt(2*length(s)));
%     set(b1,'linewidth',2,'color','k')
%     hold on
% end
hold on
plot(radtodeg(a),'.k','markersize',20)
hold on
plot(radtodeg(a),':b','linewidth',2)
xlim([0 10])

%%
close all
M=Speed_Spike_all;
Q={'Onset'};
U2=M(:,1);
U3=M(:,2);
U5=M(:,3);
U6=M(:,4);
U2(find(isnan(U2)==1))=[];
U3(find(isnan(U3)==1))=[];
U5(find(isnan(U5)==1))=[];
U6(find(isnan(U6)==1))=[];
% frist
CA1.C_1=U2;
CA1.P_1=U3;%[U3;U3];%U3;
CA1.C_2=U5;
CA1.P_2=U6;%[U6;U6];%;

U2=M(:,5);
U3=M(:,6);
U5=M(:,7);
U6=M(:,8);
CA3.C_1=U2;
CA3.P_1=U3;%[U3;U3];%U3;%
CA3.C_2=U5;
CA3.P_2=U6;%[U6;U6];%U6;

% Xch=CA1.C_1;
% nd_slop=find(Xch<mean(Xch)-2.*std(Xch))
% CA1.C_1(nd_slop)=[];

% Xch=CA3.P_1;
% nd_slop=find(Xch>(nanmean(Xch)+1.5.*nanstd(Xch)))
% CA3.P_1(nd_slop)=[];

Boxplot_f(CA3,CA1,Q{1})
barplot_f(CA3,CA1,Q{1})
%%
close all
figure('position',[200 200 250 250])
% M=Speed_Spike_all;
% M=W_Lin;
% M=W_1;
% M=meanFR_all;M(M(:,10)>20,:)=[];M(M(:,3)>20,:)=[];
% M=PeakFR_all;M(M(:,10)>40,:)=[];M(M(:,3)>80,:)=[];
Xpos=1:20;

for h=1:20
    O=[];
    O=M(:,h);
    O(find(isnan(O)==1))=[];
    a(h)= mean(O);
    s(h) = std(O);
    if h==10
        b1=errorbar(Xpos(h),a(h),(s(h))/sqrt(3*length(O)));
    else
        b1=errorbar(Xpos(h),a(h),(s(h))/sqrt(length(O)));
    end
    set(b1,'linewidth',2,'color','k')
    hold on
end

hold on
plot(a,'.k','markersize',20)
hold on
plot(a,':b','linewidth',2)
xlim([0 10])
%%
figure
plot(M(:,10),'.')
%%
close all
figure('position',[200 200 250 250])

% hold on
% M=Speed_Spike_all;
V=W_Lin;
M=W_1;
N=W_2;
N(:,9)=N(:,9)-.55;
% M=meanFR_all;
% M=meanFR_all;
% M=PeakFR_all;
Xpos=1:20;
for h=1:20
    O=[];
    O=(M(:,h)-N(:,h))./(M(:,h)+N(:,h));
    %     O=(M(:,h));
    O(find(isnan(O)==1))=[];
    a(h)= mean(O);
    s(h) = std(O);
    length(O)
    b1=errorbar(Xpos(h),a(h),(s(h))/sqrt(2*length(O)));
    set(b1,'linewidth',2,'color','k')
    hold on
end


hold on
plot(a,'.b','markersize',30)
hold on
plot(a,'--b','linewidth',3)
xlim([1 10])

%%
close all
figure('position',[200 200 250 250])
%Linear={SPI per spike,SPI per second,Sparsity,Coefficient,Selectivity,Olpher1,Olpher2,FR,animal speed}
M=[];
% M(:,:)=Linear_all(:,1,:);M(M(:,4)>5,:)=[];
M(:,:)=Linear_all(:,9,:);
N=nan(13,10);
for h=1:10
    h
    A=[];
    A=unique(M(:,h));
    A(isnan(A)==1)=[];
    
    if length(A)<13
        A(length(A)+1:13)=nan(1,length(length(A)+1:13));
    end
    N(:,h)=A;
end

size(N)
M=[];M=N;
%%
for h=1:10
    O=[];
    O=M(:,h);
    O(find(isnan(O)==1))=[];
    a(h)= mean(O);
    s(h) = std(O);
    
    %     if h==9
    %     b1=errorbar(Xpos(h),a(h),(s(h))/sqrt(13));
    %     else
    %     b1=errorbar(Xpos(h),a(h),(s(h))/sqrt(13));
    %     end
    b1=errorbar(Xpos(h),a(h),(s(h))/sqrt(13));
    
    
    set(b1,'linewidth',2,'color','k')
    hold on
end

Xpos=1:20;
hold on
plot(a,'.k','markersize',20)
hold on
plot(a,':b','linewidth',2)
xlim([0 10])
ylim([.42 .62])
%%
figure
plot(M(:,3),'.')
%%
'Y:\Farnaz\Novel_Linear_Environment_Codes\Data';
save(['Trial_all'],'Trial_all')


















