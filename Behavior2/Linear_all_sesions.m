
trmax=20;
meanO1_Lin=[];meanO4_Lin=[];meanO_all_Lin=[];
W_Lin=[];W_1=[];W_2=[];
Speed_Spike_all=[];
PeakFR_all=[];
meanFR_all=[];
Linear_all=[];
for ses=1:size(Trial_all,2)%[1 3 4 6 8]%
    CellG=[];
    CellG=Trial_all{2, ses};
    G=Trial_all{1, ses}{1, 1}.G; % First Segment
    %% find ndx of the cell groups in each session
    nd=[];
    for C=1:length(CellG)
        nd=[nd;find(G==CellG(C))];
    end
    
    %%
    meanO1=nan(length(CellG),trmax);meanO4=nan(length(CellG),trmax);meanO_all=nan(length(CellG),trmax);
    W=nan(length(CellG),trmax);W1=nan(length(CellG),trmax);W2=nan(length(CellG),trmax);
    Speed_Spike=nan(length(CellG),trmax);PeakFR=nan(length(CellG),trmax);
    meanFR=nan(length(CellG),trmax);Linear=nan(length(CellG),9,trmax);
    
    for Seg=1:(size(Trial_all{1,ses},2)-1)
        tr=Trial_all{1, ses}{1,Seg};
        
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


















