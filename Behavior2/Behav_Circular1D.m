%% Creat Behav file for Circular track
clear
close all
filename='Achilles_10252013';
load([filename '_sessInfo.mat'])

txvt(:,1)=sessInfo.Position.TimeStamps;
txvt(:,2)=sessInfo.Position.OneDLocation;

win=1;
efps=20000./512;
x=sessInfo.Position.TwoDLocation(:,1);%sessInfo.Position.OneDLocation;
y=sessInfo.Position.TwoDLocation(:,2);%zeros(length(x),1);
spd = sqrt(sum(diff([NaN,NaN;[x,y]],1,1).^2,2))*efps;
spd = smoothdata(spd,'gaussian',2*win*efps);
spd(1) = spd(2) - diff(spd(2:3));
txvt(:,3)=spd;

MazeEpoch=find(txvt(:,1)>=sessInfo.Epochs.MazeEpoch(:,1) & txvt(:,1)<=sessInfo.Epochs.MazeEpoch(:,2));
txvt=txvt(MazeEpoch,:);

TXVT=[];
T=txvt(1,1):0.01:txvt(end,1);
X=interp1(txvt(:,1),txvt(:,2),T);
V=interp1(txvt(:,1),txvt(:,3),T);
TXVT=[T' X' V'];
Matrix=TXVT;

%%
close all
C=[];
C=Matrix(:,2);
C(isnan(C)==1)=0;
C1=diff(C);
C1(length(C))=C1(end);
figure
plot(C1,'.')
figure
plot(C,'.')

win1=2.7;%
win2=-2.7;%
ndx_1=find(C1>win1);
ndx_2=find(C1<win2);

figure
plot(C1)
hold on
plot(ndx_1,C1(ndx_1),'or','markersize',5)
hold on
plot(ndx_2,C1(ndx_2),'ok','markersize',5)

ndx_1=ndx_1+1;
ndx_2=ndx_2+1;

figure
plot(C,'.')
hold on
plot(ndx_1,C(ndx_1),'or','markersize',5)
hold on
plot(ndx_2,C(ndx_2),'ok','markersize',5)

%%

thr_low=0.05;
for i=1:length(ndx_1)
    alfa=[];
    alfa=find(Matrix(ndx_1(i):end,2)<thr_low);
    if isempty(alfa)==1
        alfa=nan;
    end
    ndx_1(i,2)=alfa(1)+ndx_1(i,1)-1;
end
ndx_1(isnan(ndx_1(:,2))==1,:)=[];

ND_all=[];A=[];
A=unique(ndx_1(:,2));
for  j=1:length(A)
    nd=find(ndx_1(:,2)==A(j));
    [a,b]=max(ndx_1(nd,1));
    ND_all(j,:)=[a A(j)];
end

XX=[];B=[];win=100;
for i=1:length(ND_all)
    
    XX=[];
    XX=Matrix(ND_all(i,1):ND_all(i,2),2);
    XX(isnan(XX)==1)=[];
    B(i,1)=length(XX);
    if length(XX) > win
        B(i,2)=mean(XX(length(XX)-100:length(XX)));
    else
        B(i,2)=-1;
    end
end

ND_del=find(B(:,2)==-1 | B(:,2)>1.8);
ND_all(ND_del,:)=[];

for j=1:length(ND_all)
    Matrix(ND_all(j,1):ND_all(j,2),4)=j;
end


figure
for i=1:length(ND_all)
    XX=[];
    XX=Matrix(Matrix(:,4)==i,1:4);
    plot(XX(:,2),XX(:,4),'.b')
    hold on
end

close all
for i=1:length(ND_all)
    figure
    nd_tr=find(Matrix(:,4)==i);
    XX=Matrix(nd_tr,1:4);
    nd_interpolate1=find(XX(:,2) >= 0 );%nd_interpolate2=find(XX(:,2) <= 3);
    nd_interpolate=nd_interpolate1(1):nd_interpolate1(end);
    YY=fillmissing(XX(nd_interpolate,2),'pchip');
    plot(XX(nd_interpolate,1),YY,'.k')
    hold on
    plot(XX(:,1),XX(:,2),'.')
    Matrix(nd_tr(nd_interpolate),2)=YY;
    
    hold on
    plot(XX(:,1),XX(:,2),'.b')
    Matrix(nd_tr(nd_interpolate),2)=YY;
    
    if isnan(XX(nd_interpolate,3))
        V_nan=[Vnan  j]
    end
    
%     hold on
%     plot(XX(:,1),flip(XX(:,2),1),'.r')
%     hold on
%     plot(XX(:,1),flip(XX(:,2),1),'.k')
%     hold on
%     plot(XX(:,1),win1-XX(:,2),'.r')
end

%%

thr_high=2.8;
for i=1:length(ndx_2)
    alfa=[];
    alfa=find(Matrix(ndx_2(i):end,2)>thr_high);
    if isempty(alfa)==1
        alfa=nan;
    end
    ndx_2(i,2)=alfa(1)+ndx_2(i,1)-1;
end
ndx_2(isnan(ndx_2(:,2))==1,:)=[];

ND_bad=[];A=[];
A=unique(ndx_2(:,2));
for  j=1:length(A)
    nd=find(ndx_2(:,2)==A(j));
    [a,b]=max(ndx_2(nd,1));
    ND_bad(j,:)=[a A(j)];
end

XX=[];B=[];win=100;
for i=1:length(ND_bad)
    
    XX=[];
    XX=Matrix(ND_bad(i,1):ND_bad(i,2),2);
    XX(isnan(XX)==1)=[];
    B(i,1)=length(XX);
    if length(XX) > win
        B(i,2)=mean(XX(length(XX)-100:length(XX)));
    else
        B(i,2)=-1;
    end
end
ND_del=[];
ND_del=find(B(:,2)==-1 | B(:,2) <1.5);
ND_bad(ND_del,:)=[];

for j=1:size(ND_bad,1)
    Matrix(ND_bad(j,1):ND_bad(j,2),4)=-j;
end

XX=[];
close all
figure
for i=1:size(ND_bad,1)
    XX=Matrix(Matrix(:,4)==-i,1:4);
    plot(XX(:,2),XX(:,4),'.')
    hold on
end

for i=1:size(ND_bad,1)
    XX=[];YY=[];
    figure
    nd_tr=find(Matrix(:,4)==-i);
    XX=Matrix(nd_tr,1:4);
    nd_interpolate1=find(XX(:,2) >= 0 );%nd_interpolate2=find(XX(:,2) <= 3);
    nd_interpolate=nd_interpolate1(1):nd_interpolate1(end);
    YY=fillmissing(XX(nd_interpolate,2),'pchip');
    plot(XX(nd_interpolate,1),YY,'.k')
    hold on
    plot(XX(:,1),XX(:,2),'.r')
    Matrix(nd_tr(nd_interpolate),2)=YY;
    if isnan(XX(nd_interpolate,3))
        V_nan=[Vnan  j]
    end
end


%%
% Matrix=[];
% Matrix=behav.TXVt;
Matrix((Matrix(:,3)<0),3)=0;
nd_nantr=find(Matrix(:,4)==0)
Matrix(nd_nantr,2)=nan;

close all
C=unique(Matrix(:,4));C(C==0)=[];
figure('position',[100 100 1800 400])
subplot(2,1,1)
plot(Matrix(:,1),Matrix(:,2),'.r')
title(filename,'fontsize',10)
hold on
for i=1:length(C)
    subplot(2,1,1)
    XX=[];mean_V=[];
    nd_tr=find(Matrix(:,4)==C(i));
    XX=Matrix(nd_tr,1:4);
    
    Matrix(nd_tr,2)=max(TXVT(:,2))-XX(:,2);
    if C(i)>0
        plot(Matrix(nd_tr,1),Matrix(nd_tr,2),'.b')
    else
        plot(Matrix(nd_tr,1),Matrix(nd_tr,2),'.k')
    end
    xlabel('time','fontsize',15)
    ylabel('Position','fontsize',15)
    hold on
    subplot(2,1,2)
    plot(nanmean(XX(:,1)),nanmean(XX(:,3))*100,'ob','linewidth',2)
    hold on
    xlabel('time','fontsize',15)
    ylabel('Speed','fontsize',15)
end



figure
plot(Matrix(nd_nantr,2))
%%
% filename=sessions{ses};
behav.TXVt=Matrix;
% behav.TXVt_Phas=[];
% behav.TXVt_rate=[];
save([filename '.mat'],'-append','behav')







