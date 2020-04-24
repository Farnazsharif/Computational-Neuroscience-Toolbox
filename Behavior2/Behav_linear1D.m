%% Creat behave file for linear track
clc
clear
filename= 'Achilles_10252013';
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
plot(C1)
figure
plot(C,'.')

win1=1.57;%1.57;%'Achilles_11012013';1.56;1.98%
win2=-1.57;%1.57;%-1.57%-1.57;1.98%
ndx_1=find(C1>win1);
ndx_2=find(C1<win2);

figure
plot(C1)
hold on
plot(ndx_1,C1(ndx_1),'or','markersize',5)
hold on
plot(ndx_2,C1(ndx_2),'ok','markersize',5)

ndx_1=ndx_1;
ndx_2=ndx_2;

figure
plot(C,'.')
hold on
plot(ndx_1,C(ndx_1),'or','markersize',5)
hold on
plot(ndx_2,C(ndx_2),'ok','markersize',5)

%%
thr=0.0148;
for i=1:length(ndx_2)
    alfa=[];
    alfa=find(Matrix(ndx_2(i):end,2)<thr);
    if isempty(alfa)==1
        alfa=nan;
    end
    ndx_2(i,2)=alfa(1)+ndx_2(i,1)-1;
end
ndx_2(isnan(ndx_2(:,2))==1,:)=[];


ND_even=[];A=[];
A=unique(ndx_2(:,2));
for  j=1:length(A)
    nd=find(ndx_2(:,2)==A(j));
    [a,b]=max(ndx_2(nd,1));
    ND_even(j,:)=[a A(j)];
end

close all
for j=1:length(ND_even)
    YY=[];nd_tr=[];XX=[];
    figure
    nd_tr=ND_even(j,1):ND_even(j,2);
    XX=Matrix(nd_tr,1:3);
    nd_interpolate1=find(XX(:,2) <= win1 );nd_interpolate2=find(XX(:,2) >= 0.05);
    nd_interpolate=nd_interpolate1(1):nd_interpolate2(end);
    YY=fillmissing(XX(nd_interpolate,2),'pchip');
    plot(XX(nd_interpolate,1),YY,'.k')
    Matrix(nd_tr(nd_interpolate),2)=YY;
    Matrix(nd_tr(nd_interpolate),4)=2*j;
    
%     VV=fillmissing(XX(nd_interpolate,3),'pchip');  
%     Matrix(nd_tr(nd_interpolate),3)=VV;
%     figure
%     plot(XX(nd_interpolate,1),VV*100,'.r')
%     hold on
%     plot(XX(nd_interpolate,1),XX(nd_interpolate,3)*100,'.')
    if isnan(XX(nd_interpolate,3))
        V_nan=[Vnan  j]
    end
end

tr=unique(Matrix(:,4));tr(tr==0)=[];
figure
for i=1:length(tr)
    XX=Matrix(Matrix(:,4)==tr(i),1:4);
    plot(XX(:,2),XX(:,4),'.')
    hold on
end


%%

thr=0.0148;
for i=1:length(ndx_1)
    alfa=[];
    alfa=find(Matrix(1:ndx_1(i),2)<thr);
    if isempty(alfa)==1
        alfa=nan;
    end
    ndx_1(i,2)=alfa(end);
end
ndx_1(isnan(ndx_1(:,2))==1,:)=[];
% ndx_1(1:3,:)=[];%'Cicero_09012014'

ND_odd=[];A=[];
A=unique(ndx_1(:,2));

for  j=1:length(A)
    nd=find(ndx_1(:,2)==A(j));
    [a,b]=max(ndx_1(nd,1));
    ND_odd(j,:)=[A(j) a ];
end

close all
for j=1:length(ND_odd)
    YY=[];nd_tr=[];XX=[];nd_interpolate1=[];nd_interpolate2=[];nd_interpolate=[];
    figure
    nd_tr=ND_odd(j,1):ND_odd(j,2);
    XX=Matrix(nd_tr,1:3);
    nd_interpolate1=find(XX(:,2) >= 0.04 );nd_interpolate2=find(XX(:,2) >= win1 );
    if isempty(nd_interpolate2)==1
       [~,nd_interpolate2]=max(XX(:,2)); 
    end
    nd_interpolate=nd_interpolate1(1):nd_interpolate2(1);
    YY=fillmissing(XX(nd_interpolate,2),'pchip');
    plot(XX(nd_interpolate,1),YY,'.k')
    Matrix(nd_tr(nd_interpolate),2)=YY;
    Matrix(nd_tr(nd_interpolate),4)=2*j-1;
    if isnan(XX(nd_interpolate,3))
        V_nan=[Vnan  j]
    end
  end

tr=unique(Matrix(:,4));tr(tr==0)=[];
figure
for i=1:length(tr)
    XX=Matrix(Matrix(:,4)==tr(i),1:4);
    plot(XX(:,2),XX(:,4),'.')
    hold on
end

%%
% velo
% Matrix=[];
% Matrix=behav.TXVt;

close all
nd_nantr=find(Matrix(:,4)==0)
Matrix(nd_nantr,2)=nan;

C=unique(Matrix(:,4));C(C==0)=[];
figure('position',[100 100 1800 400])
subplot(2,1,1)
plot(Matrix(:,1),Matrix(:,2),'.k')
title(filename,'fontsize',10)
hold on
for i=1:length(C)
    
    XX=[];mean_V=[];
    nd_tr=find(Matrix(:,4)==C(i));
    XX=Matrix(nd_tr,1:4);
    
    if mod(C(i),2)==0
       
        Matrix(nd_tr,2)=max(TXVT(:,2))-XX(:,2);
        subplot(2,1,1)
        plot(Matrix(nd_tr,1),Matrix(nd_tr,2),'.b')
        xlabel('time','fontsize',15)
        ylabel('Position','fontsize',15)
        hold on
        subplot(2,1,2)
        plot(mean(XX(:,1)),mean(XX(:,3))*100,'ob','linewidth',2)
        hold on
        xlabel('time','fontsize',15)
        ylabel('Speed','fontsize',15)
    else
        subplot(2,1,1)
        plot(Matrix(nd_tr,1),Matrix(nd_tr,2),'.r')
        xlabel('time','fontsize',15)
        ylabel('Position','fontsize',15)
        hold on
        subplot(2,1,2)
        plot(mean(XX(:,1)),mean(XX(:,3))*100,'or','linewidth',2)
        hold on
        xlabel('time','fontsize',15)
        ylabel('Speed','fontsize',15)
    end
        
end

figure
plot(Matrix(nd_nantr,2))

behav.TXVt=Matrix;
behav.TXVt_Phas=[];
behav.TXVt_rate=[];
save([filename '.mat'],'-append','behav')







