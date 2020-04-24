
clearvars;clc;
dirData = 'Y:\OptoMECLEC\';
animal  = {'OML18','OML18','OML18','OML18','OML19','OML19'};
sessions  = {'day1','day2','day4','day5','day2','day3'};
ses=5;
dirses = [dirData animal{ses} '\'  sessions{ses}];
cd(dirses);

%% Add Laser time
L=[];
L=LoadBinary('analogin.dat','nChannels',1,'Frequency',30000);

%%
L2=[];
L2=(L - 32768) * 0.0003125;
L2(L2>-10)=min(L2);

%%
L2(1:226557952)=min(L2);
figure;
plot(L2(1:100:end),'.')

%%
tic
Laser=bz_getAnalogPulses('data',L2,'manualThr',true);
toc
save('Laser','Laser')

%% Get trials

load([sessions{ses} '.tracking.behavior.mat'])
load([sessions{ses} '.DigitalIn.events.mat'])
load('Laser')

txvtl(:,1)=tracking.timestamps;
txvtl(:,2)=tracking.position.x;

%%
figure;plot(txvtl(:,2));
%%
%OML19\day3 

% nd=find([txvtl(:,1)>9550 & txvtl(:,1)<10020]);
% txvtl(nd,:)=[];
% digitalIn.timestampsOn{1, 8}(20)=[];
%%

for jj=1:size(Laser.intsPeriods,1)
    nd=[];
    nd=find(txvtl(:,1)>=Laser.intsPeriods(jj,1) & txvtl(:,1)<=Laser.intsPeriods(jj,2) );
    txvtl(nd,5)=1;
end

ndx_On=[];
for i=1:length(digitalIn.timestampsOn{1, 7})
    [~,ndx_On(i,1)]=min(abs(digitalIn.timestampsOn{1, 7}(i)-(tracking.timestamps)));
end

ndx_Off=[];
for i=1:length(digitalIn.timestampsOn{1, 8})
    [~,ndx_Off(i,1)]=min(abs(digitalIn.timestampsOn{1, 8}(i)-(tracking.timestamps)));
end

nd_noise=[];
nd_noise=find(txvtl(:,2)>max(tracking.position.x(ndx_Off)));
txvtl(nd_noise,2)=max(tracking.position.x(ndx_Off));
nd_noise=[];
nd_noise=find(txvtl(:,2)<min(tracking.position.x(ndx_On)));
txvtl(nd_noise,2)=min(tracking.position.x(ndx_On));

%% sensor ndx to be deletted
%OML18\day1
k=39;
ndx_Off(k)=[];
digitalIn.timestampsOn{1,8}(k)=[];
ndx_On(k+1)=[];
digitalIn.timestampsOn{1,7}(k+1)=[];

%%
load([sessions{ses} '.tracking.behavior.mat'])
% load('OML18day5_linear_part2_190921_124957.tracking.behavior.mat')
% load('digitalIn.events.mat')
% Time_shift=+11;
% Time_shift=0;
tracking.timestamps=tracking.timestamps+Time_shift;
ndx_On=[];
for i=1:length(digitalIn.timestampsOn{1, 7})
    [~,ndx_On(i,1)]=min(abs(digitalIn.timestampsOn{1, 7}(i)-(tracking.timestamps)));
end

ndx_Off=[];
for i=1:length(digitalIn.timestampsOn{1, 8})
    [~,ndx_Off(i,1)]=min(abs(digitalIn.timestampsOn{1, 8}(i)-(tracking.timestamps)));
end

figure
plot(tracking.timestamps,tracking.position.x)
ylim([-3 3])
hold on
plot(digitalIn.timestampsOn{1, 7},tracking.position.x(ndx_On),'ok')
hold on
plot(digitalIn.timestampsOn{1, 8},tracking.position.x(ndx_Off),'or')
hold on

nd=find(txvtl(:,5)==1);
size(nd)

plot(txvtl(:,1),txvtl(:,2),'b')
hold on
plot(txvtl(nd,1),txvtl(nd,2),'.')
% ndx_Off(32)=[];
%%
% save('OML18day5_linear_part2_190921_124957.tracking.behavior.mat','tracking')
%%

figure
plot(ndx_On(1:145)-ndx_Off(1+1:145+1))
%%

%OML18\day2 
% ndx_Off(end)=[];
% ndx_Off=[ndx_Off];
% ND_all=[ndx_Off ndx_On];

%%
% OML18\day4
ndx_Off(16)=[];

% OML19\day2 
ndx_Off=[1; ndx_Off];


% OML19\day3 
ND_all=[ndx_Off ndx_On];


%OML19\day2 

ndx_Off=[1; ndx_Off];

%%
ND_all=[];
ND_all=[ndx_Off ndx_On];
for i=1:size(ND_all)
    
    txvtl(ND_all(i,1):ND_all(i,2),4)=(2*i-1); % Laser Off
 
    if i~=size(ND_all)
        
        txvtl(ND_all(i,2)+1:ND_all(i+1)-1,4)=(2*i); % Laser On
        txvtl(ND_all(i,2)+1:ND_all(i+1)-1,5)=1;
    else
        txvtl(ND_all(i,2)+1:size(txvtl,1),4)=(2*i); % Laser On
    end
end



win=1;
efps=30000./120;
x=txvtl(:,2);%sessInfo.Position.OneDLocation;
y=txvtl(:,2);%zeros(length(x),1);
spd = sqrt(sum(diff([NaN,NaN;[x,y]],1,1).^2,2))*efps;
spd = smoothdata(spd,'gaussian',2*win*efps);
spd(1) = spd(2) - diff(spd(2:3));
txvtl(:,3)=spd;

%% Plot the behaviour figure
filename=[animal{ses} sessions{ses}];
% Matrix=txvtl;
Matrix=behav.txvtl;
Matrix(:,2)=Matrix(:,2)-min(Matrix(:,2));
Matrix(:,2)=Matrix(:,2)./2;
% close all
C=unique(Matrix(:,4));
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
%         Matrix(nd_tr,2)=max(txvtl(:,2))-XX(:,2);

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
%%
behav.txvtl=Matrix;
behav.txvtl(:,5)=0;

for jj=1:size(Laser.intsPeriods,1)
    nd=[];
    nd=find(behav.txvtl(:,1)>=Laser.intsPeriods(jj,1) & behav.txvtl(:,1)<=Laser.intsPeriods(jj,2) );
    behav.txvtl(nd,5)=1;
end

%% Save behaviour file
save([filename '.mat'],'behav')




%% Plot PFs

trN=5;speed_thr=0.04;
delta=0.3;smooth_rate=10;smooth_phase=10;
Rate_Matrix=map.xtvtls2;
Phase_Matrix=map.xtvtlph;
BinSize=(max(behav.txvtl(:,2))./max(xtvtls2(:,1)))
smooth=10;

%%
cellN=1;  
plot_PF_1D(Rate_Matrix,smooth_rate,cellN);

%%
trN=5;speed_thr=0.04;
delta=0.1;smooth_rate=1;smooth_phase=1;
Rate_Matrix=map.xtvtls2;
Phase_Matrix=map.xtvtlph;
BinSize=(max(behav.txvtl(:,2))./max(map.xtvtls2(:,1)));

TR=[];
TR=unique(behav.txvtl(:,4));
TR(TR==0)=[];
r0=[];TRi=[];TR_all=[];TR_odd=[];
TR_odd=find(mod(TR,2)==1);
TR_all=TR(TR_odd);Steps=2;% 1 D linear
r0=[1 trN];
Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
TRi=TR_all(1:Seg*trN);
TRi=reshape(TRi,[trN Seg])';
tr=TRi;

cellN=1;  
plot_PF_1D(Rate_Matrix,smooth_rate,cellN,5);

%%
[Field_edge,rows,PF] = Feild_width_linearMaze(1,smooth_rate,delta,Rate_Matrix,TR_Off,5);
M(:,:)=PF(1,:,:);
figure
imagesc(matnorm(M,2))
figure
imagesc(M)

%%
M=[];
M=tiral{1, 12}.In_Spike{1, 1};
figure
plot(M(:,2),M(:,6),'.')
%%
TR=[];
TR=unique(behav.txvtl(:,4));
TR(TR==0)=[];
r0=[];TRi=[];TR_all=[];TR_even=[];
TR_even=find(mod(TR,2)==0);
TR_On=TR(TR_even);

TR=[];
TR=unique(behav.txvtl(:,4));
TR(TR==0)=[];
r0=[];TRi=[];TR_all=[];TR_Odd=[];
TR_Off=find(mod(TR,2)==1);
TR_Off=TR(TR_Off);

%%
for h = 1:size(spikes.times,2)

cellN = h;
TXVtlPh=[];
TXVtlPh=spikes.maze{1,cellN};
TXVtlPh(isnan(TXVtlPh(:,2))==1,:)=[];
TXVtlPh(find(TXVtlPh(:,3)<speed_thr),:)=[];%.04
TXVtlPh(TXVtlPh(:,4)==(0),:)=[];

tr=TR_On;
        trials_port=tr;
        nd_int=[];
        for i=1:length(trials_port)
            nd_int=[nd_int;find(TXVtlPh(:,4)==trials_port(i))];
        end

      

subplot(2,1,1)
plot([TXVtlPh(nd_int,2);TXVtlPh(nd_int,2)],[TXVtlPh(nd_int,6)+2*pi;TXVtlPh(nd_int,6)+4*pi],'.k','markersize',5)
    title(['Cell# = ' num2str(h) ' G=' num2str(spikes.UID(h)) '   ON'])  
xlim([0 1.2])
tr=TR_Off;
        trials_port=tr;
        nd_int=[];
        for i=1:length(trials_port)
            nd_int=[nd_int;find(TXVtlPh(:,4)==trials_port(i))];
        end

      

subplot(2,1,2)
plot([TXVtlPh(nd_int,2);TXVtlPh(nd_int,2)],[TXVtlPh(nd_int,6)+2*pi;TXVtlPh(nd_int,6)+4*pi],'.k','markersize',5)
xlim([0 1.2])
    title(['Cell# = ' num2str(h) ' G=' num2str(spikes.UID(h)) '   Off'])  

    cd PFs
    print('-djpeg',['Phasecell_' num2str(h)]) 
    cd ..
    
    
end

%%

hold on
for ii=2
    far=slope.*t+phi0+360*(ii-1);
    plot(t+min(M(:,2)),far,'Color','b','linewidth',2)
    hold on
end
hold on
binrange = [0:2*pi/100:4*pi];
B2=3;
x = binrange;
y1 = -cos(x)/10+B2;
x=radtodeg(binrange);
plot([y1],[x+B2],'r','linewidth',2)














