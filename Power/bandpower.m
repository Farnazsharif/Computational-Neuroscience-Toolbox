close all
% clear
filename='FM05_1';
load([filename '.mat']);
[trial]=Trialfinder1(filename);
c=[6]
d=[8]
theta=[6,7,8];
% c=[27,34,60,80,100,150];
% d=[32,60,80,100,150,250];
frq1=1;
frq2=15;
y1=300000;
% noise=141;% for rest PowerT3(161,:)=[];
[restT,runT]=allpower('FM05_1',10,2,frq1,frq2);
T=runT;


% load(filename1); 
%% Making Trial vector from Time
m=0;
for i=1:length(T) 
II = find(behav.TXDTS(trial(:,1),1)< T(i,2)& T(i,2)<=behav.TXDTS(trial(:,2),1));
TF = isempty(II);
if TF==1
    m=m+1;
    ndx_del(m)=i;
    if i < floor(length(T)/2)
      I(i)=1;
    else
      I(i)=length(trial); 
    end
else
    I(i)=II;
end
end
%% speed
d1=1;
speed=x2v(behav.TXDTS(:,3)*d1,behav.TXDTS(:,1));
for k= 1:length(T)
 an=find(behav.TXDTS(:,1)==T(k,1));
 bn=find(behav.TXDTS(:,1)==T(k,2));
 v=(speed(an:bn));
 V(k)=mean(v);
end

%% noise to delet
V(noise)=[];
I(noise)=[];
%% if noise exist
load(filename1);
PowerT3(noise,:)=[];
%% a - plot whole power
load(filename1);
% -1) plain power
clf
U=1;
subplot('Position', [0.21, .085, .62, .9])
imagesc(frq1:frq2,I,abs(PowerT3)) %1- plain power
% imagesc(I,1:150,abs(matnorm(PowerT3,1))) % for rest and run
ylabel('Trial Number ','fontsize',12)
xlabel('Frequency (Hz) ','fontweight','bold')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
colorbar('Position',[0.073  0.64 0.02 0.3])
cbar_handle = findobj(get(gcf,'Children'),'Tag','Colorbar')
% set(get(cbar_handle,'title'),'string', 'Normalized Power','fontsize',10)
set(get(cbar_handle,'title'),'string', ' Power','fontsize',10)
hold on
subplot('Position', [0.84, .085, .08, .9]);
imagesc(1,I,V')
xlabel('Speed ','fontweight','bold')
set(gca,'xtick',[1],'box','on')
set(gca,'ytick',[],'box','off')
colorbar('Position',[0.95 0.085 0.02 0.4])
cbar_handle = findobj(get(gcf,'Children'),'Tag','Colorbar')
set(get(cbar_handle(1),'title'),'string','cm/Sec','fontsize',10)
print('-depsc',['whole Power_' num2str(frq2) '_' num2str(U)  ])


% -2)
clf
U=2;
subplot('Position', [0.21, .085, .62, .9])
imagesc(frq1:frq2,I,abs(matnorm(PowerT3,2))) %1- plain power
ylabel('Trial Number ','fontsize',12)
xlabel('Frequency (Hz) ','fontweight','bold')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
colorbar('Position',[0.073  0.64 0.02 0.3])
cbar_handle = findobj(get(gcf,'Children'),'Tag','Colorbar')
% set(get(cbar_handle,'title'),'string', 'Normalized Power','fontsize',10)
set(get(cbar_handle,'title'),'string', ' Power','fontsize',10)

hold on
subplot('Position', [0.84, .085, .08, .9]);
imagesc(1,I,V')
xlabel('Speed ','fontweight','bold')
set(gca,'xtick',[1],'box','on')
set(gca,'ytick',[],'box','off')
colorbar('Position',[0.95 0.085 0.02 0.4])
cbar_handle = findobj(get(gcf,'Children'),'Tag','Colorbar')
set(get(cbar_handle(1),'title'),'string','cm/Sec','fontsize',10)
print('-depsc',['whole Power_' num2str(frq2) '_' num2str(U)  ])

% -3)
clf
U=3;
subplot('Position', [0.21, .085, .62, .9])
imagesc(frq1:frq2,I,abs(matnorm(PowerT3,1))) %1- plain power
% imagesc(I,1:150,abs(matnorm(PowerT3,1))) % for rest and run
ylabel('Trial Number ','fontsize',12)
xlabel('Frequency (Hz) ','fontweight','bold')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
colorbar('Position',[0.073  0.64 0.02 0.3])
cbar_handle = findobj(get(gcf,'Children'),'Tag','Colorbar')
% set(get(cbar_handle,'title'),'string', 'Normalized Power','fontsize',10)
set(get(cbar_handle,'title'),'string', ' Power','fontsize',10)
hold on
subplot('Position', [0.84, .085, .08, .9]);
imagesc(1,I,V')
xlabel('Speed ','fontweight','bold')
set(gca,'xtick',[1],'box','on')
set(gca,'ytick',[],'box','off')
colorbar('Position',[0.95 0.085 0.02 0.4])
cbar_handle = findobj(get(gcf,'Children'),'Tag','Colorbar')
set(get(cbar_handle(1),'title'),'string','cm/Sec','fontsize',10)
print('-depsc',['whole Power_' num2str(frq2) '_' num2str(U)  ])

%% b) Power-speed
clf
load(filename1);    
for t=1:length(c)   
 [hAx,h1,h2] =plotyy(I,mean(PowerT3(:,c(t):d(t)),2),I,V')
   
set(h2,'LineStyle','--')
title([ num2str(c(t)) '-' num2str(d(t)) ' HZ Power - Speed tr=20'],'fontsize',12) 
ylabel(hAx(1),'Wavelet Power dB','fontsize',12) % left y-axis
ylabel(hAx(2),'Speed cm/s','fontsize',12')
xlabel('Trial Number','fontsize',12)
hold on 
x=find(trial(:,4)==2,1)*ones(y1,1);
y=1:1:y1;
plot(x,y,'--k','LineWidth',2)
hold on
x=find(trial(:,4)==3,1)*ones(y1,1);
plot(x,y,'--k','LineWidth',2)
set(gcf,'color','w');
print('-depsc',['Power_speed_' num2str(c(t)) '_' num2str(d(t)) ])
close 
end
%%
%clf
figure(3)
theta=7;
for t=1:length(theta)   
[hAx,h1,h2] =plotyy(I,PowerT3(:,theta(t)),I,V')
   
set(h2,'LineStyle','--')
title([ num2str(theta(t)) ' HZ Power - Speed tr=10'],'fontsize',12) 
ylabel(hAx(1),'Wavelet Power dB','fontsize',12) % left y-axis
ylabel(hAx(2),'Speed cm/s','fontsize',12')
xlabel('Trial Number','fontsize',12)
hold on 
xi=find(trial(:,4)==2,1)*ones(y1,1);
y=1:1:y1;
plot(xi,y,'--k','LineWidth',2)
hold on
xj=find(trial(:,4)==3,1)*ones(y1,1);
plot(xj,y,'--k','LineWidth',2)
% 
% hold on
% x2=1:1:122;
% y2=mean(PowerT3(1:122,7),1);
% plot(x2,y2,'.r','LineWidth',10)

set(gcf,'color','w');
print('-depsc',['Power_speed_' num2str(theta(t)) ])
% close 
end
%%
plot(PowerT3(1:122,7))
hold on
x2=1:1:122;
y2=mean(PowerT3(1:122,7),1);
plot(x2,y2,'.k','LineWidth',10)
tLtrial
% first right lick :tr=7

6:50
%peyda kardane inke che trial hayi ro hazf mikone


%% behnam
% ndxj=[];
% DD=[];
% FDD=[];
% DD=diff(I)==0;
% FDD=diff(DD);
% [~ ,ndxj]=find(FDD==1); % here we have found reapetead elements
% [~ ,ndxe]=find(FDD==-1);
% ndxj1=ndxj+1;
% ndxj2=ndxe+1;
% if FDD(1)==-1;
%     ndxj11(2:size(ndxj1)+1)=ndxj1;
%     ndxj11(1)=1;
%     ndxj22(2:size(ndxj2)+1)=ndxj2;
%     ndxj22(1)=2;
% end
%     
% ndxt=[];
% ndxt(:,1)=ndxj1;
% ndxt(:,2)=ndxj2;

% for i=1:length(a3)
% j2(a3(i)-(i-1),:)=[];
% end



