clear
clc
clear
filename='FM05_1';
load([filename '.mat']);
[trial]=Trialfinder1(filename);
frq1=1;
frq2=150;
durationthreshold=2;
speedthreshold=10;
% [restT,runT]=allpower('FM05_1',speedthreshold,durationthreshold,frq1,frq2);
[restT1,runT1]=findrunrestT(behav.TXDTS(:,1),behav.TXDTS(:,3),speedthreshold,durationthreshold,1);
T=restT;
%% delet noise
noise=161;% rest
% noise=[];
T(noise,:)=[];
size(T)
filename1='PowerT3_n_150_rest';
load([filename1 '.mat']);
size(PowerT3)
PowerT3(noise,:)=[];
size(PowerT3)
%% Making Trial vector from Time
m=0;
ndx_del=[];
I=[];
II=[];
for i=1:length(T) 
II = find(behav.TXDTS(trial(:,1),1)< T(i,2)& T(i,2)<=behav.TXDTS(trial(:,2),1));
TF = isempty(II);
if TF==1
    m=m+1;
    ndx_del(m)=i;
    I(i)=NaN;
else
    I(i)=II;
end
end
%this section delet the out band endices

for i=1:length( ndx_del)
I( ndx_del(i)-(i-1))=[];
end
size(I)

%% shrink the size of the power and speed matrix
ndxj=[];
bl=[];
sndx=[];
ndxj1=[];
ndxj2=[];
[~ ,ndxj]=find(diff(I)==0); % here we have found reapetead elements
for i= 1:length(ndxj)
    sndx=find(I==I(ndxj(i)));
    %ndxj1(i)=ndxj(i); % or ndxj1(i)=sndx(1);
    ndxj1(i)=sndx(1);
    ndxj2(i)=sndx(end);
    bl=cat(2,bl,sndx);
end
ndxt=[];
ndxt(:,1)=ndxj1;
ndxt(:,2)=ndxj2;
%% delet same rows
ndxt1=[];
[C,ia,ic] = unique(ndxt(:,1));
ndxt1(:,1)=C;
ndxt1(:,2)=ndxt(ia,2);
%% Making related Power
% filename1='PowerT3_n_15';
% load([filename1 '.mat']);
size(PowerT3)
PowerT3( ndx_del,:)=[];
size(PowerT3)

%% shrink power to corelat with lick
ndxj4=[];
theta=7;
jt=[];
jt=PowerT3(:,:);
size(jt)

for i=1:length(ndxt1)
ndxj4(i,:)= mean(jt(ndxt1(i,1):ndxt1(i,2),:),1);
end
size(jt)

jt(bl,:)=0;
size(jt)

for  i=1:length(ndxt1)
jt(ndxt1(i,1),:)=ndxj4(i,:);
end

xx=find(jt(:,1)==0);
jt(xx,:)=[];
size(jt)
%% speed
V=[];
d1=1;
speed=x2v(behav.TXDTS(:,3)*d1,behav.TXDTS(:,1));
for k= 1:length(T)
 an=find(behav.TXDTS(:,1)==T(k,1));
 bn=find(behav.TXDTS(:,1)==T(k,2));
 v=(speed(an:bn));
 V(k)=mean(v);
end
V=V';
size(V)
V(ndx_del)=[];

size(V)

% apply to speed
ndxj4=[];
for i=1:length(ndxt)
ndxj4(i)= mean(V(ndxt(i,1):ndxt(i,2)),1);
end
V(bl)=0;
for  i=1:length(ndxt)
V(ndxt(i,1))=ndxj4(i);
end
V(find(V==0))=[];
size(V)
%% calculate Lick numbers
xbinNumber=100;
[li,Lt,lc]= licknumber(filename,xbinNumber);
%%
pow=[j1 unique(I)'];

lcpow=(lc(pow(:,2)));

powli=[j1 lcpow];

%%
[R,P] = corrcoef(j1,lc(pow(:,2)))
[R1,P1] = corrcoef(j1,V)
[R2,P2] = corrcoef(V,lc(pow(:,2)))

% [R1,P1] = corrcoef(powli(1:118,:))
% [R1,P1] = corrcoef(j1,li(pow(:,2)))
% [R2,P2] = corrcoef(j1,Lt(pow(:,2)))

%%
plotyy(pow(:,2),j1,pow(:,2),lcpow)
%%
[B1,I1] = sort(j1);
Spowli=[B1 lc(I1)];
%%
subplot(3,1,1)
plot(j1)
title(['Theta Power=' num2str(theta) 'HZ  Speed tr=' num2str(Speedthreshold)  ],'fontsize',12)
dim = [.15 .6 .3 .3];
str = [ 'Theta-lick-Cor=' num2str(R(1,2)) ' ; '  'Theta-Speed-Cor=' num2str(R1(1,2))];
annotation('textbox',dim,'String',str,'FitBoxToText','on');

hold on
y1=200000;
xtr1=find(trial(:,4)==1, 1, 'last');
[s1 s2]=min(abs(xtr1-unique(I)))
xi=s2*ones(y1,1);
y=1:1:y1;
plot(xi,y,'--k','LineWidth',2)
hold on

xtr2=find(trial(:,4)==2, 1,'last');
[s1 s2]=min(abs(xtr2-unique(I)));
xj=s2*ones(y1,1);
plot(xj,y,'--k','LineWidth',2)


subplot(3,1,2)
plot(lcpow)
title('Correct Lick percent')

hold on
y1=80;
xtr1=find(trial(:,4)==1, 1, 'last');
[s1 s2]=min(abs(xtr1-unique(I)))
xi=s2*ones(y1,1);
y=1:1:y1;
plot(xi,y,'--k','LineWidth',2)
hold on

xtr2=find(trial(:,4)==2, 1,'last');
[s1 s2]=min(abs(xtr2-unique(I)));
xj=s2*ones(y1,1);
plot(xj,y,'--k','LineWidth',2)


subplot(3,1,3)
plot(V)
title([ 'Speed tr=' num2str(Speedthreshold)  ],'fontsize',12)
dim = [.15 0 .3 .3];
str = [ 'Speed-lick-Cor=' num2str(R2(1,2)) ];
annotation('textbox',dim,'String',str,'FitBoxToText','on');hold on
y1=40;
xtr1=find(trial(:,4)==1, 1, 'last');
[s1 s2]=min(abs(xtr1-unique(I)))
xi=s2*ones(y1,1);
y=1:1:y1;
plot(xi,y,'--k','LineWidth',2)
hold on

xtr2=find(trial(:,4)==2, 1,'last');
[s1 s2]=min(abs(xtr2-unique(I)));
xj=s2*ones(y1,1);
plot(xj,y,'--k','LineWidth',2)
set(gcf,'color','w');