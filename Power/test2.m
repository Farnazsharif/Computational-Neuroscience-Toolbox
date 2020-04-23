%% power 
%close all
%1-make correlation btween time and trials, 2-make correlation for
%frequency
clc
clear
ch=115;
speedthreshold=5;
durationthreshold=2;
filename='FM05_1';
load([filename '.mat']);
eeg= readmulti([filename '.lfp'], 128,ch);
EEGsamplerate=(spkinfo.samplerate./25);
Teeg = (1:length(eeg))/EEGsamplerate;
[trial]=Trialfinder(filename);
% [restT,runT]=findrunrestT(behav.TXDTS(:,1),behav.TXDTS(:,3),speedthreshold,durationthreshold);
T=behav.runT;
%T=behav.restT;
%T=[behav.TXDTS(trial(:,1)) behav.TXDTS(trial(:,2))];
Ndx=[];
for j=1:length(T)
  [~, Ndx(j,1)]=min(abs(T(j,1)-Teeg));
  [~, Ndx(j,2)]=min(abs(T(j,2)-Teeg));
end
Ndx;
aa=Ndx(:,2)-Ndx(:,1);
%Senv = envelop(S)
%% Making Trial vector from Time
for i=1:length(T) 
I(i) = find(behav.TXDTS(trial(:,1),1)< T(i,2)& T(i,2)<=behav.TXDTS(trial(:,2),1));
end
I=I+5;
%% speed calculation

d = 0.5; %wheel increment ~ 0.5cm
speed=[];
for k= 1:length(T)
 an=find(behav.TXDTS(:,1)==T(k,1));
 bn=find(behav.TXDTS(:,1)==T(k,2));
 t=behav.TXDTS(an:bn,1); 
 x=behav.TXDTS(an:bn,3);
v=x2v(x*d,t,1);
speed=cat(1,v,speed);
V(k)=mean(v,1);
end
imagesc(1,I,V')

%% Method1 sacling based on frequency
power=[];
for g=1:5%length(Ndx)
L = floor(0.5*(Ndx(g,2)- Ndx(g,1)));
[Pxx F1]= pwelch(eeg(Ndx(g,1):Ndx(g,2)),L,[],[1:150],EEGsamplerate);
 power(:,g)=Pxx;
end
%% ploting
% normaxe==1 normalize the columns
% normaxe==2 normalize the rows
figure(1)
imagesc(F1,I,matnorm(power,1)')
xlim([ 4 10])
set(gca,'XMinorTick','on','YMinorTick','on')

figure(2)
scale=50;
B=imresize(power,scale,'lanczos3');
I2=imresize(I,[1 length(I)*scale],'lanczos3');
F2=imresize(F1,[1 length(F1)*scale],'lanczos3');
imagesc(F2,I2,abs(matnorm(B,1))');
set(gca,'XMinorTick','on','YMinorTick','on')
xlim([ 4 10])
%% Method 2 Scalin based on interpolation
% [~, mi]=min(aa)
% [~, ma]=max(aa)
% [Pxxma F]= pwelch(eeg(Ndx(ma,1):Ndx(ma,2)),L,[],[],EEGsamplerate);

%% Method 2 Scalin based on interpolation
% power2=[];
% scale2=4000;
% for g=1:length(Ndx)    
% % g=ma;
% L = floor(0.5*(Ndx(g,2)- Ndx(g,1)));
% [Pxx F]= pwelch(eeg(Ndx(g,1):Ndx(g,2)),L,[],[],EEGsamplerate);
% c=imresize(Pxx,[scale2 1],'lanczos3');
% power2(:,g)=c;
% end
% imagesc(matnorm(power2,1)')
% xlim([ 0 500])
%%
% F3=imresize(F,[1 scale2],'lanczos3');
% imagesc(F3,I,matnorm(power2,1)')
% set(gca,'XMinorTick','on','YMinorTick','on')
% %

% hold on
% plot ( matnorm(b,1) , 'r')
% xlim([ 0 800])
% hold on
%%

plot(aa./EEGsamplerate)
set(gcf,'color','w');
set(gca,'fontsize',18,'fontweight','bold')
% title('frequency ranges from 40 to 60','fontweight','bold')
ylabel('Running Window Duration (sec)','fontweight','bold')
xlabel('Runing Window #','fontweight','bold')
set(gcf,'color','w');

%% wrong power
tic
power3=[];

for g=1:5%length(Ndx)    
% g=ma;

[Pxx F]= pwelch(eeg(Ndx(g,1):Ndx(g,2)),floor(0.5*EEGsamplerate),[],[],EEGsamplerate);

power3(:,g)=Pxx;
end
toc
%%
figure(1)
imagesc(F,I,matnorm(power3,1)')
xlim([ 4 10])

ylabel('Trial Number ','fontweight','bold')
xlabel('Frequency (Hz) ','fontweight','bold')
set(gcf,'color','w');

figure(2)
scale=10;
B3=imresize(power3,scale,'lanczos3');
I3=imresize(I,[1 length(I)*scale],'lanczos3');
F3=imresize(F,[1 length(F)*scale],'lanczos3');
imagesc(F3,I3,abs(matnorm(B3,1))');
%imagesc(abs(matnorm(B3,1))');
xlim([ 4 10])

%% subplot power and speed
close all
clc

%subplot(1,2,1)
subplot('Position', [0.21, .085, .62, .9])
% imagesc(F2,I2,abs((B))');
 imagesc(F3,I3,abs((B3))');

% imagesc(F3,I3,abs(matnorm(B3,1))');
% imagesc(F2,I2,abs(matnorm(B,1))');

%imagesc(F1,I,matnorm(power,1)')
%imagesc(F,I,power3')
%xlim([ 0 80])
%title('power 20-40')
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
%[left bottom width height])
% left + width < 1 and that bottom + height < 1 (for the first subplot)
%% test plot of individual power
close all
hold on
for ss=1:5

plot(F,power3(:,ss),'color',[((ss)*9+6)/100  .5 1])
hold on
end
%% eeg concating

 lfp2=[];
for xx=1:2
    
 eegp=eeg(Ndx(xx,1):Ndx(xx,2));
 lfp2=cat(1,eegp,lfp2);
    
end
%%
hold on
[pcc F]= pwelch(cc,floor(0.5*EEGsamplerate),[],[],EEGsamplerate);

L = floor(0.5*(length(cc)));
[pcc2 F2]= pwelch(cc,L,[],[1:150],EEGsamplerate);
hold on
plot(F2,pcc2,'color',[.5  .5 .5])
%%
plot(F,pcc,'color',[.3  .2 .5])
hold on
plot(F2,pcc2,'color',[.5  .6 .1])
hold on
pmean=power3(:,1:13);
plot(F,mean(pmean,2),'r')
hold on
pmean=power(:,1:13);
plot(F1,mean(pmean,2),'b')
xlim([ 0 80])
set(gcf,'color','w');





