

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
%% ploting
% normaxe==1 normalize the columns
% normaxe==2 normalize the rows
figure(1)
imagesc(F1,I,matnorm(PowerT1,1)')

 xlim([ 1 80])
set(gca,'XMinorTick','on','YMinorTick','on')

figure(2)
scale=50;
B=imresize(PowerT1,scale,'lanczos3');
Ire1=imresize(I,[1 length(I)*scale],'lanczos3');
Fre1=imresize(F1,[1 length(F1)*scale],'lanczos3');
imagesc(Fre1,Ire1,abs(matnorm(B,1))');
% imagesc(F2,I2,abs((B))');
set(gca,'XMinorTick','on','YMinorTick','on')
 xlim([ 1 80])
%%
figure(1)
imagesc(F2,I,matnorm(PowerT2,1)')
xlim([ 1 80])

ylabel('Trial Number ','fontweight','bold')
xlabel('Frequency (Hz) ','fontweight','bold')
set(gcf,'color','w');
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


