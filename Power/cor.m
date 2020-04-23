%% correlation for theta matrix only

% noise=141 % trial
% noise=161 % rest
% 
% V(noise)=[];
% I(noise)=[];
% size(V)
% size(I)

%%
% tr1=[1,1,134]; % for trials 
% tr2=[143,129,143];
% tr1=[2,2,134]; % for run I start from 2
% tr2=[144,129,143];
I(1)
tr1=[1,1,134]; % for rest I start from 2
tr2=[144,129,143];

%%
figure
% filename1='PowerT3_n_15_run.mat';
% load(filename1);
% PowerT3(noise,:)=[];  
for i=1%:length(W);
   w1=1;
   w2=128;
% w1=min(find(I==tr1(i)));
% w2=max(find(I==tr2(i)));
PV=[PowerT3(w1:w2,:) V(w1:w2)];
[R,P] = corrcoef(PV);
imagesc(R)
title(['Trials ' num2str(tr1(i)) ' - ' num2str(tr2(i)) '   Frequency 1-15-ret'] ,'fontsize',18)
set(gcf,'color','w');
colorbar
% print('-dtiff',['Correlation_1-15_' num2str(tr1(i)) '_' num2str(tr2(i))])
end

%%
% pp=PowerT3;
% pp(135,:)=[];
% pp(138,:)=[];
% pp(139,:)=[];
% sor=V;
% sor(:,135)=[];
% sor(:,138)=[];
% sor(:,139)=[];

%%
close all

% W=[1 6;7 45;46 52;53 94;95 101;102 128];
W=[1 6;7 36;37 45;46 52;53 86;87 94;95 101;102 121;122 128;129 I(end);1 128;1 I(end)];
ndx=[];1
for i=1: length(W)
   a= find(unique(I)>=W(i,1), 1, 'first');
   b= find(unique(I)>W(i,2), 1, 'first')-1;
   if isempty(b)==1
      b=length(unique(I));
   end
   ndx(i,1)=a;
   ndx(i,2)=b;
end
ndx
aaa=unique(I)';
% L=lc(aaa);
L=Lt(aaa);
close all

for i=1:length(W)
VV=[];
PVf=[];
Rf=[];
LL=[];
% VV=repmat(V(ndx(i,1):ndx(i,2)),1,20);
LL=repmat(L(ndx(i,1):ndx(i,2)),1,20);
PVf=[jt(ndx(i,1):ndx(i,2),:) LL];
[Rf,Pf] = corrcoef(PVf);
figure
imagesc(Rf)
title(['Rest ' num2str(W(i,1)) ' - ' num2str(W(i,2)) 'Freq 1-150-Lick#'] ,'fontsize',18)
set(gcf,'color','w');
colorbar
figure
imagesc(Pf)
title(['Rest ' num2str(W(i,1)) ' - ' num2str(W(i,2)) 'Frequ 1-150-Lick#'] ,'fontsize',18)
set(gcf,'color','w');
colorbar
% print('-dtiff',['Correlation-rest-1-150' num2str(tr1(i)) '_' num2str(tr2(i))])
end

%%
% tr1=[1,14,46,96,1,1,134,];
% tr2=[13,45,95,133,129,143,143];
% tr1=[1,14,47,97,1,1,134,];%rest
% tr2=[13,45,95,133,129,143,144];
% tr1=[2,14,47,97,2,2,134,];%run
% tr2=[13,45,95,133,129,143,143];




















