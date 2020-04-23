function plot_selectT(filename,Celln,smooth,Fr_max,tr_Select)
%plot_selectT('DJ55_S2_3T_known',17,10,Fr_max,tr_Select)

%%
%put the information for manipulation
% Qx_beg=[110,NaN,10];
% Qx_end=[130,NaN,25];
% Set=[1 2 3 4];
% Se=[1 11 31 52 tr(end)]



% filename='DJ55_S2_1T'
% Celln=107;
% smooth=10;
% Fr_max=max(rr);
% tr_Select=[60 80];

load([filename '.mat'])
load(['BeltInfo.mat'])

nPFbin=100;
% disp([  'initial trial = ' num2str(xttsc(1,3))   '        last trial = ' num2str(xttsc(end,3))  ]);
% tr_Select=input('tr_beg ,tr_end = ');
xttsc=xttsc(find(xttsc(:,3)==tr_Select(1) & xttsc(:,1)==1):find(xttsc(:,3)==tr_Select(2) & xttsc(:,1)==100),:);
smoothT=smooth1D(xttsc(:,2),smooth,1);
smoothC=smooth1D(xttsc(:,Celln+4),smooth,1);
ncell=length(smoothC(1,:));
rate=smoothC./repmat(smoothT,1,ncell);
xtsr=[xttsc(:,[1 3 4]) rate]; 
rr=xtsr(:,4)./Fr_max;
matC=reshape(rr,nPFbin,length(rr)/nPFbin)';
%  matC=matnorm(matC,2);
imagesc(matC.^(1))
% imagesc(matC)











% %% draw objects
% [y,x]=size(matC);
% hold on
%     yspan=y;
%     for jj=1:length(BeltInfo.object_bgn)
%         plot([BeltInfo.object_bgn(jj) BeltInfo.object_bgn(jj) BeltInfo.object_end(jj) BeltInfo.object_end(jj) BeltInfo.object_bgn(jj)]/BeltInfo.Length*nPFbin,[yspan+1 -1 -1 yspan+1 yspan+1],BeltInfo.objectC{jj})
%     end  
% 
% %% ploting added cues
% hold on 
% 
% Se=[];
% tr=[];
% 
% Set=unique(behav.txlrts(:,6));
% for j=1:length(Set)
% tr_session=find(behav.txlrts(:,6)==Set(j));
% tr(j)=behav.txlrts(tr_session(1),5);
% end
% Se=[tr behav.txlrts(end,5)];
% 
% 
%  if length(Set)>1
%  for k=1:length(Set)
%      
% plot([0 x],[Se(k)-1 Se(k)-1],'color','w')
% plot([BeltInfo.Qx_beg(k) BeltInfo.Qx_beg(k)]/BeltInfo.Length*nPFbin,[Se(k)-1 Se(k+1)-1],'color','w')
% plot([BeltInfo.Qx_end(k) BeltInfo.Qx_end(k)]/BeltInfo.Length*nPFbin,[Se(k)-1 Se(k+1)-1],'color','w')
%  end
%  end
%  
% % set(gca,'XTickLabel',[])
% % set(gca,'Xtick',[0:20:100])
% % set(gca,'XTickLabel',[0:40:Beltinfo.Length])

