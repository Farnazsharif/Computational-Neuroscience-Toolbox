function plotBeltRew_chronic(filename)

% filename='SFA3_S1_TRD2';
% filename2='SFA3_S1_TRD2PlaceField_Total';
% smooth=10;
% Celln=26;

%%
load([filename '.mat'])
load([filename 'PlaceField.mat'])
Celln=1;
smooth=1;
%%
nPFbin=100;
smoothT=smooth1D(xttsc(:,2),smooth,1);
smoothC=smooth1D(xttsc(:,Celln+4),smooth,1);
ncell=length(smoothC(1,:));
rate=smoothC./repmat(smoothT,1,ncell);
xtsr=[xttsc(:,[1 3 4]) rate];
rr=xtsr(:,4);
matC=reshape(rr,nPFbin,length(rr)/nPFbin)';
%% draw objects

[y,x]=size(matC);

hold on
yspan=y;
for jj=1:length(Beltinfo.object_bgn)
    plot([Beltinfo.object_bgn(jj) Beltinfo.object_bgn(jj) Beltinfo.object_end(jj) Beltinfo.object_end(jj) Beltinfo.object_bgn(jj)],[yspan+1 -1 -1 yspan+1 yspan+1],Beltinfo.objectC{jj},'linewidth',2)
end
% xlim([0 101])
% ylim([-1.2 13.2])
% set(gca,'Color','black')
% set(gca,'YTickLabel',[])
% %     set(gca,'Color',[0.8,0.8,0.8])
%     set(gca,'XTickLabel',[])
%% plotting reward positions

hold on
Vacume=15;
Set=unique(behav.txlrts(:,6));

if length(Set)>1
    
    for k=1:length(Set)
        
        y_ndx=find(xtsr(:,3)==Set(k));
        y1=floor(y_ndx(1)/100);
        y2=floor(y_ndx(end)/100);

        
        tmp=behav.txlrts(find(behav.txlrts(:,6)==Set(k)),:);
        Rew_x=tmp(find(tmp(:,4)==1),2);
        Rew_x(1);
        
        plot([Rew_x(1) Rew_x(1) Rew_x(1)+Vacume Rew_x(1)+Vacume Rew_x(1)],[y2 y1 y1 y2 y2],'r','linewidth',2)
        hold on
    end
    
    
else
    
    
    Rew_ndx=find(behav.txlrts(:,4)==1);
    Rew_x=behav.txlrts(Rew_ndx,2);
    plot([Rew_x(1) Rew_x(1) Rew_x(1)+Vacume Rew_x(1)+Vacume Rew_x(1)],[yspan+1 -1 -1 yspan+1 yspan+1],'r','linewidth',2)
    
end

axis ij
% title(['Cell # =  ' num2str(G_C(Celln)) ])

