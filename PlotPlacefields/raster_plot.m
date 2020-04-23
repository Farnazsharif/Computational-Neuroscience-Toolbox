 
function raster_plot(filename,filename2,Celln)
% raster_plot('DJ57_S2_4T_known',16)
%%
% figure
% filename='DJ57_S2_1T';
% Celln=33;

load([filename '.mat'])
load([filename2 '.mat'])
xttsc=xttscT; 
nPFbin=100;

% xttsc=xttsc(find(xttsc(:,3)==tr_Select(1) & xttsc(:,1)==1):find(xttsc(:,3)==tr_Select(2) & xttsc(:,1)==100),:);%find initial and last trials
matC=reshape(xttsc(:,Celln+4),nPFbin,length(xttsc)/nPFbin)';
matCN=matnorm(matC,1);
[a,b]=size(matC);
% imagesc(matC)
matCN(matCN==0)=NaN;
hold on
for i=a:-1:1
plot(matCN(i,:)-i,'.k')
hold on
end

xlim([0 100])
set(gca,'YTickLabel',[])
ylabel('tr # ','fontsize',11)

