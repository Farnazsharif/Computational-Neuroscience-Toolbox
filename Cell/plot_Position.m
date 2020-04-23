

filename='SFA4_S3_TRD1';
filename2='SFA4-S3-TRD2';
CellVector=TRD1_single;%Rew2;

% ndx=[find(fix(G_C(CellVector')/100)==2)]; 
% CellVector(ndx)=[];
load([filename 'PlaceField_Total_G_C.mat'])

probeID=3;
load([filename '.mat']);
CellN=[];
[sortndx] = Sort_Cell(filename,CellVector,xttscT);
CellN(:,1)=[CellVector(sortndx)];
for i=1:length(CellVector)
    CellN(i,2)=find(G==G_C(CellN(i,1)));
end
CA=1;
channel_order_tip=[10 9 8 7 6 5 4 3 2 1];

%%
 Cell_vector=TRD1_single_cue;
 [C shank_Avg_C]=Shank_Avg(filename, Cell_vector,probeID,channel_order_tip,Rf_ch_shank_pow(:,6),6);
 Cell_vector=TRD1_single_P;
 [P shank_Avg_P]=Shank_Avg(filename, Cell_vector,probeID,channel_order_tip,Rf_ch_shank_pow(:,6),6);
 Cell_group.C=C;
 Cell_group.P=P;
 Cell_group.shank_Avg_C=shank_Avg_C;
Cell_group.shank_Avg_P=shank_Avg_P;
save([filename '.mat'],'-append','Cell_group')
% figure
% hist((shank_Avg_C-shank_Avg_p))
%%
figure('Position', [100 100 1000 500] ) 
Clayout_2(Cellinfo.PrePostG_T,filename,CellN([22],2),Rf_ch_shank_pow,CA,probeID ,channel_order_tip,{ 'r*','bo','ko'});
ylim([0 250])

 %%
 figure
 cellposition=Position(filename,CellN(:,2),probeID,channel_order_tip,Rf_ch_shank_pow(:,6));
 plot (cellposition,'--*r')

%%
plot(ones(1, length(C_1(:,5))),C_1(:,5),'*')

 %%
 figure
 V1= cellposition([1:15 24]);
 V2= cellposition([16:23]);
 V3= cellposition;
 
 subplot(1,2,1)

 plot(sort(V1),(1:length(V1))/length(V1),'Linewidth',3,'color','g')
 hold on
 plot(sort(V2),(1:length(V2))/length(V2),'Linewidth',3,'color','k')
 hold on
 plot(sort(V3),(1:length(V3))/length(V3),'Linewidth',3,'color','r')
  title([filename2])
  
subplot(1,2,2)
bar(2,mean(V1),'FaceColor',[99/256 209/256 62/256],'EdgeColor','r','LineWidth',2)
set(gca,'XTickLabel',[])
hold on
b2=errorbar(2,mean(V1),std(V1)/sqrt(length(V1)));
set(b2,'linewidth',2,'color','k')
hold on
bar(3,mean(V2),'FaceColor',[88/256 88/256 88/256],'EdgeColor','r','LineWidth',2)
set(gca,'XTickLabel',[])
hold on
b3=errorbar(3,mean(V2),std(V2)/sqrt(length(V2)));
set(b3,'linewidth',2,'color','k')

hold on
bar(1,mean(V3),'FaceColor','r','EdgeColor','r','LineWidth',2)
set(gca,'XTickLabel',[])
hold on
b1=errorbar(1,mean(V3),std(V3)/sqrt(length(V3)));
set(b1,'linewidth',2,'color','k')

clc
C={'k','r'};
[Pa,Pk] = Stest( V1,V2 )
if Pa < 0.1 || Pk < 0.1
    i=2;
else
    i=1;
end
title(['Pa=' num2str(round(Pa,4)) '  '  'Pk=' num2str(round(Pk,4))],'color',C{i})