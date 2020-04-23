function allcellplot (filename)
% filename='FM05_1'
%% first cancate G1 and G3
% G3=cat(1,PrePostG3,PrePostG1);
% PrePostG4=unique(G3,'rows'); 
%%


load(filename);

CellID(:,1)=fix((G/100));
CellID(:,2)=mod(G,100);
CellID=[CellID G];
for i=1:length(CellID)

cellplot(CellID(i,1),CellID(i,2),CellID(i,3))
% print('-dtiff',['shank' num2str(shankID) 'cell' num2str(CellN) ])
% print('-depsc',['shank' num2str(CellID(i,1)) 'cell' num2str(CellID(i,2)) ])
close all
end
end
