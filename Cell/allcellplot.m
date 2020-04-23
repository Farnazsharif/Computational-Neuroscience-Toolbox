function allcellplot (filename,Cell_Start,Winrange,CA)
%%
%  R_ndex=find(ismember(Ripple,RippleT_Single)==1);194
filename='SFA4_S4_TRD2';
Cell_Start=14;
Cell_end=486;
CA=1;
load([filename '.mat']);
Winrange=100;
CellID(:,1)=fix((G/100));
CellID(:,2)=mod(G,100);
CellID=[CellID G];
 mkdir('Cue_zone')
 mkdir('Path_zone')
%%
for g=1:16
N(g)=length(find(fix((G/100))==g));
end
%%

Cell_vector=[Cell_group.C(:,3)];
for h= 1:length(Cell_vector);%Cell_Start:Cell_end%length(CellID)
    
Cell_vector(h)                              
  
% % plot phase
figure ('position', [1600 0 1000  1000])

cellplotN2(filename,Cell_vector(h),CellID,1000,64)
close all

% plot place fields 
% figure ('position', [1600 0 1000  1000])
% cellplotN(CellID(h,1),CellID(h,2),CellID(h,3),filename,CA)
% close all

% % plot Ripples 
% % figure
% % Cellripple(filename,Winrange,CellID,h)
% % close all

end

end
%  shankID=9;
%  CellN=3;
%  Cellg=903;
% cellplot_rest(CellID(i,1),CellID(i,2),CellID(i,3))