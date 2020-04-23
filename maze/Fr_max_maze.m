
function [fmap,TimeSpent, Fr_max]=Fr_max_maze(filename,celln,Speed_threshold,smooth)
%%  plot_maze('DJ57_S2_2Oc',10,1000,0,20)
% filename='DJ55_S2_2Oc';
% 
% Speed_threshold=2;
% Speed_noise=70;
% celln=2;
% smooth=10;
% % close all

%
b=4;
load([filename '.mat'])

Speedndx=find(TXYVC(:,4) <= Speed_threshold);
TXYVC(Speedndx,:)=[];

[fmap,TimeSpent]=PF2D(TXYVC(:,1:3),TXYVC(:,celln+b),100,smooth);
 Fr_max=max(max(fmap))
