% clc
% clear

apparatus_order={'SFA5_S2_TRD2','SFA5_S2_maze','SFA5_S2_TRD1'};
% plot_chronic(115,3,0,apparatus_order,Mouse_filename)
% plot_chr(Celln,smooth,Speed_threshold_maze,apparatus_order,Mouse_filename)
%%
clc
clear
apparatus_order={'SFA3_S6_TRD2','SFA3_S6_maze','SFA3_S6_TRD1'};
Mouse_filename='SFA3_S6';
 mkdir([Mouse_filename '_Pics'])
load([Mouse_filename '.mat'])
%%
for i=123:length(G_C)
    i
    
    plot_chronic(i,3,0,apparatus_order,Mouse_filename)
%      cd ([filename '_Pics'])
%     print('-djpeg',['Cell#' num2str(i)])
%     cd ..
    close all
end
%%
