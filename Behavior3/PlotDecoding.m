
%% Plot decoding
Q={'mse_rate','mse_rate_chance','mse_phase','mse_phase_chance','mse_phase_cos','mse_phase_cos_chance', ...
    'mse_phase_sin','mse_phase_sin_chance','mse_phase_all','mse_phase_all_chance'};
U_C=[];
for k=1:5
   
U_C=[U_C;(table2array(positionDecodingMaxCorr_C{1,k}.results(:,1:7)))];
end

U_P=[];
for k=1:5
U_P=[U_P;(table2array(positionDecodingMaxCorr_P{1,k}.results(:,1:7)))];
end

%%
close all
Colomn=2;
CA1.C_1=U_C(:,1);%SPT_C_all(:,Colomn);%SPT_C_all(:,Colomn);%(SPT_C_all-SPT_P_all)./(SPT_C_all+SPT_P_all);
CA1.P_1=U_P(:,1);%SPT_P_all(:,Colomn);%[U3;U3];%U3;
CA1.C_2=U_C(:,2);%[U2(:,Colomn); U5(:,Colomn)];
CA1.P_2=U_P(:,2);%[U3(:,Colomn);U6(:,Colomn)];%[U6;U6];%;

CA3.C_1=U_C(:,3);%SPT_C_all(:,Colomn);%SPT_C_all(:,Colomn);%(SPT_C_all-SPT_P_all)./(SPT_C_all+SPT_P_all);
CA3.P_1=U_P(:,3);%SPT_P_all(:,Colomn);%[U3;U3];%U3;
CA3.C_2=U_C(:,7);%[U2(:,Colomn); U5(:,Colomn)];
CA3.P_2=U_P(:,7);%[U3(:,Colomn);U6(:,Colomn)];%[U6;U6];%;



Boxplot_f(CA3,CA1,Q{Colomn})
barplot_f(CA3,CA1,Q{Colomn})
% ylim([0 300])

