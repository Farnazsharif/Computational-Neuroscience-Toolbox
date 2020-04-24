function [matC_right,matC_left]=plot_PF_1D(Matrix,smooth,cellN,b)
% Matrix=xtvts1;
% % smooth=10;
% cellN=135;
% nPFbin=50;
%smooth
% Matrix=xtvtls2;
% left=even=up
% b=5;
%%
nPFbin=max(Matrix(:,1));
smoothT=smooth1D(Matrix(:,2),smooth,1);
smoothC=smooth1D(Matrix(:,cellN+b),smooth,1);
rate=smoothC./smoothT;

Nd_even=find(mod(Matrix(:,4),2)==0);
Nd_odd=find(mod(Matrix(:,4),2)==1);
rate_left=rate(Nd_even,:);
rate_right=rate(Nd_odd,:);
%%
matC_left=reshape(rate_left,nPFbin,length(rate_left)/nPFbin)';
matC_right=reshape(rate_right,nPFbin,length(rate_right)/nPFbin)';

% plot PF
matC=[matC_left;matC_right];
matC=matnorm(matC,2);
imagesc(matC)
hold on
plot([0 size(matC_left,2)],[size(matC_left,1) size(matC_left,1)],'--w','linewidth',2)

%%

% trials_port=tr;
% nd_int=[];
% for i=1:length(trials_port)
%     nd_int=[nd_int;find(Rate_Matrix(:,4)==trials_port(i))];
% end
% 
% nPFbin=max(Rate_Matrix(:,1));
% 
% MRate=Rate_Matrix(nd_int,:);
% Mphase=Phase_Matrix(nd_int,:);
% 
% cellN=Cellvector(i);
% smoothT=smooth1D(MRate(:,2),smooth_rate,1);
% smoothC=smooth1D(MRate(:,cellN+4),smooth_rate,1);
% rate=smoothC./smoothT;
% 







