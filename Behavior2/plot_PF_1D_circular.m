function [matC]=plot_PF_1D_circular(Matrix,smooth,nPFbin,cellN)
% Matrix=xtvts;
% smooth=10;
% cellN=122;
% nPFbin=50;
%smooth
%%
smoothT=smooth1D(Matrix(:,2),smooth,1);
smoothC=smooth1D(Matrix(:,cellN+4),smooth,1);
rate=smoothC./smoothT;
matC=reshape(rate,nPFbin,length(rate)/nPFbin)';
% plot PF
matC=matnorm(matC,2);
imagesc(matC)
