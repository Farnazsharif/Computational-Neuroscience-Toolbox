load('FM05_1.mat')
[thetaI,thetaf,nbcount]=thetaIndex(acg.acg10t,acg.acg10(:,:));
% I need to save the output
%%
[B,I]=sort(thetaI);

figure
acgN=matnorm(acg.acg10,1);
imagesc(acgN(:,I)')
