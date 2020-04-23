%% 1- making new xttsc for rest v<5
xttsc=makexttsc(behav.TXDTS,spk,100,5);
save('FM05_1PlaceField2S.mat','xttsc');
%% 2- making new spatial caverege
[SpatialCrest]=spatialcoverage2('FM05_1','FM05_1PlaceField2S');
save('SpatialCrest.mat','SpatialCrest');
%% Burst index & refractory period
[burstIndex_rest,refractoryT_rest]=burst_refT('FM05_1',acg.acg1rest);
save('burstIndex_rest.mat','burstIndex_rest');
save('refractoryT_rest.mat','refractoryT_rest');
[burstIndex_run,refractoryT_run]=burst_refT('FM05_1',acg.acg1run);
save('burstIndex_run.mat','burstIndex_run');
save('refractoryT_run.mat','refractoryT_run');
%% Theta index
[thetaI_rest,thetaf_rest]=thetaIndex(acg.acg10t,acg.acg10rest(:,:));
save('thetaI_rest.mat','thetaIrest');
save('thetaf_rest.mat','thetafrest');
%%
[thetaI_run,thetaf_run]=thetaIndex(acg.acg10t,acg.acg10run(:,:));
save('thetaI_run.mat','thetaI_run');
save('thetaf_run.mat','thetaf_run');
%%
