
function [cellposition]=Position(filename,CellGroup,probeID,channel_order_tip,Rf_ch_shank_pow)
%  CellGroup=CellN(:,2);
%  channel_order_tip =[10 9 8 7 6 5 4 3 2 1];
% filename='SFA4_S3_TRD2';

%%
load([filename '.mat']);
% load('Ripple.mat')
Origin=125;
[nshank,shank_spacing,shank_contour,nsite,site_position,site_contour]=SHKlayout(probeID);

for i=1:length(CellGroup)
shank=fix(Cellinfo.CellXYZg(CellGroup(i),4)/100);
CellY=Cellinfo.CellXYZg(CellGroup(i),2);
Re_sh=Rf_ch_shank_pow(shank);
cellposition(i,1)=CellY-site_position{1, 1}(channel_order_tip(Re_sh),2);
cellposition(i,2)=Cellinfo.CellXYZg(CellGroup(i),4);
end