%% Plot selected
minTime=0.01;
mapFRca1Std = ca1pStdFrM{21}{3}; infoFRca1Std = ca1pStdFrS{21}{3}.info1;
mapPhca1Std = ca1pStdPhM{21}{3}; infoPhca1Std = ca1pStdPhS{21}{3}.oly1;
mapFRca3Std = ca3pStdFrM{34}{4}; infoFRca3Std = ca3pStdFrS{34}{4}.info1;
mapPhca3Std = ca3pStdPhM{34}{4}; infoPhca3Std = ca3pStdPhS{34}{4}.oly1;

mapFRca1Obj = ca1pObjFrM{4}{4}; infoFRca1Obj = ca1pObjFrS{4}{4}.info1;
mapPhca1Obj = ca1pObjPhM{4}{4}; infoPhca1Obj = ca1pObjPhS{4}{4}.oly1;
mapFRca3Obj = ca3pObjFrM{2}{3}; infoFRca3Obj = ca3pObjFrS{2}{3}.info1;
mapPhca3Obj = ca3pObjPhM{2}{3}; infoPhca3Obj = ca3pObjPhS{2}{3}.oly1;

figure;
subplot(2,4,1);PlotColorMap(mapFRca1Std.z,mapFRca1Std.time,'threshold',minTime); 
set(gca,'XTickLabel',[],'YTickLabel',[]);title(['info =' num2str(infoFRca1Std,2)]);
ax1=subplot(2,4,2);PlotColorMap(mapPhca1Std.z,mapPhca1Std.time,'threshold',minTime); 
set(gca,'XTickLabel',[],'YTickLabel',[]);title(['info =' num2str(infoPhca1Std,2)]);clim([0 2*pi]);
subplot(2,4,3);PlotColorMap(mapFRca3Std.z,mapFRca3Std.time,'threshold',minTime); 
set(gca,'XTickLabel',[],'YTickLabel',[]);title(['info =' num2str(infoFRca3Std,2)]);
ax2=subplot(2,4,4);PlotColorMap(mapPhca3Std.z,mapPhca3Std.time,'threshold',minTime); 
set(gca,'XTickLabel',[],'YTickLabel',[]);title(['info =' num2str(infoPhca3Std,2)]);clim([0 2*pi]);

subplot(2,4,5);PlotColorMap(mapFRca1Obj.z,mapFRca1Obj.time,'threshold',minTime); 
set(gca,'XTickLabel',[],'YTickLabel',[]);title(['info =' num2str(infoFRca1Obj,2)]);
ax3=subplot(2,4,6);PlotColorMap(mapPhca1Obj.z,mapPhca1Obj.time,'threshold',minTime); 
set(gca,'XTickLabel',[],'YTickLabel',[]);title(['info =' num2str(infoPhca1Obj,2)]);clim([0 2*pi]);
subplot(2,4,7);PlotColorMap(mapFRca3Obj.z,mapFRca3Obj.time,'threshold',minTime); 
set(gca,'XTickLabel',[],'YTickLabel',[]);title(['info =' num2str(infoFRca3Obj,2)]);
ax4=subplot(2,4,8);PlotColorMap(mapPhca3Obj.z,mapPhca3Obj.time,'threshold',minTime); 
set(gca,'XTickLabel',[],'YTickLabel',[]);title(['info =' num2str(infoPhca3Obj,2)]);clim([0 2*pi]);

colormap(ax1,'hsv');colormap(ax2,'hsv');colormap(ax3,'hsv');colormap(ax4,'hsv');
