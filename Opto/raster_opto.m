filename='DJ52_S3_1T';

load([filename '.mat'])


%%acute
%eventTime=opto.laserON_OFF(:,1);
% TimeWindow=[-4 5];

%chronic
eventTime=laserT;
TimeWindow=[-2 3];
binsize=0.1;




for Celln=1:length(G);
    
spikeTime=selectgroup(spk.ts,spk.g,G(Celln));
Tperievent(spikeTime,eventTime,TimeWindow,binsize);

 title(['Cell# =' num2str(G(Celln)) '(' num2str(Celln) ')'  ] , 'FontSize', 15)
cd laser_raster

print([filename '_' num2str(G(Celln))], '-depsc' )
cd ..
close all

% print(['Laser_raster_' num2str(G(Celln))], '-depsc' )
end
%%
