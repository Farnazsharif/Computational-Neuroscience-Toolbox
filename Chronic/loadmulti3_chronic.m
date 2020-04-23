function loadmulti3_chronic(filename)
%%
filename='SFA5_S4_TRD1';
load([filename '.mat'])
smplrate=spkinfo.samplerate;

Elnb=6;nbchperEl=10;
[G,spk,spkinfo.waveform] = SpikeTimeWave3_new(filename,1:Elnb,nbchperEl);
spk.ts = spk.i/spkinfo.samplerate;
'step 1 done'

%%
%2) spike duration/asymmetry

[spkinfo.duration,spkinfo.asymmetry]=WaveformFeature(spkinfo.waveform,smplrate);

'step 2 done'

%%
%3) spk autocorrelogram , speed threshold=5
[restT,runT]=findrunrestT(behav.TXDTS(:,1),behav.TXDTS(:,3),5,2,1);


[acg.acg10,acg.acg10t]=ACG(spk.i,spk.g,G,smplrate,10,700);
[acg.acg1,acg.acg1t]=ACG(spk.i,spk.g,G,smplrate,1,50);

[tt,gg]=selectt(spk.i,spk.g,restT,smplrate);
acg.acg10rest=ACG(tt,gg,G,smplrate,10,700);
acg.acg1rest=ACG(tt,gg,G,smplrate,1,50);

[tt,gg]=selectt(spk.i,spk.g,runT,smplrate);
acg.acg10run=ACG(tt,gg,G,smplrate,10,700);
acg.acg1run=ACG(tt,gg,G,smplrate,1,50);

'step 3 done'

%4) refractoryT, burst index
for ii=1:length(G)
    
    if sum(acg.acg1(:,ii))>100
        
        [maxacg,Imaxacg]=max(acg.acg1(41:50,ii));
        
        edge=mean(acg.acg1(1:10,ii));
		
        amp=maxacg-edge;
		
        if amp>0
			burstIndex(ii)=amp/maxacg;
		else
			burstIndex(ii)=amp/edge;
		end
		
		refractoryT(ii)=ACGrefractoryT(acg.acg1(:,ii));
		
	else
	
		burstIndex(ii)=nan;
		refractoryT(ii)=nan;
	end
	
end

acg.burstIndex = burstIndex;
acg.refractoryT = refractoryT;

'step 4 done'

%5) maximum firing rate
%%
[count,tbin]=binevents(spk.i,spk.g,1,smplrate);
Crest=selectwindow(count,tbin,restT);
Crun=selectwindow(count,tbin,runT);

m=[];
for ii=1:60:length(Crest(:,1))-60
    m=[m;mean(Crest(ii:ii+59,:),1)];
end
Frate.rest=max(m);
Frate.restmean=mean(m,1);
m=[];
for ii=1:60:length(Crun(:,1))-60
    m=[m;mean(Crun(ii:ii+59,:),1)];
end
Frate.run=max(m);
Frate.runmean=mean(m,1);
        
'step 5 done'
%%
%6) save into .mat file
behav.restT=restT;
behav.runT=runT;

save([filename '.mat'],'-append','spk','G','spkinfo','behav','acg','Frate')

'all done'


