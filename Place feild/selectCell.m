
filename='FM05_1';
shankID=1;
CellN=101;

chorder=[1 8 2 7 3 6 4 5];
nsamples=40;
nchannels = length(chorder);
fp = fopen([filename '.spk.' num2str(shankID)], 'r');
spkW = fread(fp, [nchannels, inf], 'short');
nspike = size(spkW, 2)/nsamples;
spkW = reshape(spkW, [nchannels, nsamples, nspike]);
spkW = spkW(chorder,:,:);

ndex=find(spk.g==CellN);
Cellspk=spkW(:,:,ndex);

% CCellspk=spkW(:,:,find(spk.g==CellN));

%% claculate SPKts
% spki=Loadres(['FM05_1' '.res.' num2str(16)]);
% [spki,isort]=sort(spki);
% 
% tstamps=TDTtimestampconv(['FM05_1' '_shk1TimeStamp.txt']);
% tstamps=tstamps-tstamps(1);
% 
% nbSample = 2048;
% smplrate = nbSample/mean(diff(tstamps));
% indexTS=(0:length(tstamps)-1)*nbSample + 1;
% indexTS=[indexTS indexTS(end)+nbSample];
% tstamps=[tstamps;tstamps(end)+nbSample/smplrate];
% spk.ts=interp1(indexTS,tstamps,spki);
% 
% [r ,~]=find(abs(spk.ts-0.0028) < 0.0001)
