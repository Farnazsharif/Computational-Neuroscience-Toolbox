function [Com]=CFC_linear(filename,Amp,fase,tr,Fs)

% Colomn={LFP ndx,LFP-filtered Amplitude,LFP theta Phase,Position X,speed V,trial t }
% tr
% Fs=1250;
load(filename)
load([filename '_sessInfo.mat'])
%% For comodulation calculation (only has to be calculated once)
nbin = 18;
position =zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin
    position(j) = -pi+(j-1)*winsize;
end

for ghj=1:size(tr,1)
    ghj
    
    tr(ghj,1)
    tr(ghj,end)
    M=[];
    nd1=find(behav.TXVt(:,4)==tr(ghj,1));
    nd2=find(behav.TXVt(:,4)==tr(ghj,end));
    M=behav.TXVt(nd1(1):nd2(end),:);
    [restT,runT]=findrunrestT_Antoni(M(:,1),M(:,3),0.03,1);
    
    NDX=[];NDX_sync=[];
    NDX=round(runT*Fs);
    a1=sessInfo.Epochs.MazeEpoch(1)*Fs+(Fs/2);
    NDX_sync=NDX-a1;
    
    PhaseFreqVector=[];
    AmpFreqVector=[];
    for i=1:size(NDX_sync,1)
        AmpFreqVector=[AmpFreqVector,Amp(:,NDX_sync(i,1):NDX_sync(i,2))];
        PhaseFreqVector=[PhaseFreqVector,fase(:,NDX_sync(i,1):NDX_sync(i,2))];
    end
    
    
    %%
    Comodulogram=single(zeros(size(PhaseFreqVector,1),size(AmpFreqVector,1)));
    
    MI=[];
    Comodulogram=[];
    counter1=0;
    for ii=1:size(PhaseFreqVector,1)
        counter1=counter1+1;
        
        counter2=0;
        for jj=1:size(AmpFreqVector,1)
            counter2=counter2+1;
            
            [MI,MeanAmp]=ModIndex_v2(PhaseFreqVector(ii, :), AmpFreqVector(jj, :), position);
            Comodulogram(counter1,counter2)=MI;
        end
        
    end
    
    
    Com{1,ghj}=Comodulogram;
end





