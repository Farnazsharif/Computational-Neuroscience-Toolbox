

function [Cellinfo]=bz_CellMetricsSimple_chronicFarnaz(filename,Celvector)
% filtWaveforms=[];
% Celvector=1:length(G);

load([filename '.mat'])
% Cellinfo=rmfield(Cellinfo,'ACG');
smplrate=spkinfo.samplerate;

for CellN=Celvector
    
    for ChN = 1:size(spkinfo.waveform.average,1)
        if isnan(sum(spkinfo.waveform.average(ChN,:,CellN)))
            peak1(ChN)=0;
            peak2(ChN)=0;
            trough(ChN)=0;
            peak1time(ChN)=0;
            peak2time(ChN)=0;
            troughtime(ChN)=0;
            
        else
            [peak1(ChN),peak2(ChN),trough(ChN),peak1time(ChN),peak2time(ChN),troughtime(ChN)]=spikeinfo(spkinfo.waveform.average(ChN,:,CellN),smplrate);
        end
    end
    
    [crap,Maxindex]=max([peak1-trough]);
    filtWaveforms(CellN,:)=(spkinfo.waveform.average(Maxindex,:,CellN));
end


%%

wave=[]; t_before=[]; t_after=[]; peakA=[]; peakB=[]; trough=[];
for m = 1:size(filtWaveforms,1)
    if isnan(sum(filtWaveforms(m,:)))==1
        t_after(m) = nan;
        peakA(m) = nan;
        peakB(m) = nan;
        trough(m) = nan;
    else
        wave = interp1([1:size(filtWaveforms,2)],zscore(filtWaveforms(m,:)),[1:0.5:size(filtWaveforms,2),size(filtWaveforms,2)],'spline');
        midPoint = round((length(wave)/2));
        [MIN2,I2] = min(wave(midPoint-10:midPoint+10));
        [MAX3,I3] = max(wave(1:midPoint));
        [MAX4,I4] = max(wave((I2+midPoint-11):end));
        t_before(m) = (I2+midPoint-11)-I3;
        t_after(m) = I4;
        peakA(m) = MAX3;
        peakB(m) = MAX4;
        trough(m) = MIN2;
    end
    clear MIN2 MAX3 MAX4 I2 I3 I4
end

Cellinfo.PeaktoTrough = (t_before/(smplrate*2))';
Cellinfo.TroughtoPeak = (t_after/(smplrate*2))';
Cellinfo.AB_ratio = ((peakB-peakA)./(peakA+peakB))';
Cellinfo.trough = (trough)';

%%

for i=1:length(Celvector)
    Cellg=G(Celvector(i));
    nd=find(spk.g==Cellg);
    spiketime{i,1}=spk.i(nd)./(spkinfo.samplerate);
end

for i = 1:length(spiketime)
  
    if length(spiketime{i})>1
    [ccg,time] = CCG(spiketime{i},ones(length(spiketime{i}),1),'binSize',0.001,'duration',0.100); %100ms wide CCG with 0.5ms bins
    ACG_mat=1000*ccg;
    ACG_mat(51) = 0;
    [fmodel,~,~,paut] = fitpyrint(ACG_mat',0:50,0,20);
    Cellinfo.ACG(i,1) .fmodel= fmodel;
    Cellinfo.ACG(i,1) .ydata= ACG_mat;
    Cellinfo.ACG(i,1).xdata = linspace(-0.050,0.05, length(ACG_mat));
    Cellinfo.ACG(i,:).doubleExponentialACG =  paut;
    else
    Cellinfo.ACG(i,1) .fmodel= nan(1,50); 
    Cellinfo.ACG(i,1) .ydata= nan(101,1);
    Cellinfo.ACG(i,1).xdata = nan(1,101);
    Cellinfo.ACG(i,:).doubleExponentialACG = nan(1,5);
    end
    
    
    
    spkTmp=spiketime{i};
    if length(spkTmp)>5
        for jj = 2 : length(spkTmp) - 1
            bursty(jj) =  any(diff(spkTmp(jj-1 : jj + 1)) < 0.006);
        end
        Cellinfo.burstIndexKenji(i) = length(find(bursty > 0))/length(bursty);
    else
        Cellinfo.burstIndexKenji(i)=nan;
    end
    
    
    clear s ccg ACG_mat spkTmp bursty
    
end

%% FR
for i = 1:length(spiketime)
    Cellinfo.FR(i,1) = length(spiketime{i})/(spiketime{i}(end)-spiketime{i}(1));
end

%% Preliminary pyr/int classification
% avg FR > 6Hz or narrow waveform (trought to peak < 5 ms) = INT
for i = 1:length(spiketime)
    if Cellinfo.FR(i,1) > 5 || Cellinfo.TroughtoPeak(i) < 5.5e-4
        Cellinfo.putativeClass(i,1) = 2; % int
    else
        Cellinfo.putativeClass(i,1) = 1; % pyr
    end
end

