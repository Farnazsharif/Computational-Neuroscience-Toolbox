
function [spikes] = getCellFeatures(varargin)
%% This function computes the features of the spikes. It does this for the entire length of the session and then also does it for a specified time period (baseline) which is pulled from the maxTime inputted in getWave
%spike info should be pulled from data.mat NOT spikes.cellinfo!!

basepath = pwd;
try [sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', false);
    fs = sessionInfo.rates.wideband;
    nChannels = sessionInfo.nChannels;
catch
    warning('SessionInfo file not found.');
end
fs = sessionInfo.rates.wideband;
spikesFile=dir('*spikes.cellinfo.mat*');
load(spikesFile.name,'spikes');
wfWin = 0.008; % Larger size of waveform windows for filterning
wfWin = round((wfWin * fs)/2);

% spike measures
disp('Computing spike features... ');
for ii = 1 : size(spikes.times,2)
    % firing rate
    spikes.firing_rate(ii) = size(spikes.times{ii},1)/(spikes.times{ii}(end) - spikes.times{ii}(1));
            
    % peak (neg) to peak (pos) duration, as Senzai et al 2017
    [~,tmp] = max(spikes.rawWaveform{ii}(size(spikes.rawWaveform{ii},1)/2:end));
    spikes.spk_duration(ii) = (tmp)/fs; % peak (negative) to peak (second positive) duration
    tmp = tmp + size(spikes.rawWaveform{ii},1)/2 - 1;
            
    % half width
    interpFac = 50;
    mean_spike = interp1(linspace(-wfWin,wfWin,length(spikes.rawWaveform{ii})), ...
        spikes.rawWaveform{ii} - mean(spikes.rawWaveform{ii}),...
        linspace(-wfWin,wfWin,interpFac * length(spikes.rawWaveform{ii})));
    mean_spike = mean_spike - min(mean_spike); mean_spike = abs(mean_spike - max(mean_spike)/2);
    [~,cutpoint_1] = min(mean_spike(1:length(mean_spike)/2));
    [~,cutpoint_2] = min(mean_spike(length(mean_spike)/2:end));
    cutpoint_2 = cutpoint_2 + length(mean_spike)/2;
    spikes.half_width(ii) = ((cutpoint_2 - cutpoint_1)/ interpFac)/fs;
            
    % asymmetry
    spkTemp = spikes.rawWaveform{ii};
    if size(spkTemp,1) < size(spkTemp,2)
        spkTemp = spkTemp';
    end
    [pks, locs] = findpeaks(spkTemp);
    if isempty(locs(locs < size(spikes.rawWaveform{ii},1)/2)) 
        [~,peak1Loc] = min(abs(spkTemp(1:size(spkTemp,1)/2))); 
    else
        peak1Loc = locs(find(max(locs(locs < size(spkTemp,1)/2))==locs));
    end
    peak1 = spkTemp(peak1Loc);
            
    if isempty(pks(locs > size(spkTemp,1)/2))
        peak2Loc = size(spkTemp,1);
    else
        peak2Loc = locs(find(max(pks(locs > size(spkTemp,1)/2)) == pks));
    end
    peak2 = spkTemp(peak2Loc);
                
    spikes.asymmetry(ii) = (peak2 - peak1)/ (peak1 + peak2);
            
    % AUTOCORRELOGRAM FEATURES & DOUBLE EXPONENTIAL FITTING MODEL
    ACG_mat = 1000*CCG(spikes.times{ii},ones(length(spikes.times{ii}),1),'binSize',0.001,'duration',0.1)/length(spikes.times{ii});
    [fmodel,~,~,paut] = fitpyrint(ACG_mat',0:50,0,20);
    spikes.ACG.fmodel{ii} = fmodel;
    spikes.ACG.ydata{ii} = ACG_mat;
    spikes.ACG.xdata{ii} = linspace(-0.050,0.05, length(ACG_mat));
    spikes.doubleExponentialACG(ii,:) =  paut;
            
    % Burstiness (As in Mizuseki et al, 2011). Fraction of spikes
    % with a ISI for following or preceding spikes < 0.006

    bursty = [];
    for jj = 2 : length(spikes.times{ii}) - 1
        bursty(jj) =  any(diff(spikes.times{ii}(jj-1 : jj + 1)) < 0.006);
    end 
    spikes.burstIndex(ii) = length(find(bursty > 0))/length(bursty);
end

if isfield(spikes,'maxTime')
    disp('Computing features for maxTime spikes... ');
    for ii = 1 : size(spikes.times,2)
        spkTmp = spikes.times{ii}(spikes.times{ii}<=spikes.maxTime.value);
        if length(spkTmp) > 5
            % firing rate
            spikes.maxTime.firing_rate(ii) = size(spkTmp,1)/(spkTmp(end) - spkTmp(1));

            % peak (neg) to peak (pos) duration, as Senzai et al 2017
            [~,tmp] = max(spikes.maxTime.rawWaveform{ii}(size(spikes.maxTime.rawWaveform{ii},1)/2:end));
            spikes.maxTime.spk_duration(ii) = (tmp)/fs; % peak (negative) to peak (second positive) duration
            tmp = tmp + size(spikes.rawWaveform{ii},1)/2 - 1;

            % half width
            interpFac = 50;
            mean_spike = interp1(linspace(-wfWin,wfWin,length(spikes.maxTime.rawWaveform{ii})), ...
                spikes.maxTime.rawWaveform{ii} - mean(spikes.maxTime.rawWaveform{ii}),...
                linspace(-wfWin,wfWin,interpFac * length(spikes.maxTime.rawWaveform{ii})));
            [~,cutpoint_1] = min(abs(mean_spike(1 : wfWin * interpFac) - (mean_spike(wfWin * interpFac)/2)));
            [~,cutpoint_2] = min(abs(mean_spike(wfWin * interpFac + 1 : end) - (mean_spike(wfWin * interpFac)/2)));
            cutpoint_2 = cutpoint_2 + wfWin * interpFac;
            spikes.maxTime.half_width(ii) = ((cutpoint_2 - cutpoint_1)/ interpFac)/fs;


            % asymmetry
            spkTemp = spikes.rawWaveform{ii};
            [pks, locs] = findpeaks(spkTemp);
            if isempty(locs(locs < size(spikes.maxTime.rawWaveform{ii},1)/2)) 
                [~,peak1Loc] = min(abs(spkTemp(1:size(spkTemp,1)/2))); 
            else
                peak1Loc = locs(find(max(locs(locs < size(spkTemp,1)/2))==locs));
            end
            peak1 = spkTemp(peak1Loc);

            if isempty(pks(locs > size(spkTemp,1)/2))
                peak2Loc = size(spkTemp,1);
            else
                peak2Loc = locs(find(max(pks(locs > size(spkTemp,1)/2)) == pks));
            end
            peak2 = spkTemp(peak2Loc);

            spikes.maxTime.asymmetry(ii) = (peak2 - peak1)/ (peak1 + peak2);

        % AUTOCORRELOGRAM FEATURES & DOUBLE EXPONENTIAL FITTING MODEL
        ACG_mat = 1000*(CrossCorr(spkTmp,spkTmp,.001,100)/length(spkTmp));
        ACG_mat(51) = 0;
        [fmodel,~,~,paut] = fitpyrint(ACG_mat',0:50,0,20);
        spikes.maxTime.ACG.fmodel{ii} = fmodel;
        spikes.maxTime.ACG.ydata{ii} = ACG_mat;
        spikes.maxTime.ACG.xdata{ii} = linspace(-0.050,0.05, length(ACG_mat));
        try spikes.maxTime.doubleExponentialACG(ii,:) =  paut;
        catch
            keyboard;
        end

        % Burstiness (As in Mizuseki et al, 2011). Fraction of spikes
        % with a ISI for following or preceding spikes < 0.006

        bursty = [];
        for jj = 2 : length(spkTmp) - 1
            bursty(jj) =  any(diff(spkTmp(jj-1 : jj + 1)) < 0.006);
        end 
        spikes.maxTime.burstIndex(ii) = length(find(bursty > 0))/length(bursty);
    
        else
        spikes.maxTime.firing_rate(ii) = NaN;
        spikes.maxTime.spk_duration(ii) = NaN;
        spikes.maxTime.half_width(ii) = NaN;
        spikes.maxTime.burstIndex(ii) = NaN;
        spikes.maxTime.asymmetry(ii) = NaN;
        spikes.maxTime.doubleExponentialACG(ii,:) = NaN(1,5);
        spikes.maxTime.ACG.fmodel{ii}= NaN;
        spikes.maxTime.ACG.ydata{ii}=NaN;
        spikes.maxTime.ACG.xdata{ii}=NaN;
    end
    end
end

% answ = input('Do you want to overwrite data.mat file with spikes data?: ','s');
answ = 'y';
if strcmpi(answ,'y')
    save('dataLocal.mat','spikes');
end


end
