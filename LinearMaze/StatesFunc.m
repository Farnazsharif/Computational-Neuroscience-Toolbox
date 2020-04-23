function [STATES] = StatesFunc(dir,filesRUN,filesSLEEP,nChans,thetaChan,threshold)
% [STATES] = StatesFunc(dir,filesRUN,filesSLEEP,thetaChan,threshold)
%   dir = directory with all eegfiles
%   filesRUN = .eeg files of running / behaviour
%   filesSLEEP = .eeg files of sleep
%   nChans =  number of channels in files
%   thetaChan = channel to calculate theta/delta power ratio
%   threshold = theta/delta power ratio threshold 
%   Antonio FR, 2015

    lfpRUN = [];
    lfpSLEEP = [];
    for i=1:length(filesRUN)
        temp = readmulti([dir filesRUN{i} '.eeg'],nChans,thetaChan);
        lfpRUN = cat(1,lfpRUN,temp);
    end
    for i=1:length(filesSLEEP)
        temp = readmulti([dir filesSLEEP{i} '.eeg'],nChans,thetaChan);
        lfpSLEEP = cat(1,lfpSLEEP,temp);
    end
    
    Fs = 1250;
    t1 = size(lfpSLEEP,1); % length sleep session
    LFPtheta = cat(1,lfpSLEEP,lfpRUN);

%% Find theta epochs
TDratio_thresh = threshold; 

t = (0:(size(LFPtheta,1)-1))'/Fs;
Fs_spec = Fs/10;
t_spec = downsample(t,10);
[thetaSpec,theta_f] = long_wavespec(downsample(LFPtheta, 10), t_spec, 10, [1 20], 1);
smooth_sec = 2;
theta_Finds = find_inds(5, theta_f):find_inds(8, theta_f);
delta_Finds = find_inds(1, theta_f):find_inds(4, theta_f);
thetaDeltaRatio = runavg([zeros(round(smooth_sec*Fs_spec/2),1); mean(thetaSpec(theta_Finds,:),1)'./mean(thetaSpec(delta_Finds,:),1)'; ...
    zeros(round(smooth_sec*Fs_spec/2),1)], round(smooth_sec*Fs_spec/2)*2+1);
thetaDeltaInterp = [lin_interp(thetaDeltaRatio, 10); ones(9,1)*thetaDeltaRatio(end)];
thetaBools = thetaDeltaInterp > TDratio_thresh;
thetaInds = find(thetaBools);
if length(thetaBools) > length(LFPtheta)
thetaBools=thetaBools(1:length(LFPtheta));
end
NOthetaInds = find(not(thetaBools));

%% plot para comprobar resultado
    params = struct ('tapers',{[3 5]},'pad',{0},'Fs',{1250},'fpass',{[1 20]},'err',{[1 0.95]},'trialave',{1});
    [spect2,fbins2,Serr2] = mtspectrumc(LFPtheta(NOthetaInds,:),params);
    [spect1,fbins1,Serr1] = mtspectrumc(LFPtheta(thetaInds,:),params);
    figure;plot_vector(spect1,fbins1,'n',Serr1,'b',2);hold on;
           plot_vector(spect2,fbins2,'n',Serr2,'r',2);hold on;
    set(gca,'XLim',[params.fpass(1) params.fpass(2)]);
    
%% Create STATES variable 
    tsleep = t1; % sep sleep session / run session
    tsT = find(thetaInds >= tsleep);
    tsT = tsT(1);
    tsNT = find(NOthetaInds >= tsleep);
    tsNT = tsNT(1);

    STATES = struct;
    STATES.REM = thetaInds(1:tsT);
    STATES.REM = STATES.REM(1:end-1); % quito el ultimo punto pq a veces da problemas
    STATES.SWS = NOthetaInds(1:tsNT);
    STATES.SWS = STATES.SWS(1:end-1);
    STATES.RUN = thetaInds(tsT:end);
    STATES.RUN = STATES.RUN(1:end-1);
    STATES.QUIET = NOthetaInds(tsNT:end);
    STATES.QUIET = STATES.QUIET(1:end-1);
end

