function    [SleepState]=Theta_states(basePath)
% Append theta state info to the sleepscoremaster output
% F.Sharif 2020            
            [~,recordingname] = fileparts(basePath);
            savefolder =basePath;
            bz_sleepstatepath = fullfile(savefolder,[recordingname,'.SleepState.states.mat']);
            load([recordingname '.SleepState.states.mat'])
 
            thratio = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.thratio;
            ththr = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.histsandthreshs.THthresh;
            emg = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.EMG;
            emgthr = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.histsandthreshs.EMGthresh;
            ts = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.t_clus;
            theta_run = thratio>ththr & emg>emgthr;
            Theta_NDX=find(thratio>ththr & emg>emgthr);
            
            A(:,1)=SleepState.idx.states;
            A(Theta_NDX,2)=7;  
            Non_thetaNDX=find(A(:,1)==1& A(:,2)==0);       
            A(Non_thetaNDX,2)=9;  
            
            SleepState.idx.theta_states.states=A(:,2);
            SleepState.idx.theta_states.timestamps=SleepState.idx.timestamps;
            SleepState.idx.theta_states.statenames{7} = 'Theta';
            SleepState.idx.theta_states.statenames{9} = 'Non_theta';
            [INT] = bz_IDXtoINT(SleepState.idx.theta_states);
            SleepState.ints.Thetastate=INT.Thetastate;
            SleepState.ints.Non_thetastate=INT.Non_thetastate;
            save(bz_sleepstatepath,'SleepState');            
end

