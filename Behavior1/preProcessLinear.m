%% Pre processing linear track rat
clearvars;clc;
dirData = 'A:\OptoMECLEC\OML';
    animal  = {'19','19','19','19','19','19'};
    sessions  = {'day7','day8','day9','day10'};


%% concatenate files and run kilosort 
for ses = 1:length(sessions)
    dirses = [dirData animal{ses} '\'  sessions{ses}];
    cd(dirses);
    
            bz_ConcatenateDats(pwd,0,1)
            basename = bz_BasenameFromBasepath(pwd);
            if ~exist([basename '.lfp'],'file')
               ResampleBinary([basename '.dat'],[basename '.lfp'],64,1,24); % change for LFPfromDAT
            end
            
            KiloSortWrapper;
            disp(['session' num2str(ses) ' DONE']);
end

%% Get positions 

% read digital input channel
[digitalIn] = bz_getDigitalIn;

% 1- copy csv to indiv subfolder
% 2-write script to automat run:
        % manual (go to subfolder...)
        basename = bz_BasenameFromBasepath(pwd);
        bz_processConvertOptitrack2Behav(basename,'syncSampFq',30000);
%In day folder:
        [digitalIn] = bz_getDigitalIn;
        [tracking] = catPositionsDay(pwd);
        
%% Get spikes after manual clustering
spikes = bz_LoadPhy('noPrompts',true);


%% 
analogin = LoadBinary('analogin.dat','nChannels',1,'freqeucny',30000);

