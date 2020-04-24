%% OMLlinear ICA + CFC

clearvars;
addpath('A:\Antonio\MATLAB\fieldtrip-20150209'); 
ft_defaults;

path = 'A:\OptoMECLEC\OML18\day1\';
cd(path);

SR = 30000;
chOrd = [29 35 25 39 21 43 17 47 16 46 19 45 18 44 20 42]+1;

    xmlfile = ([path 'day1.xml']);
    eegfile = ([path 'day1.lfp']);
    hdr = ft_read_header((xmlfile), 'headerformat', 'neuroscope_xml');

    % gamma filtered lfp
    cfg = []; cfg.datafile = eegfile ; 
    cfg.dataformat = 'neuroscope_bin'; cfg.headerformat = 'neuroscope_bin';
    cfg.hpfilter = 'yes'; cfg.hpfreq = 25;   
    cfg.lpfilter = 'yes'; cfg.lpfreq = 200;     
    lfptemp = ft_preprocessing(cfg);
    lfp = lfptemp;
    lfp.trial{1} = lfptemp.trial{1}(chOrd,:); 
    lfp.label = lfptemp.label(chOrd); clear lfptemp;
    save([dirData sessions{ses} '\lfpG.mat'],'lfp','-v7.3'); 
    clear lfp xmlfile eegfile hdr;

    % ICA
    cfg = [];cfg.method='runica';             
    cfg.numcomponent = 8;
    ICA = ft_componentanalysis(cfg, lfp);  
    save('icaG8.mat' ,'ICA', '-v7.3');
%     clear lfp ICA      
    
    figure;plot(ICA.topo);
    ICs = [8 4 5 1 3 2 7];
    
    % theta epochs
    load('day1.MergePoints.events.mat')
    basename = bz_BasenameFromBasepath(pwd);
    [thetaEpochs,~] = thetaStates(basename,64,47,1);
    thetaTask = ExcludeIntervals(thetaEpochs.thetaInt,(MergePoints.timestamps_samples([1 4],:)/30000)*1250);
    thetTaskInd = intvl_to_inds(thetaTask);
    
    % 60Hz nothch
    cfg = []; cfg.bsfilter = 'yes'; cfg.bsfreq = [58 62];   
    ICAf = ft_preprocessing(cfg,ICA);    
    
    
%% CFC
load('LasreInt.mat');
thetaON = Restrict(thetaTask,LasreInt.On_Intv*1250);
thetaONind = intvl_to_inds(thetaON);
thetaOFF = Restrict(thetaTask,LasreInt.Off_Intv*1250);
thetaOFFind = intvl_to_inds(thetaOFF);
laserONind = round(intvl_to_inds(LasreInt.On_Intv*1250));
laserOFFind = round(intvl_to_inds(LasreInt.Off_Intv*1250));

Fs = 1250;
freqs1 = 5:0.5:10;
freqs2 = 70:5:150;
Nbins=18; binedges=linspace(-pi,pi,Nbins+1); binedges(1)=-pi-0.1; binedges(Nbins+1)=pi+0.1;

    lfpRef = LoadBinary([basename '.lfp'],'frequency',1250,'Channels',47,'nChannels',64);
    
    tortICAoff = CFCtort2(lfpRef(thetaOFFind),ICA.trial{1}(:,thetaOFFind )',freqs1,freqs2,1250,1);
    tortICAon = CFCtort2(lfpRef(thetaONind),ICA.trial{1}(:,thetaONind )',freqs1,freqs2,1250,1);
    colormap jet;

    tortICAoff = CFCtort2(lfpRef(laserOFFind(1:750000)),ICA.trial{1}(:,laserOFFind(1:750000))',freqs1,freqs2,1250,1);colormap jet;    
    tortICAon = CFCtort2(lfpRef(laserONind(1:750000)),ICA.trial{1}(:,laserONind(1:750000))',freqs1,freqs2,1250,1);colormap jet;    
 
    
    