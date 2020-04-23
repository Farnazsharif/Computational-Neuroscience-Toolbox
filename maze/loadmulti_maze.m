%run loadbehav (_RHD) before

function loadmulti_maze(samplerate,shankN,chN,noise_threshold)

%loadmulti_maze(30000,6,10,5);
% filename='DJ56_S2_3O';
% samplerate=30000;
% shankN=6;
% chN=10;

name='uncombined';
dirinfo1=dir(name);
[ds,~]=size(dirinfo1); %%check the number of folders


fln1=0;
for hh=1:ds
    if isempty(strfind(dirinfo1(hh).name,'SFA'))==1 %% detect folder name 'DJ'
        fln1=fln1+1;
    end
end
cd(name)

for j=1:ds-fln1
    

    if isempty(strfind(dirinfo1(j+fln1).name,'m'))==0
        
        cd(dirinfo1(j+fln1).name);
        filename=dirinfo1(j+fln1).name
        

  
        TextfilenameTTL=dir('Hardware*.txt');
        TextfilenameTracking=dir('Track*.txt');
        
        %1) making TXYDV
        
        loadmaze(TextfilenameTTL(1).name,TextfilenameTracking(1).name,filename,samplerate);
        
        if exist('Empty.mat')==0
            
        load([filename '.mat']);
        
        'step 1 done'
        
        % 2) making samplerate
        spkinfo.samplerate = samplerate;
        smplrate=spkinfo.samplerate;
        
        'step 2 done'
        
        
        %3) load spike time, group and waveform
        
        [G,spk,spkinfo.waveform]=SpikeTimeWave3(filename,1:shankN,chN);
        
        spk.ts = spk.i/spkinfo.samplerate;
        
        'step 3 done'
        
        
        %4) spike duration/asymmetry
        [spkinfo.duration,spkinfo.asymmetry]=WaveformFeature(spkinfo.waveform,smplrate);
        
        'step 4 done'
        
        
        %5) spk autocorrelogram
        [restT,runT]=findrunrestT_maze(TXYDV(:,1),TXYDV(:,5),2,1);% speed and duration
        
        [acg.acg10,acg.acg10t]=ACG(spk.i,spk.g,G,smplrate,10,700);
        [acg.acg1,acg.acg1t]=ACG(spk.i,spk.g,G,smplrate,1,50);
        
        [tt,gg]=selectt(spk.i,spk.g,restT,smplrate);
        acg.acg10rest=ACG(tt,gg,G,smplrate,10,700);
        acg.acg1rest=ACG(tt,gg,G,smplrate,1,50);
        
        [tt,gg]=selectt(spk.i,spk.g,runT,smplrate);
        acg.acg10run=ACG(tt,gg,G,smplrate,10,700);
        acg.acg1run=ACG(tt,gg,G,smplrate,1,50);
        
        'step 5 done'
        
        
        %6) refractoryT, burst index
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
        
        'step 6 done'
        
        
        %7) maximum firing rate
        [count,tbin]=binevents(spk.i,spk.g,1,smplrate);
        Crest=selectwindow(count,tbin,restT);
        Crun=selectwindow(count,tbin,runT);
        
        m=[];
        for ii=1:60:length(Crest(:,1))-60
            m=[m;mean(Crest(ii:ii+59,:),1)];
        end
        
        [a,~]=size(m);
        if a==1
            Frate.rest=m;
        else
            Frate.rest=max(m);
            Frate.restmean=mean(m,1);
        end
        
        m=[];
        for ii=1:60:length(Crun(:,1))-60
            m=[m;mean(Crun(ii:ii+59,:),1)];
        end
        
        [a,~]=size(m);
        if a==1
            Frate.run=m;
        else
            Frate.run=max(m);
            Frate.runmean=mean(m,1);
        end
        'step 7 done'
        
        
        %8) make behave
        behav.restT=restT;
        behav.runT=runT;
        
        'step 8 done'
        
        % 9) make XYC bin
        %%
       
        T1=[];
        X1=[];
        Y1=[];
        
        % if isnan(TXYDV(end,1))==1;
        %     TXYDV(end,:)=[];
        % end
        
        NaNndx_t=find(isnan(TXYDV(:,1)));
        TXYDV(NaNndx_t,:)=[];
        
        NaNndx_x=find(isnan(TXYDV(:,2)));
        TXYDV(NaNndx_x,:)=[];
        
        NaNndx_y=find(isnan(TXYDV(:,3)));
        TXYDV(NaNndx_y,:)=[];
        
        NaNndx_D=find(isnan(TXYDV(:,4)));
        TXYDV(NaNndx_D,:)=[];
        

        % removing noise of track
        
        for i=1:10000
            ndx=find(abs(diff(TXYDV(:,4)))>mean(abs(diff(TXYDV(:,4))))*noise_threshold);
            if isempty(ndx)==0
                TXYDV(ndx(1):ndx(1)+1,:)=[];
            else
                break;
            end
        end
        
      
        max_x=max(abs(TXYDV(:,2)));
        max_y=max(abs(TXYDV(:,3)));
        if max_x >= max_y
            max_arena = max_x ;
        else
            max_arena = max_y ;
        end
        
        
        T1=(TXYDV(1,1):0.01:TXYDV(end,1))';
        X1=(interp1(TXYDV(:,1),TXYDV(:,2),T1)+max_arena)/(max_arena*2);
        Y1=(interp1(TXYDV(:,1),TXYDV(:,3),T1)+max_arena)/(max_arena*2);
        V1=interp1(TXYDV(:,1),TXYDV(:,5),T1);
%       D1=interp1(TXYDV(:,1),TXYDV(:,4),T1);
        
        TXYVC=maketxdtsc_DJ([T1,X1,Y1,V1],spk.ts,spk.g,G_C);
        
        'step 9 done'
        
     
     
        % 10) save into .mat file
    
        save([filename '.mat'],'-append','spk','G','spkinfo','behav','acg','Frate','TXYVC')
        
        'all done'
    
        end
        cd ..
    end
end
cd ..
%}














