%% 1 first step is doing Ripplefinder
%% 2 find the Ripple Refrence channels
% clear
% clc
% Mouse_filename='SFA3_S3_TRD2';
Mouse_filename=filename;
% save([Mouse_filename '.mat'],'-append','Beltinfo')
%filename='SFA4_S2_TRD2';
NumberProbs=1;%2;
% load('Rt4.mat')
probeID=3;
Winrange=100;
ShankN=6;%16;
load([Mouse_filename '.mat']);
smplrate=spkinfo.samplerate;
Elnb=6;nbchperEl=10;

%% 4 Refractory time rest and run

for ii=1:length(G) 
    
    if sum(acg.acg1run(:,ii))>100
         
    refractoryTrun(ii)=ACGrefractoryT(acg.acg1run(:,ii));
    
    else
         
	refractoryTrun(ii)=nan;
    end
    
    if sum(acg.acg1rest(:,ii))>100 
         
    refractoryTrest(ii)=ACGrefractoryT(acg.acg1rest(:,ii));
    
    else
         
	refractoryTrest(ii)=nan;
    end
    
end

acg.refractoryTrun=refractoryTrun;
acg.refractoryTrest=refractoryTrest;

save([Mouse_filename '.mat'],'-append','acg')
 
%% 5 Cell_Position

Starting_shk=1;
%  CellXYZg_1=CellPosition(G,spkinfo.waveform.average,[5 4 6 3 7 2 8 1],spkinfo.samplerate,2,Starting_shk);
CellXYZg_1=CellPosition(G,spkinfo.waveform.average,[10 9 8 7 6 5 4 3 2 1],spkinfo.samplerate,probeID,Starting_shk);
Cellinfo.CellXYZg_1=CellXYZg_1;
save([Mouse_filename '.mat'],'-append','Cellinfo')

 % ******* if we have 2 probs
 CellXYZg_2=[];
 if NumberProbs==2
 Starting_shk=9;
 CellXYZg_2=CellPosition(G,spkinfo.waveform.average,[5 4 6 3 7 2 8 1],spkinfo.samplerate,probeID,Starting_shk);
 Cellinfo.CellXYZg_2=CellXYZg_2;
 save([Mouse_filename '.mat'],'-append','Cellinfo')
 end
 
 CellXYZg=([CellXYZg_1; CellXYZg_2]);
 Cellinfo.CellXYZg=CellXYZg;
 save([Mouse_filename '.mat'],'-append','Cellinfo')
 
%% 6 Cell_Channel %% wroonggggggggg!!!! this is for acute
 ProbN=1;
%filename='CellXYZg_1';

 [Cellchn_1]=Cellsite(ProbN,CellXYZg_1,probeID);
 Cellinfo.Cellchn_1=Cellchn_1;
 save([Mouse_filename '.mat'],'-append','Cellinfo')
 
 % ******* if we have 2 probs
 Cellchn_2=[];
 if NumberProbs==2 
 ProbN=2;
 filename='CellXYZg_2';
 [Cellchn_2]=Cellsite(ProbN,CellXYZg_2);
 Cellinfo.Cellchn_2=Cellchn_2;
 save([Mouse_filename '.mat'],'-append','Cellinfo')
 end  
 
 Cellchn_T=([Cellchn_1; Cellchn_2]);
 Cellinfo.Cellchn_T=Cellchn_T;
 save([Mouse_filename '.mat'],'-append','Cellinfo')
 
 
%%
[~,~,spkinfo.waveform] = SpikeTimeWave3_new(Mouse_filename,1:Elnb,nbchperEl);
[spkinfo.duration,spkinfo.asymmetry]=WaveformFeature(spkinfo.waveform,smplrate);
'step 2 done'

% spkinfo = rmfield(spkinfo,'refractoryTrun');
% spkinfo = rmfield(spkinfo,'refractoryTrest');
save([Mouse_filename '.mat'],'-append','spkinfo')

%% 3 PrePostG

PrePostG=DetectSynapticInteraction(spk.i,spk.g,spkinfo.samplerate,5,1);

Cellinfo.PrePostG_T=PrePostG;


%%
open Ripplefinder

%% Add Pre postG

load('Cellinfo.mat')
PrePostG=Cellinfo.PrePostG_T;
Cellinfo.PrePostG_T=PrePostG;
save([Mouse_filename '.mat'],'-append','Cellinfo')

% remove ripples V > 5
% clear
% Mouse_filename='SFA4_S1_TRD1';
% load([Mouse_filename '.mat']);
%% Modify Ripple
load('Ripple.mat')
speedsm=x2v(behav.TXDTS(:,3),behav.TXDTS(:,1),1);
TXDTSV=[behav.TXDTS speedsm];
rippleT=Ripple_Time';
for i=1:length(rippleT)
[~,b(i)]=min(abs(rippleT(i)-behav.TXDTS(:,1)));
end

rippleT=[rippleT TXDTSV(b,6)];
ndx=find(rippleT(:,2)>5);
length(ndx)
%%
rippleT(ndx,:)=[];
Ripple_Time=rippleT(:,1);

save('Ripple','-append','Ripple_Time')
Rt(ndx,:)=[];
save(['Ripple.mat'],'-append','Rt')

 %% 7 Cell MI Overal, rest, run
mkdir('MI_Probs')
Cell_Start=86;       
CellMI(Mouse_filename,Cellinfo.Cellchn_T,Cell_Start,64);

%% 8 plot place fields
xttsc = makexttsc(behav.TXDTS,spk,100,5,'c',Mouse_filename);
save([ Mouse_filename 'PlaceField.mat'],'xttsc')

xttscT = makexttsc(behav.TXDTS,spk,100,0,'c',Mouse_filename);
save([ Mouse_filename 'PlaceField_Total.mat'],'xttscT')

xttscR = makexttsc_rest(behav.TXDTS,spk,100,5,'c',Mouse_filename);
save([ Mouse_filename 'PlaceField_Rset.mat'],'xttscR')

%% 9 Creat folders and subfolders

for i=1:ShankN
 mkdir('Ripples_images',(['Shank' num2str(i)]))
 mkdir('Placefield_images',(['Shank' num2str(i)]))  
 mkdir('Phase_images',(['Shank' num2str(i)]))  
end
 
%% 10 First Plot Phase and place fields
Rt=Rt4;
Cell_Start=3;
allcellplot (Mouse_filename,Cell_Start,Rt,Winrange,ch_Ref_Prob_1,ch_Ref_Prob_2)
%%

% files = dir('FM04-2-new')
dirinfo = dir('FM04')
cd 'FM04'
cd(dirinfo(8).name)




