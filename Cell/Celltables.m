%% Saving in mat_file

filename='FM05_1';

Cellinfo.thetaI_run=thetaI_run;
Cellinfo.thetaI_rest=thetaIrest;

Cellinfo.thetaf_run=thetaf_run;
Cellinfo.thetaf_rest=thetafrest;

Cellinfo.thetaI=thetaI;
Cellinfo.thetaf=thetaf;

Cellinfo.burstIndex_rest=burstIndex_rest;
Cellinfo.burstIndex_run=burstIndex_run;

Cellinfo.refractoryT_rest=refractoryT_rest;
Cellinfo.refractoryT_run=refractoryT_run;

Cellinfo.CellXYZg=CellXYZg;

Cellinfo.PrePostG4=PrePostG4;
Cellinfo.PrePostG1=PrePostG1;

Cellinfo.SpatialC=SpatialC;
Cellinfo.SpatialCrest=SpatialCrest;
Cellinfo.SpatialC=SpatialC;

Cellinfo.Ndx=Ndx;

Cellinfo.Cellchn=Cellchn;

Cellinfo.Cellphase=zeros(length(G),1);

Cellinfo.spkN=zeros(length(G),1);

Cellinfo.Position=repmat({blanks(1)},[425 1]);
Cellinfo.MI=repmat({blanks(1)},[425 1]);
Cellinfo.AcgSim=repmat({blanks(1)},[425 1]);
Cellinfo.Pattern=repmat({blanks(1)},[425 1]);
Cellinfo.Celltype=repmat({blanks(1)},[425 1]);
Cellinfo.Note=repmat({blanks(1)},[425 1]);

Cellinfo.SpeedM=repmat({blanks(1)},[425 1]);
save(filename, '-append', 'Cellinfo')

%%
LastName = {'Smith';'Johnson';'Williams';'Jones';'Brown'};
Age = [38;43;38;40;49];
Height = [71;69;64;67;64];
Weight = [176;163;131;133;119];
BloodPressure = [124 93; 109 77; 125 83; 117 75; 122 80];
Position={' ', '', '' , '' , ''}';
% A= table(Age,Height,Weight,BloodPressure,'RowNames',LastName)
A= table(Age,Height,Weight,BloodPressure,Position)
%%

T={'cell#' 'FRrun' 'FRrest' 'Bur' 'BurRun' 'BurRest'  'RT' 'RTRun' 'RTRest' 'Thet' 'ThetRest' 'ThetRun' 'Synap' 'MI' 'AcgSim' 'Pattern' 'Position' 'ch#' 'celltype'  'Note'};

%%
% str=reshape(blanks(6),1,6);
Cellinfo.str=repmat({blanks(4)},[425 6]);
for j= 1:length(G)
 Cellg=G(j);
 [a,b]=find(Cellg==Cellinfo.PrePostG4) ;  
if b==1
 Cellinfo.str{j,1}=('Pre');
 
    for ii=1:length(a)
    Cellinfo.str{j,1+ii}=([ 'post' num2str(ii) '=' num2str(Cellinfo.PrePostG4(a(ii),2)) ]);
    end
    
elseif b==2
 Cellinfo.str{j,1}=('Post');
 
    for ii=1:length(a)
    Cellinfo.str{j,1+ii}=([ 'pre' num2str(ii) '=' num2str(Cellinfo.PrePostG4(a(ii),1)) ]);
    end
 
end
end
% str
save('FM05_1', '-append', 'Cellinfo')

%%
T=[];
CellN=(G);

FRun=Frate.run';
FRest=Frate.rest';

BI=acg.burstIndex';
BIRun=Cellinfo.burstIndex_run';
BIRest=Cellinfo.burstIndex_rest';

Ref=acg.refractoryT';
RefRun=Cellinfo.refractoryT_run';
RefRest=Cellinfo.refractoryT_rest';

Thet=(Cellinfo.thetaf')*1000;
ThetRun=(Cellinfo.thetaf_run')*1000;
ThetRest=(Cellinfo.thetaf_rest')*1000;
Ch=Cellinfo.Cellchn(:,2);
Synap=Cellinfo.str(:,1);
MI=repmat({blanks(1)},[425 1]);
AcgSim=repmat({blanks(1)},[425 1]);
Pattern=repmat({blanks(1)},[425 1]);
% Position=repmat({blanks(1)},[425 1]);
Celltype=repmat({blanks(1)},[425 1]);
% Note=repmat({blanks(1)},[425 1]);
T= table(CellN,FRun,FRest,BI,BIRun,BIRest,Ref,RefRun,RefRest,Thet,ThetRun,ThetRest,MI,Ch,Synap,AcgSim,Pattern,Celltype);
%%
close all
CellN=1001
ss=find(G==CellN);
cd MI
load(['MI_' num2str(ss) '.mat'])
plot(MI)
cd ../
%%

phaseprocession(filename,PhaseFreqVector_a,cellN,ProbN,ch)
[a,b]=max(MI(1:10));
MI_precent=b*100
%%

[a,b]=max(MI(1:10))
MI_precent=a*100
h=389;
phaseprocession('FM05_1',b,Cellinfo_Prob_1.Cellchn(h,1),1,Cellinfo_Prob_1.Cellchn(h,2));
%%
close all
clc
CellN=1577;
ss=find(G==CellN);
ii=Cellinfo.Cellchn(ss,2)
%%
 i=1;
 Frq=23;
 ii=Frq;
 Powermap( 'FM05_1','po_', i,ii)
 
set(gcf,'color','w');
 %%
close all
for Frq=[1 5 6 7 8 10 27 44 71 103 250];
%     figure
i=Frq;
Powermap( 'FM05_1','po_', i,ii)
end
%%
clc
CellN=1210;
ss=find(G==CellN);
FR_rest=Frate.rest(ss)
FR_run=Frate.run(ss)
%%
Cellinfo.Position{ss}={'SupPyr'};
Cellinfo.MI{ss}={'NO'};
Cellinfo.AcgSim{ss}={'same'};
Cellinfo.Pattern{ss}={'NO'};
Cellinfo.Celltype{ss}={''};
Cellinfo.Note{ss}={'ThMI='};
Cellinfo.SpeedM{ss}={'Both'};
% Cellinfo.SpeedM{ss}={'Yes'};
save('FM05_1', '-append', 'Cellinfo')
%%
clc
ss=313;
Cellinfo.Position{ss}
Cellinfo.MI{ss}
Cellinfo.AcgSim{ss}
Cellinfo.Pattern{ss}
Cellinfo.Celltype{ss}
Cellinfo.Note{ss}
Cellinfo.SpeedM{ss}










