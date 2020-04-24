function [Linear] = linear_maze_PFinfo(CellN,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,tr,Boundreis,Spiketxvtlph,b)
%%
%Linear={SPI per spike,SPI per second,Sparsity,Coefficient,Selectivity,Olpher1,Olpher2,FR,Peak_FR,FanoFactor,StabilityR,StabilityP,MeanVectorlength_r,MeanVectorlength_Pval,Width,animal speed}

%%
trials_port=tr;
nd_int=[];

for i=1:length(trials_port)
    nd_int=[nd_int;find(Rate_Matrix(:,4)==trials_port(i))];
end

nPFbin=max(Rate_Matrix(:,1));
MRate=Rate_Matrix(nd_int,:);
Mphase=Phase_Matrix(nd_int,:);

if Boundreis(1,1)<0.005    
       BoundreisBin=[1 round(Boundreis(1,2).*nPFbin)];
elseif Boundreis(1,2)>1
       BoundreisBin=[round(Boundreis(1,1).*nPFbin) nPFbin];
else   
       BoundreisBin=round(Boundreis.*nPFbin);
end

smoothT=smooth1D(MRate(:,2),smooth_rate,1);
smoothC=smooth1D(MRate(:,CellN+b),smooth_rate,1);
rate=smoothC./smoothT;

Mphase=Mphase(:,CellN+b); % calculate phase


    matC=[];
    matC=reshape(rate,nPFbin,length(rate)/nPFbin)';
    snSpk=reshape(smoothC,nPFbin,length(smoothC)/nPFbin)';
    sTime=reshape(smoothT,nPFbin,length(smoothT)/nPFbin)';

    Rate_Map=matC(:,BoundreisBin(1,1):BoundreisBin(1,2));
    snSpikes=snSpk(:,BoundreisBin(1,1):BoundreisBin(1,2));
    sTimeSpent=sTime(:,BoundreisBin(1,1):BoundreisBin(1,2));

% PF Spatial info #########################################################

    [Information_1,Information_2,Sparsity,Coefficient,Selectivity] = PlaceCellInfo_Antonio(Rate_Map, snSpikes, sTimeSpent);

% PF Firingrate info ######################################################

    FR=mean(mean(Rate_Map));
    Peak_FR=max(mean(Rate_Map));
    Width=[Boundreis(1,2)-Boundreis(1,1)];

% PF Corelation Factor ####################################################

    HalfMatrix=floor(size(Rate_Map,1)/2);
    PF1=mean(Rate_Map(1:HalfMatrix,:),1);
    PF2=mean(Rate_Map(HalfMatrix+1:end,:),1);
    [StabilityR,StabilityP] = corrcoef(PF1,PF2);

% PF Fano Factor ##########################################################

RM=[];
for n = 1:size( Rate_Map,1)
    [~,RM(n,1)] = max( Rate_Map(n,:));
end
FanoFactor= var(RM)/ mean(RM);

% PF Tunity ###############################################################

if isempty(Spiketxvtlph)==1
    MeanVectorlength_r=nan;
    MeanVectorlength_Pval=nan;
else
    theta_PF(:,1)=Spiketxvtlph(:,2)*(2*pi);
    MeanVectorlength_r=[meanphase(theta_PF(:,1))];
    MeanVectorlength_Pval=[circ_rtest(theta_PF(:,1))];
end

% PF Phase info ###########################################################

    Phase_Map=[];
    Phase_Map=reshape(Mphase,nPFbin,length(Mphase)/nPFbin)';
    Phase_Map=Phase_Map(:,BoundreisBin(1,1):BoundreisBin(1,2));
    [track_info_ph,pos_info_val_ph] = bz_olypherInfo_Farn(Phase_Map,2,smooth_phase,'phase');

% #########################################################################
    
    animal_speed=nanmean(MRate(:,3));

    Linear(1,1)=  Information_1;
    Linear(1,2)=  Information_2;
    Linear(1,3)=  Sparsity;
    Linear(1,4)=  Coefficient;
    Linear(1,5)=  Selectivity;
    Linear(1,6)=  max(track_info_ph);
    Linear(1,7)=  max(max(pos_info_val_ph));
    Linear(1,8)=  FR;
    Linear(1,9)=  Peak_FR;
    Linear(1,10)= FanoFactor;
    Linear(1,11)= StabilityR(1,2);
    Linear(1,12)= StabilityP(1,2);
    Linear(1,13)= MeanVectorlength_r;
    Linear(1,14)= MeanVectorlength_Pval;
    Linear(1,15)= Width;
    Linear(1,16)= animal_speed;
    

end