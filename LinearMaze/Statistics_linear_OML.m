%% Generate Trials  #############################################################################################################

clc
clear

dirData = 'Y:\OptoMECLEC\';
animal  = {'OML18','OML18','OML18','OML19','OML19'};
sessions  = {'day1','day2','day4','day2','day3'};


for ses=1:5;
    clc
    dirses = [dirData animal{ses} '\'  sessions{ses}];
    cd(dirses);
    filename=[animal{ses} sessions{ses}]
    load([filename '.mat'])
    
    %% #############################################################################################################
    %Initial conditions
    
    trN=5;speed_thr=0.04;
    delta=0.1;smooth_rate=1;smooth_phase=1;b=5;
    Rate_Matrix=map.xtvtls2;
    Phase_Matrix=map.xtvtlph;
    BinSize=(max(behav.txvtl(:,2))./max(map.xtvtls2(:,1)));
    
    TR=unique(behav.txvtl(:,4));
    TR(TR==0)=[];
    r0=[1 trN];
    TR_odd=find(mod(TR,2)==1);
    TR_even=find(mod(TR,2)==0);
    
    %% segregate trials ##############################################################################################################
    
    TR_all=TR(TR_odd);
    Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
    TRi=TR_all(1:Seg*trN);
    TRi=reshape(TRi,[trN Seg])';
    tr=TRi;
    
    tiral=[];
    tic
    SpikeMatrix=spikes.CellMetrics.cellinfo_Infield_off;
    for j=1:(Seg)
        %     [tiral{1,j}]=Linearmaze_stats_OML(filename,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,BinSize,speed_thr,delta,tr(j,:),b);
        [tiral{1,j}]=Linearmaze_MultipleField_Trials(SpikeMatrix,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,speed_thr,delta,tr(j,:),b);
    end
    [tiral{1,(Seg+1)}]=Linearmaze_MultipleField_Trials(SpikeMatrix,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,speed_thr,delta,TR_all,b);
    spikes.CellMetrics.cellinfo_Infield_off.trials=tiral;
    toc
    
    TR_all=TR(TR_even);
    Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
    TRi=TR_all(1:Seg*trN);
    TRi=reshape(TRi,[trN Seg])';
    tr=TRi;
    
    tiral=[];
    tic
    SpikeMatrix=spikes.CellMetrics.cellinfo_Infield_on;
    for j=1:(Seg)
        [tiral{1,j}]=Linearmaze_MultipleField_Trials(SpikeMatrix,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,speed_thr,delta,tr(j,:),b);
    end
    toc
    [tiral{1,(Seg+1)}]=Linearmaze_MultipleField_Trials(SpikeMatrix,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,speed_thr,delta,TR_all,b);
    spikes.CellMetrics.cellinfo_Infield_on.trials=tiral;
    toc
    save([filename '.mat'],'-append','spikes')
    
end


%% Merge Trials ###########################################################################################################################################
clc
clear
dirData = 'Y:\OptoMECLEC\';
animal  = {'OML18','OML18','OML18','OML19','OML19'};
sessions  = {'day1','day2','day4','day2','day3'};

Trial_all=[];
for ses=1:5%length(sessions)
    ses
    dirses = [dirData animal{ses} '\' sessions{ses}];
    cd(dirses);
    filename=[animal{ses} sessions{ses}]
    load([filename '.mat'])

        cd Metrics
        tiral=[];
        load('tiral_1.mat')
        Trial_all{1,ses}=tiral;
        tiral=[];
        load('tiral_2.mat')
        Trial_all{2,ses}=tiral;
        load('PhasePrec.mat')
        statsPP{1,ses}=statsPPoff;
        statsPP{2,ses}=statsPPon;
        cd ..        
        spikes_all{ses}=spikes;
        cd Metrics
        tiral=[];
        load('tiral_1.mat')
        Trial_all{2,ses}=tiral;
        tiral=[];
        load('tiral_2.mat')
        Trial_all{1,ses}=tiral;
        load('PhasePrec.mat')
        statsPP{2,ses}=statsPPoff;
        statsPP{1,ses}=statsPPon;
        cd ..
        

end



%% Plot for all Sesions
cellinfo_Infield=[];
   
spikes=[];
Mode=1;
Cell_type=[];Cell_Position=[];num_PF=[];PhaseBoundiers=[];regression_Info=[];SubfieldSpikes=[];ALLfieldSpikes=[];AllTrials=[];
sesNumber=[1:5];
for ses=sesNumber
    ses
    Cell_type=       [Cell_type;       spikes_all{1, ses}.CellMetrics.cellinfo_Infield_off.Cell_type];
    Cell_Position=   [Cell_Position;   spikes_all{1, ses}.CellMetrics.cellinfo_Infield_off.Cell_Position];
    num_PF=          [num_PF;          spikes_all{1, ses}.CellMetrics.cellinfo_Infield_off.num_PF];
    PhaseBoundiers=  [PhaseBoundiers;  spikes_all{1, ses}.CellMetrics.cellinfo_Infield_off.PhaseBoundiers];
    regression_Info= [regression_Info; spikes_all{1, ses}.CellMetrics.cellinfo_Infield_off.regression_Info];
    SubfieldSpikes=  [SubfieldSpikes;  spikes_all{1, ses}.CellMetrics.cellinfo_Infield_off.SubfieldSpikes];
    ALLfieldSpikes=  [ALLfieldSpikes;  spikes_all{1, ses}.CellMetrics.cellinfo_Infield_off.ALLfieldSpikes];
    AllTrials{ses,1}=                  spikes_all{1, ses}.CellMetrics.cellinfo_Infield_off.trials;
end
cellinfo_Infield.Cell_type{1,Mode}=Cell_type;
cellinfo_Infield.Cell_Position{1,Mode}=Cell_Position;
cellinfo_Infield.num_PF{1,Mode}=num_PF;
cellinfo_Infield.PhaseBoundiers{1,Mode}=PhaseBoundiers;
cellinfo_Infield.regression_Info{1,Mode}=regression_Info;
cellinfo_Infield.SubfieldSpikes{1,Mode}=SubfieldSpikes;
cellinfo_Infield.ALLfieldSpikes{1,Mode}=ALLfieldSpikes;
cellinfo_Infield.Alltrials{1,Mode}=AllTrials;

Mode=2;
Cell_type=[];Cell_Position=[];num_PF=[];PhaseBoundiers=[];regression_Info=[];SubfieldSpikes=[];ALLfieldSpikes=[];AllTrials=[];

for ses=sesNumber
    ses
    Cell_type=       [Cell_type;       spikes_all{1, ses}.CellMetrics.cellinfo_Infield_on.Cell_type];
    Cell_Position=   [Cell_Position;   spikes_all{1, ses}.CellMetrics.cellinfo_Infield_on.Cell_Position];
    num_PF=          [num_PF;          spikes_all{1, ses}.CellMetrics.cellinfo_Infield_on.num_PF];
    PhaseBoundiers=  [PhaseBoundiers;  spikes_all{1, ses}.CellMetrics.cellinfo_Infield_on.PhaseBoundiers];
    regression_Info= [regression_Info; spikes_all{1, ses}.CellMetrics.cellinfo_Infield_on.regression_Info];
    SubfieldSpikes=  [SubfieldSpikes;  spikes_all{1, ses}.CellMetrics.cellinfo_Infield_on.SubfieldSpikes];
    ALLfieldSpikes=  [ALLfieldSpikes;  spikes_all{1, ses}.CellMetrics.cellinfo_Infield_on.ALLfieldSpikes];
    AllTrials{ses,1}=                  spikes_all{1, ses}.CellMetrics.cellinfo_Infield_off.trials;
end

cellinfo_Infield.Cell_type{1,Mode}=Cell_type;
cellinfo_Infield.Cell_Position{1,Mode}=Cell_Position;
cellinfo_Infield.num_PF{1,Mode}=num_PF;
cellinfo_Infield.PhaseBoundiers{1,Mode}=PhaseBoundiers;
cellinfo_Infield.regression_Info{1,Mode}=regression_Info;
cellinfo_Infield.SubfieldSpikes{1,Mode}=SubfieldSpikes;
cellinfo_Infield.ALLfieldSpikes{1,Mode}=ALLfieldSpikes;
cellinfo_Infield.Alltrials{1,Mode}=AllTrials;

spikes.CellMetrics.cellinfo_Infield=cellinfo_Infield;





SPK_speed_All=[];Linear_All=[];rows_All=[];Phase_info_All=[];
for Mode=1:2

    M=[];SPK_speed=[];Linear=[];rows=[];Phase_info=[];
    
    M=spikes.CellMetrics.cellinfo_Infield.Alltrials{1,Mode};  
    
    for ses=sesNumber
        
        seg=size(M{ses},2);
        SPK_speed=   [SPK_speed; M{ses, 1}{1, seg}.SPK_speed];
        Linear=      [Linear; M{ses, 1}{1, seg}.Linear];
        rows=        [rows; M{ses, 1}{1, seg}.rows];
        Phase_info=  [Phase_info; M{ses, 1}{1, seg}.Phase_info];
        
    end
    
    SPK_speed_All{1,Mode}=SPK_speed;
    Linear_All{1,Mode}=Linear;
    rows_All{1,Mode}=rows;
    Phase_info_All{1,Mode}=Phase_info;
end


%% Do statistics

bar_Matrix=[];
M=spikes.CellMetrics.cellinfo_Infield;  
for Mode=1:2 
    clear var Pyr_CA1 int_CA1 Pyr_DG int_DG Deep_Pyr Sup_Pyr Deep_int Sup_int ... 
    Pyr_singlePF Pyr_MultiplePF Pyr_Silence Deep_Pyr_SinglePF Sup_Pyr_MultiplePF
    Pyr_CA1=find(M.Cell_type{1,Mode}==1 & M.Cell_Position{1,Mode}<3);
    int_CA1=find(M.Cell_type{1,Mode}==2 & M.Cell_Position{1,Mode}<3);
    
    Pyr_DG=find(M.Cell_type{1,Mode}==1 & M.Cell_Position{1,Mode}==3);
    int_DG=find(M.Cell_type{1,Mode}==2 & M.Cell_Position{1,Mode}==3);
    
    Deep_Pyr=find(M.Cell_type{1,Mode}==1 & M.Cell_Position{1,Mode}==1);
    Sup_Pyr=find(M.Cell_type{1,Mode}==1 & M.Cell_Position{1,Mode}==2);
    Deep_int=find(M.Cell_type{1,Mode}==2 & M.Cell_Position{1,Mode}==1);
    Sup_int=find(M.Cell_type{1,Mode}==2 & M.Cell_Position{1,Mode}==2);
    
    Pyr_singlePF=find(M.Cell_type{1,Mode}==1 & M.Cell_Position{1,Mode}<3 & M.num_PF{1,Mode}==1);
    Pyr_MultiplePF=find(M.Cell_type{1,Mode}==1 & M.Cell_Position{1,Mode}<3 & M.num_PF{1,Mode}>1);
    Pyr_Silence=find(M.Cell_type{1,Mode}==1 & M.Cell_Position{1,Mode}<3 & M.num_PF{1,Mode}==0);
    
    Deep_Pyr_SinglePF=find(M.Cell_type{1,Mode}==1 & M.Cell_Position{1,Mode}==1 & M.num_PF{1,Mode}==1);
    Sup_Pyr_SinglePF=find(M.Cell_type{1,Mode}==1 & M.Cell_Position{1,Mode}==2 & M.num_PF{1,Mode}==1);
    
    Deep_Pyr_MultiplePF=find(M.Cell_type{1,Mode}==1 & M.Cell_Position{1,Mode}==1 & M.num_PF{1,Mode}>1);
    Sup_Pyr_MultiplePF=find(M.Cell_type{1,Mode}==1 & M.Cell_Position{1,Mode}==2 & M.num_PF{1,Mode}>1);
    
    bar_Matrix(:,Mode)=[length(Pyr_CA1),length(int_CA1),length(Pyr_DG),length(int_DG),...
        length(Deep_Pyr),length(Sup_Pyr),length(Deep_int),length(Sup_int),...
        length(Pyr_singlePF),length(Pyr_MultiplePF),length(Pyr_Silence),length(Deep_Pyr_SinglePF),...
        length(Sup_Pyr_SinglePF),length(Deep_Pyr_MultiplePF),length(Sup_Pyr_MultiplePF)]
end


figure
b=bar(bar_Matrix);

somenames={'Pyr-CA1'; 'int-CA1';'Pyr-DG'; 'int-DG'; 'Deep-Pyr';'Sup-Pyr';'Deep-int';'Sup-int';'Pyr-SinglePF';'Pyr-MultiplePF';'Pyr-Silence';'Deep-Pyr-SinglePF' ;'Sup-Pyr-SinglePF' ;'Deep-Pyr-MultiplePF' ;'Sup-Pyr-MultiplePF' };
set(gca,'xticklabel',somenames,'fontsize',15)
xtickangle(45)



%%
close all
figure('position',([50 50 500 600]))
a1P=0:0.001:0.9; a2P=a1P+0.1; x_win=[a1P;a2P]';
thr=1;
% somenames_title={'Pyr-singlePF';'Deep-Pyr-SinglePF' ;'Sup-Pyr-SinglePF' };
somenames_title={'Pyr-MultiplePF-CA1';'Deep-Pyr-MultiplePF' ;'Sup-Pyr-MultiplePF'};
% somenames_title={'DG-singlePF' };
% somenames_title={'DG-MultiplePF';'DG-MultiplePF'};

SP=[1 3 5;2 4 6];
num_SubPF=1;
meanph=[];
meanph_All=[];
M=spikes.CellMetrics.cellinfo_Infield;  
for j=1:3;
    
    for Mode=1:2;
%         Cellvector=[];
%         
                if j==1
                    Pyr_SinglePF_CA1=find(M.Cell_type{1,Mode}==1 & M.num_PF{1,Mode}==1 & M.Cell_Position{1,Mode}<3);
                    Cellvector=Pyr_SinglePF_CA1;
                elseif j==2
                    Pyr_SinglePF_CA1_Deep=find(M.Cell_type{1,Mode}==1 & M.num_PF{1,Mode}==1 & M.Cell_Position{1,Mode}==1);
                    Cellvector=Pyr_SinglePF_CA1_Deep;
                elseif j==3
                    Pyr_SinglePF_CA1_Sup=find(M.Cell_type{1,Mode}==1 & M.num_PF{1,Mode}==1 & M.Cell_Position{1,Mode}==2);
                    Cellvector=Pyr_SinglePF_CA1_Sup;
                end
                CL={'m'}
        
%         if j==1
%             Pyr_MultiplePF_CA1=find(M.Cell_type{1,Mode}==1 & M.num_PF{1,Mode}>1 & M.Cell_Position{1,Mode}<3);
%             Cellvector=Pyr_MultiplePF_CA1;
%         elseif j==2
%             Pyr_MultiplePF_CA1_Deep=find(M.Cell_type{1,Mode}==1 & M.num_PF{1,Mode}>1 & M.Cell_Position{1,Mode}==1);
%             Cellvector=Pyr_MultiplePF_CA1_Deep;
%         elseif j==3
%             Pyr_MultiplePF_CA1_Sup=find(M.Cell_type{1,Mode}==1 & M.num_PF{1,Mode}>1 & M.Cell_Position{1,Mode}==2);
%             Cellvector=Pyr_MultiplePF_CA1_Sup;
%         end
%         CL={'r','b','g'}
        
        
        %         if j==1
        %             Pyr_MultiplePF_CA1=find(M.Cell_type==1 & M.num_PF==1 & M.Cell_Position==3);
        %             Cellvector=Pyr_MultiplePF_CA1;
        %         elseif j==2
        %             Pyr_MultiplePF_CA1_Deep=find(M.Cell_type==1 & M.num_PF>1 & M.Cell_Position==3);
        %             Cellvector=Pyr_MultiplePF_CA1_Deep;
        %         end
        %         CL={'r','b','g'}
        
        
        if isempty(Cellvector)==1
            continue
        end
        meanRO=[];
        PhMatrix=[];
        for i=1:length(Cellvector)
            
            CellN=Cellvector(i);
            if iscell(M.SubfieldSpikes{1,Mode}{CellN, 1})==1
                
                if num_SubPF<= size(M.SubfieldSpikes{1,Mode}{CellN, 1},1)
                    for subfield=num_SubPF%:size(M.SubfieldSpikes{CellN, 1},1)
                        
                        T=[];
                        T=M.SubfieldSpikes{1,Mode}{CellN, 1}{subfield,1};
                        X=T(:,2);
                        X=(X-min(X))/(max(X)-min(X));
                        
                        Ph=T(:,6);
                        PHbin_Subfield=phasepersesion_count( X,x_win,Ph,thr);
                        PhMatrix{subfield,1}(i,:) = PHbin_Subfield(:,2);
                        
                        O=[];
                        O=Ph;
                        O(isnan(O)==1)=[];
                        O((O)==0)=[];
                        [meanRO(i,1),meanRO(i,2)]=meanphase(O);
                        
                    end
                end
                
                
            end
            
        end
        meanph_All{j, Mode}=meanRO;

        
        
        
        subplot(3,2,SP(Mode,j))
        for CellN=1:size(PhMatrix{subfield, 1},1)
            plot(radtodeg(PhMatrix{subfield, 1}(CellN,:)),'.','color',[190/256 190/256 190/256])
            hold on
        end
        
        hold on
        
        meanRO=[];
        for i=1:size(x_win,1)
            O=PhMatrix{subfield, 1}(:,i);
            O(isnan(O)==1)=[];
            O((O)==0)=[];
            [meanRO(i,1),meanRO(i,2)]=meanphase(O);
        end
        
        meanph{j, Mode}=meanRO;
        
        plot(1:size(x_win,1),radtodeg(meanRO(:,2)),'.','color',CL{subfield},'linewidth',2)
        ylim([0 360])
        xlim([0 size(x_win,1)])
        set(gca,'xtick',[])
        title(['subfield',num2str(subfield)])
        
        
        if Mode==2
            xlabel([somenames_title{j}  '  Laser-ON'])
        else
            xlabel([somenames_title{j} '  Laser-Off'])
        end
        
    end
    
end



%% PULL ALL Feilds
clc
% close all
figure('position',([50 50 500 600]))
a1P=0:0.001:0.9; a2P=a1P+0.1; x_win=[a1P;a2P]';
thr=1;
% somenames_title={'Pyr-singlePF';'Deep-Pyr-SinglePF' ;'Sup-Pyr-SinglePF' };
somenames_title={'Pyr-ALLsubPF-CA1';'Deep-Pyr-ALLsubPF' ;'Sup-Pyr-ALLsubPF'};
% somenames_title={'DG-singlePF' };
% somenames_title={'DG-MultiplePF';'DG-MultiplePF'};

SP=[1 3 5;2 4 6];
num_SubPF=2;
meanph=[];
PhMatrix=[];
M=spikes.CellMetrics.cellinfo_Infield;  
for j=1:3;
    PhMatrix=[];
    
    for Mode=1:2;

        Cellvector=[];
        if j==1
            Pyr_MultiplePF_CA1=find(M.Cell_type{1,Mode}==1  & M.Cell_Position{1,Mode}<3 & M.num_PF{1,Mode}>0  );
            Cellvector=Pyr_MultiplePF_CA1;
             size(Cellvector)
        elseif j==2
            Pyr_MultiplePF_CA1_Deep=find(M.Cell_type{1,Mode}==1  & M.Cell_Position{1,Mode}==1 & M.num_PF{1,Mode}>0  );
            Cellvector=Pyr_MultiplePF_CA1_Deep;
             size(Cellvector)
        elseif j==3
            Pyr_MultiplePF_CA1_Sup=find(M.Cell_type{1,Mode}==1  & M.Cell_Position{1,Mode}==2 & M.num_PF{1,Mode}>0 );
            Cellvector=Pyr_MultiplePF_CA1_Sup;
            size(Cellvector)
        end
        CL={'r','b','g'}
        
        
        
        
        meanRO=[];
        PhMatrix=[];
        for i=1:length(Cellvector)
            
            CellN=Cellvector(i);
            if iscell(M.SubfieldSpikes{1,Mode}{CellN, 1})==1
                

                 if M.num_PF{1,Mode}(CellN,1)>0
                    
                        for subfield=1%:size(M.SubfieldSpikes{1,Mode}{CellN, 1},1)
                        T=[];
                        T=M.SubfieldSpikes{1,Mode}{CellN, 1}{subfield,1};
                        X=T(:,2);
                        X=(X-min(X))/(max(X)-min(X));
                        
                        Ph=T(:,6);
                        PHbin_Subfield=phasepersesion_count( X,x_win,Ph,thr);
                        PhMatrix = cat(1,PhMatrix, PHbin_Subfield(:,2)');
                        end
                   
                  end
                
                
            end
            
        end
        
        size(PhMatrix)
        
        
        
        subplot(3,2,SP(Mode,j))
        for CellN=1:size(PhMatrix,1)
            plot(radtodeg(PhMatrix(CellN,:)),'.','color',[190/256 190/256 190/256])
            hold on
        end
        
        hold on
        
        meanRO=[];
        for i=1:size(x_win,1)
            O=PhMatrix(:,i);
            O(isnan(O)==1)=[];
            O((O)==0)=[];
            [meanRO(i,1),meanRO(i,2)]=meanphase(O);
        end
        
        meanph{j, Mode}=meanRO;
        
        plot(1:size(x_win,1),radtodeg(meanRO(:,2)),'.','color',CL{1},'linewidth',2)
        ylim([0 360])
        xlim([0 size(x_win,1)])
        set(gca,'xtick',[])
        title(['All Fields'])
        
        
        if Mode==2
            xlabel([somenames_title{j}  '  Laser-ON'])
        else
            xlabel([somenames_title{j} '  Laser-Off'])
        end
        
    end
    
end



%%
close all
SP=[1 3 5;2 4 6]';
C={'k','r'};
mean_ph=[];
figure
for KM=1:2
    for j=1:3
        subplot(3,2,SP(j, KM))
        mean_ph=meanph_All{j, KM}(:,2);
        polarhistogram(mean_ph,0:2*pi/13:2*pi,'linewidth',3,'facecolor',[220./256 220/256 201/256],'edgecolor',C{KM},'FaceAlpha',.3)
    end
end

%%
close all
M=meanph_All;
Q={'meanph'};

CA1.C_1=radtodeg(M{1,1}(:,2));
CA1.P_1=radtodeg(M{1,2}(:,2));
CA1.C_2=radtodeg(M{2,1}(:,2));
CA1.P_2=radtodeg(M{2,2}(:,2));
% U2(find(isnan(U2)==1))=[];
% U3(find(isnan(U3)==1))=[];
% U5(find(isnan(U5)==1))=[];
% U6(find(isnan(U6)==1))=[];

CA3.C_1=zeros(1,10)';
CA3.P_1=zeros(1,10)';
CA3.C_2=radtodeg(M{3,1}(:,2));
CA3.P_2=radtodeg(M{3,2}(:,2));


Boxplot_f(CA3,CA1,Q{1})
barplot_f(CA3,CA1,Q{1})

%% Whole field analysis

Q={'1 SPI per spike','2 SPI per second','3 Sparsity','4 Coefficient','5 Selectivity','6 Olpher1','7 Olpher2','8 FR','9 animal speed'};
j=1;
M=spikes.CellMetrics.cellinfo_Infield;
Col=9;
if j==1
    Pyr_MultiplePF_CA1=[];
    Pyr_MultiplePF_CA1=find(M.Cell_type{1, 1}==1 & M.num_PF{1, 1}>1 & M.Cell_Position{1, 1}<3);
    Cellvector_off=Pyr_MultiplePF_CA1;
    
    Pyr_MultiplePF_CA1=[];
    Pyr_MultiplePF_CA1=find(M.Cell_type{1, 2}==1 & M.num_PF{1, 2}>1 & M.Cell_Position{1, 2}<3);
    Cellvector_on=Pyr_MultiplePF_CA1;
    
elseif j==2
    
    Pyr_MultiplePF_CA1_Deep=find(M.Cell_type{1, 1}==1 & M.num_PF{1, 1}>1 & M.Cell_Position{1, 1}==1);
    Cellvector_off=Pyr_MultiplePF_CA1_Deep;
    
    Pyr_MultiplePF_CA1_Deep=find(M.Cell_type{1, 2}==1 & M.num_PF{1, 2}>1 & M.Cell_Position{1, 2}==1);
    Cellvector_on=Pyr_MultiplePF_CA1_Deep;
    
elseif j==3
    
    Pyr_MultiplePF_CA1_Sup=find(M.Cell_type{1, 1}==1 & M.num_PF{1, 1}>1 & M.Cell_Position{1, 1}==2);
    Cellvector_off=Pyr_MultiplePF_CA1_Sup;
    
    Pyr_MultiplePF_CA1_Sup=find(M.Cell_type{1, 2}==1 & M.num_PF{1, 2}>1 & M.Cell_Position{1, 2}==2);
    Cellvector_on=Pyr_MultiplePF_CA1_Sup;
    
end


X_off=[];X_on=[];
X_off=cell2mat(Linear_All{1, 1});
X_off=X_off(Cellvector_off,:);

X_on=cell2mat(Linear_All{1, 2});
X_on=X_on(Cellvector_on,:);


CA1.C_1=X_off(:,Col);
CA1.P_1=X_on(:,Col);
CA1.C_2=zeros(1,10)';
CA1.P_2=zeros(1,10)';
CA3=CA1;



Boxplot_f(CA3,CA1,Q{Col})
barplot_f(CA3,CA1,Q{Col})




Width=M.PhaseBoundiers{1, 1}{22, 1}{1, 2}(1,1)-M.PhaseBoundiers{1, 1}{22, 1}{1, 1}(1,1);












