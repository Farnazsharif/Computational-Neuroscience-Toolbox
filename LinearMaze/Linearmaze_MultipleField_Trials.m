function [trial]=Linearmaze_MultipleField_Trials(SpikeMatrix,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,speed_thr,delta,tr,b)

min_slope=-10; max_slope=10;
nPFbin=max(Rate_Matrix(:,1));

%% ############################################################################################################################
Linear_All=[]; regression_Info_All=[]; rows_All=[]; PF_All=[]; Phase_info_All=[]; SPK_speed_All=[];

for CellN=1:size(SpikeMatrix.num_PF,1)

    clear var SPK_speed regression_Info Phase_info 
    clear var Linear PF rows 
   
    if SpikeMatrix.num_PF(CellN)~=0
        
        for Subfeild=1:SpikeMatrix.num_PF(CellN)
            
            Spiketxvtlph=[];
            Spiketxvtlph=SpikeMatrix.SubfieldSpikes{CellN, 1}{Subfeild, 1};
            Spiketxvtlph(isnan(Spiketxvtlph(:,2))==1,:)=[]; % Remove Nan positions
            Spiketxvtlph(find(Spiketxvtlph(:,3)<speed_thr),:)=[];% Remove speed threshold
            Spiketxvtlph(Spiketxvtlph(:,4)==(0),:)=[]; % Remove zero trials
            
            trials_port=tr;
            nd_int=[];
            for i=1:length(trials_port)
                nd_int=[nd_int;find(Spiketxvtlph(:,4)==trials_port(i))];
            end
            Spiketxvtlph=Spiketxvtlph(nd_int,:);
            
            SPK_speed{1,Subfeild}=[nanmean(Spiketxvtlph(:,3))];
            
            Phase_info_1=[];Phase_info_2=[];Phase_info_3=[];Phase_info_4=[];
            
            Field_edge(1,1)=SpikeMatrix.PhaseBoundiers{CellN, 1}{Subfeild, 1}(1,1);
            Field_edge(1,2)=SpikeMatrix.PhaseBoundiers{CellN, 1}{Subfeild, 2}(1,1);
            
            X1=Field_edge(1,1);
            X2=Field_edge(1,1)+(Field_edge(1,2)-Field_edge(1,1))*0.15; % the first 15% of the place fields
            nd=find(Spiketxvtlph(:,2)>= X1 & Spiketxvtlph(:,2)<= X2);
            
            O=Spiketxvtlph(nd,6);
            [meanR,meanO]=meanphase(O);
            if isempty(O)==0
                Pval = circ_rtest(O);
            else
                Pval=NaN;
            end
            Phase_info_fisrt=[meanO,meanR,Pval];
            
            
            X1=Field_edge(1,1)+(Field_edge(1,2)-Field_edge(1,1))*0.85; % the last 15% of the place fields
            X2=Field_edge(1,2);
            nd=find(Spiketxvtlph(:,2)>= X1 & Spiketxvtlph(:,2)<= X2);
            
            O=Spiketxvtlph(nd,6);
            [meanR,meanO]=meanphase(O);
            if isempty(O)==0
                Pval = circ_rtest(O);
            else
                Pval=NaN;
            end
            Phase_info_last=[meanO,meanR,Pval];
            
            
            X1=Field_edge(1,1); % All of the place fields
            X2=Field_edge(1,2);
            nd=find(Spiketxvtlph(:,2)>= X1 & Spiketxvtlph(:,2)<= X2);
            
            O=Spiketxvtlph(nd,6);
            [meanR,meanO]=meanphase(O);
            if isempty(O)==0
                Pval = circ_rtest(O);
            else
                Pval=NaN;
            end
            Phase_info_Whole=[meanO,meanR,Pval];
            Phase_info{1,Subfeild}(1,:)=Phase_info_fisrt;
            Phase_info{1,Subfeild}(2,:)=Phase_info_last;
            Phase_info{1,Subfeild}(3,:)=Phase_info_Whole;
            
            if size(Spiketxvtlph,1)>2
            [cl_t(1,1),cl_t(1,2),cl_t(1,3),cl_t(1,4),cl_t(1,5)] = cl_corr(Spiketxvtlph(:,2),radtodeg(Spiketxvtlph(:,6)),min_slope, max_slope);
           
            regression_Info{1,Subfeild}=cl_t;
            else
            regression_Info{1,Subfeild}=nan(1,5);
            end
            
            
            Boundreis          =    [Field_edge(1,1) Field_edge(1,2)];
            Linear_Infield{1,Subfeild} =    linear_maze_PFinfo(CellN,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,tr,Boundreis,Spiketxvtlph,b);   
        end
        
        
        clear var Linear rows PF Boundreis Spiketxvtlph
        
        Boundreis     =  [0 1];
        Spiketxvtlph  =  SpikeMatrix.ALLfieldSpikes{CellN};
        [Linear]      =  linear_maze_PFinfo(CellN,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,tr,Boundreis,Spiketxvtlph,b);
        [~,rows,PF]   =  Feild_width_linearMaze(CellN,smooth_rate,delta,Rate_Matrix,tr,b);
    else
        
       SPK_speed              =   nan; 
       Linear                 =   nan(1,16); 
       rows                   =   nan(1,max(Rate_Matrix(:,1)));
       PF                     =   nan(length(tr),nPFbin);
       Phase_info{1,1}        =   nan(3,3);
       Linear_Infield{1,1}    =   nan(1,16);
       regression_Info        =   nan(1,5);
       
    end
    
        SPK_speed_All{CellN,1}             =    SPK_speed;
        Linear_All{CellN,1}                =    Linear;
        rows_All{CellN,1}                  =    rows;
        PF_All{CellN,1}                    =    PF;
        Phase_info_All{CellN,1}            =    Phase_info;
        regression_Info_All{CellN,1}       =    regression_Info;
        Linear_Infield_All{CellN,1}        =    Linear_Infield;
        
end

trial.SPK_speed=SPK_speed_All;
trial.Linear=Linear_All;
trial.rows=rows_All;
trial.PF=PF_All;
trial.Phase_info=Phase_info_All;
trial.regression_Info=regression_Info_All;
trial.Linear_Infield= Linear_Infield_All;
end


