function [tiral]=Linearmaze_stats_OML(filename,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,BinSize,speed_thr,delta,tr,b)
%%
    
portion=0.5;
% load([filename '_sessInfo.mat'])
load([filename '.mat'])


%############################################################################################################################
    Linear_All=[];Field_edge_All=[];rows_All=[];PF_All=[];Phase_info_All=[];intervf_All=[];In_Spike_All=[];SPK_speed=[];%2

    for j=1:size(spikes.Cellinfo,2)
    CellN=j;
        
    Spiketxvtlph=[];
    Spiketxvtlph=spikes.Cellinfo(CellN).cell_txvtlph;
    
    Spiketxvtlph(isnan(Spiketxvtlph(:,2))==1,:)=[]; % Remove Nan positions
    Spiketxvtlph(find(Spiketxvtlph(:,3)<speed_thr),:)=[];% Remove speed threshold
    Spiketxvtlph(Spiketxvtlph(:,4)==(0),:)=[]; % Remove zero trials
        
        
        trials_port=tr;
        nd_int=[];
        for i=1:length(trials_port)
            nd_int=[nd_int;find(Spiketxvtlph(:,4)==trials_port(i))];
        end
        Spiketxvtlph=Spiketxvtlph(nd_int,:);
        
        
        SPK_speed=[SPK_speed;nanmean(Spiketxvtlph(:,3))];
      
       
        Linear=[];Field_edge=[];rows=[];PF=[];Phase_info=[];intervf=[];In_Spike=[];%1
        
        [Linear] = linear_maze_PFinfo(CellN,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,tr,b);

        [Field_edge,rows,PF] = Feild_width_linearMaze(CellN,smooth_rate,delta,Rate_Matrix,tr,b);
        
        Phase_info_1=[];Phase_info_2=[];Phase_info_3=[];Phase_info_4=[];
        
        X2=Field_edge(1,1)+(Field_edge(1,2)-Field_edge(1,1))*portion;
        X1=Field_edge(1,1);
        interv=[X1 X2];
        [Phase_info_1]=PF_sectionphase(Spiketxvtlph,behav.txvtl,interv,BinSize,'first');
        intervf(1)=X2;
        
        X1=Field_edge(1,2)-(Field_edge(1,2)-Field_edge(1,1))*portion;
        X2=Field_edge(1,2);
        interv=[X1 X2];
        [Phase_info_2]=PF_sectionphase(Spiketxvtlph,behav.txvtl,interv,BinSize,'mid');
        
        X1=Field_edge(1,2);
        X2=Field_edge(1,2)+(Field_edge(1,3)-Field_edge(1,2))*portion;
        interv=[X1 X2];
        [Phase_info_3]=PF_sectionphase(Spiketxvtlph,behav.txvtl,interv,BinSize,'mid');
        
        X1=Field_edge(1,3)-(Field_edge(1,3)-Field_edge(1,2))*portion;
        X2=Field_edge(1,3);
        interv=[X1 X2];
        [Phase_info_4]=PF_sectionphase(Spiketxvtlph,behav.txvtl,interv,BinSize,'last');
        intervf(2)=X1;
        
        X1=Field_edge(1,1);
        X2=Field_edge(1,3);
        interv=[X1 X2];
        [Phase_info_5,In_Spike]=PF_sectionphase(Spiketxvtlph,behav.txvtl,interv,BinSize,'all');
  
       
        Phase_info=[Phase_info_1;Phase_info_2;Phase_info_3;Phase_info_4;Phase_info_5];
        
        Linear_All=cat(1,Linear_All,Linear);
        Field_edge_All=cat(1,Field_edge_All,Field_edge);
        rows_All=cat(1,rows_All,rows);
        PF_All=cat(1,PF_All,PF);
        Phase_info_All=cat(3,Phase_info_All,Phase_info);
        intervf_All=cat(1,intervf_All,intervf);
        In_Spike_All{1,j}=In_Spike;
        
    end
    
    tiral.Linear=Linear_All;
    tiral.Field_edge=Field_edge_All;
    tiral.rows=rows_All;
    tiral.PF=PF_All;
    tiral.Phase_info=Phase_info_All;
    tiral.intervf=intervf_All; 
%     tiral.spikes=spikes;
    tiral.txvtl=behav.txvtl;
    tiral.In_Spike=In_Spike_All;%3
    tiral.SPK_speed=SPK_speed;
    
end
    
    
    
    
    
    
    









































