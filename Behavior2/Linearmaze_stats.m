function [tiral]=Linearmaze_stats(filename,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,BinSize,speed_thr,delta,tr,folder)
%%
    
portion=0.5;
% load([filename '_sessInfo.mat'])
load([filename '.mat'])


%############################################################################################################################
    Linear_Cric=[];Field_edge_Cric=[];rows_Cric=[];PF_Cric=[];Phase_info_Cric=[];intervf_Cric=[];In_Spike_Cric=[];SPK_speed=[];%2

    for j=1:length(G)

% j=89
        CellN=j;
       
        cd (folder)
        load(['TXVtPh_' num2str(CellN) '.mat'])
        cd ..
        
      
        TXVtPh(isnan(TXVtPh(:,2))==1,:)=[];
        TXVtPh(find(TXVtPh(:,3)<speed_thr),:)=[];%.04
        TXVtPh(TXVtPh(:,4)==(0),:)=[];
        
        trials_port=tr;
        nd_int=[];
        for i=1:length(trials_port)
            nd_int=[nd_int;find(TXVtPh(:,4)==trials_port(i))];
        end
        TXVtPh=TXVtPh(nd_int,:);
        
        
        SPK_speed=[SPK_speed;nanmean(TXVtPh(:,3))];
      
       
        Linear=[];Field_edge=[];rows=[];PF=[];Phase_info=[];intervf=[];In_Spike=[];%1
        
        [Linear] = linear_maze_PFinfo(CellN,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,tr);

        [Field_edge,rows,PF] = Feild_width_linearMaze(CellN,smooth_rate,delta,Rate_Matrix,tr);
        
        Phase_info_1=[];Phase_info_2=[];Phase_info_3=[];Phase_info_4=[];
        
        X2=Field_edge(1,1)+(Field_edge(1,2)-Field_edge(1,1))*portion;
        X1=Field_edge(1,1);
        interv=[X1 X2];
        [Phase_info_1]=PF_sectionphase(TXVtPh,behav.TXVt,interv,BinSize,'first');
        intervf(1)=X2;
        
        X1=Field_edge(1,2)-(Field_edge(1,2)-Field_edge(1,1))*portion;
        X2=Field_edge(1,2);
        interv=[X1 X2];
        [Phase_info_2]=PF_sectionphase(TXVtPh,behav.TXVt,interv,BinSize,'mid');
        
        X1=Field_edge(1,2);
        X2=Field_edge(1,2)+(Field_edge(1,3)-Field_edge(1,2))*portion;
        interv=[X1 X2];
        [Phase_info_3]=PF_sectionphase(TXVtPh,behav.TXVt,interv,BinSize,'mid');
        
        X1=Field_edge(1,3)-(Field_edge(1,3)-Field_edge(1,2))*portion;
        X2=Field_edge(1,3);
        interv=[X1 X2];
        [Phase_info_4]=PF_sectionphase(TXVtPh,behav.TXVt,interv,BinSize,'last');
        intervf(2)=X1;
        
        X1=Field_edge(1,1);
        X2=Field_edge(1,3);
        interv=[X1 X2];
        [Phase_info_5,In_Spike]=PF_sectionphase(TXVtPh,behav.TXVt,interv,BinSize,'all');
  
       
        Phase_info=[Phase_info_1;Phase_info_2;Phase_info_3;Phase_info_4;Phase_info_5];
        
        Linear_Cric=cat(1,Linear_Cric,Linear);
        Field_edge_Cric=cat(1,Field_edge_Cric,Field_edge);
        rows_Cric=cat(1,rows_Cric,rows);
        PF_Cric=cat(1,PF_Cric,PF);
        Phase_info_Cric=cat(3,Phase_info_Cric,Phase_info);
        intervf_Cric=cat(1,intervf_Cric,intervf);
        In_Spike_Cric{1,j}=In_Spike;
        
    end
    
    tiral.Linear=Linear_Cric;
    tiral.Field_edge=Field_edge_Cric;
    tiral.rows=rows_Cric;
    tiral.PF=PF_Cric;
    tiral.Phase_info=Phase_info_Cric;
    tiral.intervf=intervf_Cric; 
    tiral.G=G;
    tiral.Cellgroup=Cellgroup;
    tiral.TXVt=behav.TXVt;
    tiral.In_Spike=In_Spike_Cric;%3
    tiral.SPK_speed=SPK_speed;
    
end
    
    
    
    
    
    
    













































% % delta=0.15;
% % speed_thr=0.04;
% % Rate_Matrix=xtvts2;
% % Phase_Matrix=xtvtph;
% % smooth_rate=10;
% % smooth_phase=1;
% % BinSize=(max(behav.TXVt(:,2))./max(xtvts1(:,1)));
% portion=0.5;
% load([filename '_sessInfo.mat'])
% load([filename '.mat'])
% 
% %%
% if isempty(strfind(sessInfo.Position.MazeType,'Linear'))==0
%     
%     
% 
%     
%     Linear_odd=[];Field_edge_odd=[];rows_odd=[];PF_odd=[];Phase_info_odd=[];intervf_odd=[];In_Spike_odd=[];
%     Linear_even=[];Field_edge_even=[];rows_even=[];PF_even=[];Phase_info_even=[];intervf_even=[];In_Spike_even=[];
%     
%     tic
%     
%     for j=1:length(G)
%         CellN=j;
% 
%         cd (folder)
%         load(['TXVtPh_' num2str(CellN) '.mat'])
%         cd ..
%         
%         TXVtPh(isnan(TXVtPh(:,2))==1,:)=[];
%         TXVtPh(find(TXVtPh(:,3)<speed_thr),:)=[];%.04
%         TXVtPh(TXVtPh(:,4)==(0),:)=[];
%         
%         Nd_even=[];Nd_odd=[];
%         Nd_even=find(mod(TXVtPh(:,4),2)==0);
%         Nd_odd=find(mod(TXVtPh(:,4),2)==1);
%         
%         
%         Linear=[];Field_edge=[];rows=[];PF=[];Phase_info=[];intervf=[];In_Spike=[];%1
%         
%         [Linear] = linear_maze_PFinfo(CellN,Rate_Matrix_odd,Phase_Matrix_odd,smooth_rate,smooth_phase,tr);
%         
%         [Field_edge,rows,PF] = Feild_width_linearMaze(CellN,smooth_rate,delta,Rate_Matrix_odd,tr);
%         
%         
%         
%         Phase_info_1=[];Phase_info_2=[];Phase_info_3=[];Phase_info_4=[];
%         X2=Field_edge(1,1)+(Field_edge(1,2)-Field_edge(1,1))*portion;
%         X1=Field_edge(1,1);
%         interv=[X1 X2];
%         [Phase_info_1]=PF_sectionphase(TXVtPh(Nd_odd,:),behav.TXVt,interv,BinSize);
%         intervf(1)=X2;
%         
%         X1=Field_edge(1,2)-(Field_edge(1,2)-Field_edge(1,1))*portion;
%         X2=Field_edge(1,2);
%         interv=[X1 X2];
%         [Phase_info_2]=PF_sectionphase(TXVtPh(Nd_odd,:),behav.TXVt,interv,BinSize);
%         
%         X1=Field_edge(1,2);
%         X2=Field_edge(1,2)+(Field_edge(1,3)-Field_edge(1,2))*portion;
%         interv=[X1 X2];
%         [Phase_info_3]=PF_sectionphase(TXVtPh(Nd_odd,:),behav.TXVt,interv,BinSize);
%         
%         X1=Field_edge(1,3)-(Field_edge(1,3)-Field_edge(1,2))*portion;
%         X2=Field_edge(1,3);
%         interv=[X1 X2];
%         [Phase_info_4]=PF_sectionphase(TXVtPh(Nd_odd,:),behav.TXVt,interv,BinSize);
%         intervf(2)=X1;
%         
%         X1=Field_edge(1,1);
%         X2=Field_edge(1,3);
%         interv=[X1 X2];
%         [Phase_info_5,In_Spike]=PF_sectionphase(TXVtPh(Nd_odd,:),behav.TXVt,interv,BinSize);
%         
%         
%         Phase_info=[Phase_info_1;Phase_info_2;Phase_info_3;Phase_info_4;Phase_info_5];
%         
%         
%         
%         Linear_odd=cat(1,Linear_odd,Linear);
%         Field_edge_odd=cat(1,Field_edge_odd,Field_edge);
%         rows_odd=cat(1,rows_odd,rows);
%         PF_odd=cat(1,PF_odd,PF);
%         Phase_info_odd=cat(3,Phase_info_odd,Phase_info);
%         intervf_odd=cat(1,intervf_odd,intervf);
%         In_Spike_odd{1,j}=In_Spike;
%         
%         Linear=[];Field_edge=[];rows=[];PF=[];Phase_info=[];intervf=[];In_Spike=[];%1
%         
%         [Linear] = linear_maze_PFinfo(CellN,Rate_Matrix_even,Phase_Matrix_even,smooth_rate,smooth_phase,tr);
%         [Field_edge,rows,PF] = Feild_width_linearMaze(CellN,smooth_rate,delta,Rate_Matrix_even,tr);
%         
%         Phase_info_1=[];Phase_info_2=[];Phase_info_3=[];Phase_info_4=[];
%         intervf(1)=X2;
%         
%         X2=Field_edge(1,1)+(Field_edge(1,2)-Field_edge(1,1))*portion;
%         X1=Field_edge(1,1);
%         interv=[X1 X2];
%         [Phase_info_1]=PF_sectionphase(TXVtPh(Nd_even,:),behav.TXVt,interv,BinSize);
%         
%         X1=Field_edge(1,2)-(Field_edge(1,2)-Field_edge(1,1))*portion;
%         X2=Field_edge(1,2);
%         interv=[X1 X2];
%         [Phase_info_2]=PF_sectionphase(TXVtPh(Nd_even,:),behav.TXVt,interv,BinSize);
%         
%         X1=Field_edge(1,2);
%         X2=Field_edge(1,2)+(Field_edge(1,3)-Field_edge(1,2))*portion;
%         interv=[X1 X2];
%         [Phase_info_3]=PF_sectionphase(TXVtPh(Nd_even,:),behav.TXVt,interv,BinSize);
%         
%         X1=Field_edge(1,3)-(Field_edge(1,3)-Field_edge(1,2))*portion;
%         X2=Field_edge(1,3);
%         interv=[X1 X2];
%         [Phase_info_4]=PF_sectionphase(TXVtPh(Nd_even,:),behav.TXVt,interv,BinSize);
%         intervf(2)=X1;
%         
%         X1=Field_edge(1,1);
%         X2=Field_edge(1,3);
%         interv=[X1 X2];
%         [Phase_info_5,In_Spike]=PF_sectionphase(TXVtPh(Nd_even,:),behav.TXVt,interv,BinSize);
%         
%         
%         Phase_info=[Phase_info_1;Phase_info_2;Phase_info_3;Phase_info_4;Phase_info_5];
%         
%         
%         Linear_even=cat(1,Linear_even,Linear);
%         Field_edge_even=cat(1,Field_edge_even,Field_edge);
%         rows_even=cat(1,rows_even,rows);
%         PF_even=cat(1,PF_even,PF);
%         Phase_info_even=cat(3,Phase_info_even,Phase_info);
%         intervf_even=cat(1,intervf_even,intervf);
%         In_Spike_even{1,j}=In_Spike;
%     end
%     
%     tiral.Linear_even=Linear_even;
%     tiral.Linear_odd=Linear_odd;
%     tiral.Linear_Cric=nan(1, size(Linear_odd,2));
%     
%     tiral.Field_edge_even=Field_edge_even;
%     tiral.Field_edge_odd=Field_edge_odd;
%     tiral.Field_edge_Cric=nan(1, size(Field_edge_odd,2));
%     
%     tiral.rows_even=rows_even;
%     tiral.rows_odd=rows_odd;
%     tiral.rows_Cric=nan(1, size(rows_odd,2));
%     
%     tiral.PF_even=PF_even;
%     tiral.PF_odd=PF_odd;
%     tiral.PF_Cric=nan(1,size(PF_even,2), size(PF_even,3));
%     
%     tiral.Phase_info_even=Phase_info_even;
%     tiral.Phase_info_odd=Phase_info_odd;
%     tiral.Phase_info_Cric=nan(1,size(Phase_info_Cric,2),size(Phase_info_Cric,3));
%     
%     tiral.intervf_Cric=nan(1,2);
%     tiral.intervf_even=intervf_even;
%     tiral.intervf_odd=intervf_odd;
%     
%     tiral.In_Spike_odd=In_Spike_odd;%3
%     tiral.In_Spike_even=In_Spike_even;%3
%     
%     tiral.G=G;
%     tiral.Cellgroup=Cellgroup;
%     tiral.TXVt=behav.TXVt;
%     
% else
%     %############################################################################################################################
%     Linear_Cric=[];Field_edge_Cric=[];rows_Cric=[];PF_Cric=[];Phase_info_Cric=[];intervf_Cric=[];In_Spike_Cric=[];SPK_speed=[];%2
% 
%     for j=1:length(G)
% 
%         CellN=j;
%        
%         cd (folder)
%         load(['TXVtPh_' num2str(CellN) '.mat'])
%         cd ..
%         
%       
%         TXVtPh(isnan(TXVtPh(:,2))==1,:)=[];
%         TXVtPh(find(TXVtPh(:,3)<speed_thr),:)=[];%.04
%         TXVtPh(TXVtPh(:,4)==(0),:)=[];
% 
%         trials_port=[tr(1):tr(2)];
%         nd_int=[];
%         for i=1:length(trials_port)
%             nd_int=[nd_int;find(TXVtPh(:,4)==trials_port(i))];
%         end
%    
%         TXVtPh=TXVtPh(nd_int,:);
%         SPK_speed=[SPK_speed;nanmean(TXVtPh(:,3))];
%       
%        
%         Linear=[];Field_edge=[];rows=[];PF=[];Phase_info=[];intervf=[];In_Spike=[];%1
%         
%         [Linear] = linear_maze_PFinfo(CellN,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,tr);
%         
%         [Field_edge,rows,PF] = Feild_width_linearMaze(CellN,smooth_rate,delta,Rate_Matrix,tr);
%         
%         Phase_info_1=[];Phase_info_2=[];Phase_info_3=[];Phase_info_4=[];
%         
%         X2=Field_edge(1,1)+(Field_edge(1,2)-Field_edge(1,1))*portion;
%         X1=Field_edge(1,1);
%         interv=[X1 X2];
%         [Phase_info_1]=PF_sectionphase(TXVtPh,behav.TXVt,interv,BinSize,'first');
%         intervf(1)=X2;
%         
%         X1=Field_edge(1,2)-(Field_edge(1,2)-Field_edge(1,1))*portion;
%         X2=Field_edge(1,2);
%         interv=[X1 X2];
%         [Phase_info_2]=PF_sectionphase(TXVtPh,behav.TXVt,interv,BinSize,'mid');
%         
%         X1=Field_edge(1,2);
%         X2=Field_edge(1,2)+(Field_edge(1,3)-Field_edge(1,2))*portion;
%         interv=[X1 X2];
%         [Phase_info_3]=PF_sectionphase(TXVtPh,behav.TXVt,interv,BinSize,'mid');
%         
%         X1=Field_edge(1,3)-(Field_edge(1,3)-Field_edge(1,2))*portion;
%         X2=Field_edge(1,3);
%         interv=[X1 X2];
%         [Phase_info_4]=PF_sectionphase(TXVtPh,behav.TXVt,interv,BinSize,'last');
%         intervf(2)=X1;
%         
%         X1=Field_edge(1,1);
%         X2=Field_edge(1,3);
%         interv=[X1 X2];
%         [Phase_info_5,In_Spike]=PF_sectionphase(TXVtPh,behav.TXVt,interv,BinSize,'all');
%   
%        
%         Phase_info=[Phase_info_1;Phase_info_2;Phase_info_3;Phase_info_4;Phase_info_5];
%         
%         Linear_Cric=cat(1,Linear_Cric,Linear);
%         Field_edge_Cric=cat(1,Field_edge_Cric,Field_edge);
%         rows_Cric=cat(1,rows_Cric,rows);
%         PF_Cric=cat(1,PF_Cric,PF);
%         Phase_info_Cric=cat(3,Phase_info_Cric,Phase_info);
%         intervf_Cric=cat(1,intervf_Cric,intervf);
%         In_Spike_Cric{1,j}=In_Spike;
%         
%     end
%     
%     tiral.Linear_Cric=Linear_Cric;
%     tiral.Linear_even=nan(1, size(Linear_Cric,2));
%     tiral.Linear_odd=nan(1, size(Linear_Cric,2));
%     
%     tiral.Field_edge_Cric=Field_edge_Cric;
%     tiral.Field_edge_even=nan(1, size(Field_edge_Cric,2));
%     tiral.Field_edge_odd=nan(1, size(Field_edge_Cric,2));
%     
%     tiral.rows_Cric=rows_Cric;
%     tiral.rows_even=nan(1, size(rows_Cric,2));
%     tiral.rows_odd=nan(1, size(rows_Cric,2));
%     
%     tiral.PF_Cric=PF_Cric;
%     tiral.PF_even=nan(1,size(PF_Cric,2), size(PF_Cric,3));
%     tiral.PF_odd=nan(1,size(PF_Cric,2), size(PF_Cric,3));
%     
%     tiral.Phase_info_Cric=Phase_info_Cric;
%     tiral.Phase_info_even=nan(1,size(Phase_info_Cric,2),size(Phase_info_Cric,3));
%     tiral.Phase_info_odd=nan(1,size(Phase_info_Cric,2),size(Phase_info_Cric,3));
%     
%     tiral.intervf_Cric=intervf_Cric;
%     tiral.intervf_even=nan(1,2);
%     tiral.intervf_odd=nan(1,2);
%     
%     tiral.G=G;
%     tiral.Cellgroup=Cellgroup;
%     tiral.TXVt=behav.TXVt;
%     tiral.In_Spike_Cric=In_Spike_Cric;%3
%     tiral.SPK_speed=SPK_speed;




