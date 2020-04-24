
% delta=0.3;smooth_rate=10;speed_thr=0.04;
function [dataPP_all,statsPP_all]= PhasePrecession_CirularLinear(filename,Rate_Matrix,BinSize,tr,smooth_rate,delta,speed_thr)
%% This function find the spike in the field for a given trial nmubers and calculate linear regression phase info

load([filename '.mat'])
folder='Cellinfo3';

TX=behav.TXVt;
TX(TX(:,4)==0,:)=[];
TX(isnan(TX(:,2))==1,:)=[];
TX=TX(find(TX(:,3)>speed_thr),:);
TX(:,2)=TX(:,2)./max(TX(:,2));
TX=TX(:,1:2);


for j=1:length(G)
    
    
    % j=89
    CellN=j
    
    if Phase.probeID==3
        
        if floor(CellN/100)<7
            
            thetaP(:,2)=Phase.thetaP_1;
        else
            thetaP(:,2)=Phase.thetaP_2;
        end
        
    else
        
        if floor(CellN/100)<9
            
            thetaP(:,2)=Phase.thetaP_1;
        else
            thetaP(:,2)=Phase.thetaP_2;
        end
        
    end
    thetaP(:,1)=(0:length(thetaP)-1)/1250;

    
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
    
    
        
    [Field_edge,rows,PF] = Feild_width_linearMaze(CellN,smooth_rate,delta,Rate_Matrix,trials_port);
    interv=[Field_edge(1,1) Field_edge(1,3)-1];
    [~,In_Spike]=PF_sectionphase(TXVtPh,behav.TXVt,interv,BinSize,'all');
    dataPP=[];statsPP=[];
    
    if size(In_Spike,1)>4
        
    [dataPP,statsPP] = PhasePrecession_F(TX,In_Spike(:,1),thetaP, 'boundaries',[0. .1] );
    dataPP_all{1,j}=dataPP;
    statsPP_all{1,j}=statsPP;
    
    else
        
    dataPP_all{1,j}=nan;
    statsPP_all{1,j}=nan;
    
    end
    
end


