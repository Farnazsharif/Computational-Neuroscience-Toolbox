function  [Antonio_all_Field,Antonio_Infield,Olypher,SKW]=Antonio_cellinfo(CellVector,smooth,Matrix,Field_edge,TR,TRX,filename,folder,nPFbin)
% %%
% OUTPUTS:
% 1- Information Index 1 = bits per spike
% 2- Information Index 2 = bits per second
% 3- Sparsity =
% 4- Coefficient = index to quantify the sparsity (I am using directly the
%                  sparsity not this one)
% 5- Selectivity
% 6- Fano factor of firing peak
% 7- firing rate tunity
% 8- width
% Q={'Information Index 1 = bits per spike','Information Index 2 = bits per second','Sparsity',' Coefficient','Selectivity','Fano factor of firing peak','firing rate tunity','PVR','width','TunityR','Pvalue','MeanFR-OUT','MeanFR-All','MeanFR-IN'};

%%
% clear
% clc
% filename='SFA4_S3_TRD2';
% load([filename '.mat']);
% load([filename 'PlaceField.mat'])
% Cellvector=Cell_group.C(:,3);
% smooth=1;
% Matrix=xttsc;
% tr_TRD1=[[1 6 7];[1 9 7];[1 11 7];[1 12 7];[1 3 7];[1 8 7]];
% tr_TRD2=[[1    21    7];[ 1    19    7];[ 1    18    7];[1    11    7];[1    24    7];[1    22    7];[1    21    7]];
% TRX=[1    18    43;  44    67    86];
% tr_1=[1    21    45]
speed_thr=5;
XV=0:2*pi/(nPFbin-1):2*pi;
Antonio_Infield=[];
round_to=2;

%%

clc
[a,~]=size(CellVector);

for i=1:a

    trial=[];
    if isnan(Field_edge(i,1))==1 || isnan(Field_edge(i,3))==1
        
        Antonio_Infield(i,:)=nan(1,11);
        Antonio_all_Field(i,:)=nan(1,15);
        Olypher(i,:)=nan(1,4);
        SKW(i,1)=nan(1,1);
    else
        
        Celln=CellVector(i,3);
        smoothT=smooth1D(Matrix(:,2),smooth,1);
        smoothC=smooth1D(Matrix(:,Celln+4),smooth,1);
        rate=smoothC./repmat(smoothT,1,1);
        
        matC=reshape(rate,nPFbin,length(rate)/nPFbin)';
        snSpk=reshape(smoothC,nPFbin,length(smoothC)/nPFbin)';
        sTime=reshape(smoothT,nPFbin,length(smoothT)/nPFbin)';

        cd (folder)
        load(['TXTSVPh_' num2str(i)])
        cd ..
        [Fmap]=phasemap(filename,TXTSVPh,speed_thr,nPFbin);

        if strfind(filename,'SFA5_S4_TRD2')==1
            
            if Celln == 53||Celln ==163||Celln ==20||Celln ==26||Celln ==113 || Celln ==13 || Celln == 115||Celln ==164||Celln ==232  || Celln ==39 || Celln == 26 ||Celln ==83 ||Celln ==115 % TRD2
                trial=TRX(2,:);
                matC=matC(trial(1):trial(3),:);
                snSpk=snSpk(trial(1):trial(3),:);
                sTime=sTime(trial(1):trial(3),:);
                Fmap=Fmap(trial(1):trial(3),:,:);
            else
                trial=TRX(1,:);
                matC=matC(trial(1):trial(3),:);
                snSpk=snSpk(trial(1):trial(3),:);
                sTime=sTime(trial(1):trial(3),:);
                Fmap=Fmap(trial(1):trial(3),:,:);
            end
            
        elseif strfind(filename,'SFA3_S3_TRD1')==1 % allways, no matter eqyal or nonequal
            trial=[1 3 7];
            
            matC=matC(trial(1):trial(3),:);
            snSpk=snSpk(trial(1):trial(3),:);
            sTime=sTime(trial(1):trial(3),:);
            Fmap=Fmap(trial(1):trial(3),:,:);
        else
            trial=TR;
            
            matC=matC(trial(1):trial(3),:);
            snSpk=snSpk(trial(1):trial(3),:);
            sTime=sTime(trial(1):trial(3),:);
            Fmap=Fmap(trial(1):trial(3),:,:);
            
        end
        
       
        
        Rate_Map=[];snSpikes=[];sTimeSpent=[];
        
        MR=[matC,matC,matC];
        MS=[snSpk,snSpk,snSpk];
        MT=[sTime,sTime,sTime];
        
        Rate_Map=MR(:,Field_edge(i,1):Field_edge(i,3));
        snSpikes=MS(:,Field_edge(i,1):Field_edge(i,3));
        
        sTimeSpent=MT(:,Field_edge(i,1):Field_edge(i,3));
        
         SKW(i,1)=skewness(mean(Rate_Map));
        
        [Information_1,Information_2,Sparsity,Coefficient,Selectivity] = PlaceCellInfo_Antonio(Rate_Map, snSpikes, sTimeSpent);
        Antonio_Infield(i,1)=Information_1;
        Antonio_Infield(i,2)=Information_2;
        Antonio_Infield(i,3)=Sparsity;
        Antonio_Infield(i,4)=Coefficient;
        Antonio_Infield(i,5)=Selectivity;
        n=[];
        for n = 1:size(Rate_Map,1)
            [~,RM(n,1)] = max(Rate_Map(n,:));
        end
        Antonio_Infield(i,6)= var(RM)/ mean(RM);
        
        
        %         PFo_spk=sum(Rate_Map,1);
        %         A=sum(cos(XV).*PFo_spk);
        %         B=sum(sin(XV).*PFo_spk);
        %         Antonio_Infield(i,7)=sqrt(A^2+B^2)/sum(PFo_spk);
        %         if isnan(Antonio_Infield(i,7))==1
        %         Antonio_Infield(i,7)=0;
        %         end
        
        PF1   =  mean(Rate_Map(1:round(size(Rate_Map,1)/2),:),1);
        PF2   =  mean(Rate_Map(round(size(Rate_Map,1)/2)+1:end,:),1);
        [R,P] =  corrcoef(PF1,PF2);
        
        %         if  P(1,2)>0.08 || isnan(R(1,2))==1
        if   isnan(R(1,2))==1
            R(1,2)=0;
        end
        Antonio_Infield(i,8)=R(1,2);
        
        
        Antonio_Infield(i,9)=length(round(Field_edge(i,1)):round(Field_edge(i,3)));
        
        %         PF_xt=[];theta_run=[];
        %         [PF_xt,~] = PF_info(filename,Cellvector(i,4),Field_edge(i,:),speed_thr);  % get the spk during the run
        %         theta_run(:,1)=PF_xt(:,1)*(2*pi/230);
        %         Antonio_Infield(i,10)=[meanphase(theta_run(:,1))];
        %         Antonio_Infield(i,11)=[circ_rtest(theta_run(:,1))];
        
        XTPH=[];theta_PF=[];
        cd (folder)
        load(['XTPhT_' num2str(i)])
        cd ..
        if isempty(XTPH)==1
            Spk_tunity=[nan nan nan nan ];
        else
            theta_PF(:,1)=XTPH(:,1)*(2*pi);
            Antonio_Infield(i,10)=[meanphase(theta_PF(:,1))];
            Antonio_Infield(i,11)=[circ_rtest(theta_PF(:,1))];
        end
        Antonio_all_Field(i,14)=mean(mean(Rate_Map));
        
        Rate_Map=[];snSpikes=[];sTimeSpent=[];
        
        Rate_Map=matC;
        snSpikes=snSpk;
        sTimeSpent=sTime;
        
        [Information_1,Information_2,Sparsity,Coefficient,Selectivity] = PlaceCellInfo_Antonio(Rate_Map, snSpikes, sTimeSpent);
        Antonio_all_Field(i,1)=Information_1;
        Antonio_all_Field(i,2)=Information_2;
        Antonio_all_Field(i,3)=Sparsity;
        Antonio_all_Field(i,4)=Coefficient;
        Antonio_all_Field(i,5)=Selectivity;
        
        n=[];
        for n = 1:size( Rate_Map,1)
            [~,RM(n,1)] = max(Rate_Map(n,:));
        end
        Antonio_all_Field(i,6)= var(RM)/ mean(RM);
        
        PFo_spk=sum(Rate_Map,1);
        A=sum(cos(XV).*PFo_spk);
        B=sum(sin(XV).*PFo_spk);
        Antonio_all_Field(i,7)=sqrt(A^2+B^2)/sum(PFo_spk);
        if isnan(Antonio_all_Field(i,7))==1
            Antonio_all_Field(i,7)=0;
        end
        
        PF1   =  mean(Rate_Map(1:round(size(Rate_Map,1)/2),:),1);
        PF2   =  mean(Rate_Map(round(size(Rate_Map,1)/2)+1:end,:),1);
        
        [R,P] = corrcoef(PF1,PF2);
        %         if  P(1,2)>0.08 || isnan(R(1,2))==1
        if   isnan(R(1,2))==1
            R(1,2)=0;
        end
        
        Antonio_all_Field(i,8)=R(1,2);
        
        Antonio_all_Field(i,9)=length(round(Field_edge(i,1)):round(Field_edge(i,3)));
        
%         run_xt=[];theta_run=[];
%         [~,run_xt] = PF_info(filename,CellVector(i,4),Field_edge(i,:),speed_thr,nPFbin);  % get the spk during the run
%         theta_run(:,1)=run_xt(:,1)*(2*pi/230);
%         Antonio_all_Field(i,10)=[meanphase(theta_run(:,1))];
%         Antonio_all_Field(i,11)=[circ_rtest(theta_run(:,1))];
        
        AA=[];XX=[];
        AA=[round(Field_edge(i,4)):round(Field_edge(i,5)) round(Field_edge(i,6)):round(Field_edge(i,7))];      
        AA((AA==0))=[];
        XX=matC;
        XX(:,AA)=[];

        Antonio_all_Field(i,12)=mean(mean(XX));
        Antonio_all_Field(i,13)=mean(mean(matC));
        
        Antonio_all_Field(i,15)=max(mean(matC));

%         Phasemap=[];
%         Phasemap(:,:)=Fmap(:,:,2);
%         Phasemap(isnan(Phasemap)==1)=0;
%         [Olypher(i,1),Olypher(i,2)] = PlaceCellInfo_Antonio(Phasemap, snSpikes, sTimeSpent);

        data=matC;
        [track_info_FR,pos_info_val_FR] = bz_olypherInfo_Farn(data,round_to,smooth,'rate');
        Olypher(i,1)=max(track_info_FR);Olypher(i,2)=max(max(pos_info_val_FR));
        matP=Fmap(:,:,2);
        data=matP;
        [track_info_ph,pos_info_val_ph] = bz_olypherInfo_Farn(data,round_to,smooth,'phase');
        Olypher(i,3)=max(track_info_ph);Olypher(i,4)=max(max(pos_info_val_ph));
    end
    
    
end



