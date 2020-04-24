
%xtts=[xbin TimeSpent trialnb nspike]
%SpeedThreshold in cm/s

function [xtvtls2,xtvtlph]=make_rate_matrix_1D_bz(behav,xbinNumber,SpeedThreshold,Cellinfo)
%%
% ses=2;
% filename=sessions{ses};
% SpeedThreshold=.04;
% nBins=round(max(behav.TXVt(:,2))*100,1)/2;
% SpeedThreshold=.04;
% xbinNumber=round(max(behav.TXVt(:,2))*100,1)/2;
% SpeedThreshold=.04;
% xbinNumber=round(max(behav.TXVt(:,2))*100,1)/2;
%%
%for each belt rotation, compute rate in each space bin

% ncell=size(spikes.maze,2);
ncell=size(Cellinfo,2);
txvtlC=behav.txvtl;
nspike=[];
TimeSpent=[];
xbin=[];
trial=[];
V=[];
laser=[];
mean_Phase=[];
nspike2=[];


binT=min(diff(txvtlC(:,1)));
trialnb=unique(txvtlC(:,4));
trialnb(trialnb==0)=[];
test=[];
T=[];
b=size(txvtlC,2);

%% load Phase
tic
Tph=[];
size(Cellinfo,2)
for ii = 1:size(Cellinfo,2)

    Spiketxvtlph=[];
    Spiketxvtlph=Cellinfo(ii).cell_txvtlph;
%     if isempty(Spiketxvtlph)==0
    Spiketxvtlph(isnan(Spiketxvtlph(:,2))==1,:)=[];
    Spiketxvtlph(find(Spiketxvtlph(:,3)<SpeedThreshold),:)=[];
    Spiketxvtlph(Spiketxvtlph(:,4)==(0),:)=[];   
    
    Tph{ii}(:,1)=Spiketxvtlph(:,1);
    Tph{ii}(:,2)=Spiketxvtlph(:,size(Spiketxvtlph,2));% Define Phase channel
%     end
end
toc

%%
for ii =1:length(trialnb)
    
    ndtrial=[];T_bin=[];
    ndtrial=find(txvtlC(:,4)==trialnb(ii));
    %compute xbin edges
    xstep=max(txvtlC(ndtrial,2))/xbinNumber;
    
    if isnan(xstep)~=1
        
        binedges=0:xstep:max(txvtlC(ndtrial,2));
        xbin=[xbin;(1:xbinNumber)'];
        
        %select speed
        tmpspeed=txvtlC(ndtrial,3);
        tmp=txvtlC(ndtrial,:);
        tmpS=tmp(tmpspeed>=SpeedThreshold,:);
        test=[test; size(tmpS,1) size(tmp,1)];
        
        for jj=1:length(binedges)-1
            
            %spike count and time spent
            ndx=find(tmpS(:,2)>binedges(jj)&tmpS(:,2)<=binedges(jj+1));
            trial=[trial;trialnb(ii)];
            
            
            if ~isempty(ndx)
                
                nspike=[nspike;sum(tmpS(ndx,b+1:end),1)];

                TimeSpent=[TimeSpent;length(ndx)*binT];
                
                V=[V;nanmean(tmpS(ndx,3))];
                
                laser=[laser;sum(tmpS(ndx,5))]; % Add Laser time
                
                Ph=[];N_SP=[];
                for cell=1:ncell;
                    sp=[];
                    sp=find(Tph{cell}(:,1)>= tmpS(ndx(1),1) & Tph{cell}(:,1)< tmpS(ndx(end),1));
                    if ~isempty(sp)
                        [~,Ph(1,cell)]=meanphase(Tph{cell}(sp,2));
                        N_SP(1,cell)=length(sp);
                    else
                        Ph(1,cell)=nan;
                        N_SP(1,cell)=0;
                    end
                end
                
                mean_Phase=[mean_Phase; Ph];
                nspike2=[nspike2; N_SP];
                %             T_bin(jj,2)=tmpS(ndx(end),1);
                %             if jj==1
                %                 T_bin(jj,1)=tmpS(ndx(1),1);
                %             else
                %                 T_bin(jj,1)=T_bin(jj-1,2);
                %             end
            else
                
                nspike=[nspike;zeros(1,ncell)];
                
                TimeSpent=[TimeSpent;0.000001];
                
                V=[V;nan];
                
                laser=[laser;nan]; % Add Laser time
                
                mean_Phase=[mean_Phase; nan(1,ncell)];
                
                nspike2=[nspike2; zeros(1,ncell)];
                %             T_bin(jj,2)=nan;
                %
                %             if jj==1
                %                 T_bin(jj,1)=nan;
                %             else
                %                 T_bin(jj,1)=nan;
                %             end
                
            end
            
        end
        
        %     T=cat(1,T,T_bin);
        
    else
        xbin=[xbin;nan(xbinNumber,1)] ;
        trial=[trial;nan(xbinNumber,1)];
        laser=[laser;nan(xbinNumber,1)] ;
        V=[V;nan(xbinNumber,1)];
        TimeSpent=[TimeSpent;nan(xbinNumber,1)] ;
        mean_Phase=[mean_Phase;nan(xbinNumber,ncell)] ;
        nspike=[nspike;nan(xbinNumber,ncell)];
        nspike2=[nspike2;nan(xbinNumber,ncell)] ;
        
%         size(trial)
%         size(laser)
%         size(xbin)
%         size(V )
%         size(mean_Phase)
%         size(TimeSpent )
%         size(nspike )
%         size(nspike2 )
        
        
        
    end
end

%%
% xtvts1=[xbin TimeSpent V trial laser nspike];
xtvtls2=[xbin TimeSpent V trial laser nspike2];
xtvtlph=[xbin TimeSpent V trial laser mean_Phase];

%%
% test1=find(isnan(mean_Phase(:,61))==0);
% test2=find(nspike2(:,61)>0);
% sum(test1==test2)==length(test1)


%%


