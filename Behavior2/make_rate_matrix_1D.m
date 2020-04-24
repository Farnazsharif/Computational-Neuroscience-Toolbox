
%xtts=[xbin TimeSpent trialnb nspike]
%SpeedThreshold in cm/s

function [xtvts1,xtvts2,xtvtph]=make_rate_matrix_1D(filename,xbinNumber,SpeedThreshold,folder)
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
load([filename '.mat'])
txvtlC=behav.txvtl_rate;  
ncell=length(unique(G));

%%
%for each belt rotation, compute rate in each space bin

nspike=[];
TimeSpent=[];
xbin=[];
trial=[];
V=[];

mean_Phase=[];
nspike2=[];

binT=min(diff(txvtlC(:,1)));
trialnb=unique(txvtlC(:,4));
trialnb(trialnb==0)=[];
test=[];
T=[];
b=size(txvtlC,2)-ncell;

%% load Phase
tic
Tph=[];
cd (folder)
for ii = 1:ncell
    
    load( ['TXVtPh_' num2str(ii)])
    
    TXVtPh(isnan(TXVtPh(:,2))==1,:)=[];
    TXVtPh(find(TXVtPh(:,3)<SpeedThreshold),:)=[];
    TXVtPh(TXVtPh(:,4)==(0),:)=[]; 
    
    Tph{ii}(:,1)=TXVtPh(:,1);
    Tph{ii}(:,2)=TXVtPh(:,5);
end
cd ..
toc

%%
for ii =1:length(trialnb)
    
    ndtrial=[];T_bin=[];
    ndtrial=find(txvtlC(:,4)==trialnb(ii));
    %compute xbin edges
    xstep=max(txvtlC(ndtrial,2))/xbinNumber;
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
            clc
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

end
%%
xtvts1=[xbin TimeSpent V trial  nspike];
xtvts2=[xbin TimeSpent V trial  nspike2];
xtvtph=[xbin TimeSpent V trial  mean_Phase];

%%
% test1=find(isnan(mean_Phase(:,61))==0);
% test2=find(nspike2(:,61)>0);
% sum(test1==test2)==length(test1)











