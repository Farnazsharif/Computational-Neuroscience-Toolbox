function [midg_tr,Lowg_tr]=count_Gamma(Ndx_g_lg,Ndx_g_hg,r,TXVt,TRN)
    
% CA3_Phase=[degtorad(45) degtorad(160)];
% EC3_Phase=[degtorad(200) degtorad(360)];
%   TXVt= Trial_all{1, 1}{1, 9}.TXVt;  
    %%
    
    [~,C_L]=setdiff(Ndx_g_lg(:,1),Ndx_g_hg(:,1));%delet mid gamma
    [~,C_h]=setdiff(Ndx_g_hg(:,1),Ndx_g_lg(:,1));%delet low gamma
    Ndx_g_lg=Ndx_g_lg(C_L,:);
    Ndx_g_hg=Ndx_g_hg(C_h,:);

    d=[];M=[];
    M=Ndx_g_hg;
    for i=1:size(Ndx_g_lg,1)
        idx=rangesearch(M(:,1),Ndx_g_lg(i,1),r);
        d=[d idx{1}];
    end
    M(d,:)=[];
  
    
%     ph_g=EC3_Phase;
%     ndx_ph=find(M(:,3)>= ph_g(1,1) & M(:,3)<= ph_g(1,2));
%     M=M(ndx_ph,:);


ND_TR=[];
for jj=1:length(TRN)
    nd=[];
    nd=find(M(:,6)==TRN(jj));
    ND_TR(jj,1)=TRN(jj);
    ND_TR(jj,2)=length(nd);
    ND_TR(jj,3)=mean(M(nd,2));
    f=[];[~,f]=meanphase(M(nd,3));
    ND_TR(jj,4)=f;
    nd=[];nd=find(TXVt(:,4)==TRN(jj));
    ND_TR(jj,5)=(TXVt(nd(end),1)-TXVt(nd(1),1));
end

midg_tr=ND_TR;
    
M=[];d=[];
M=Ndx_g_lg;
for i=1:size(Ndx_g_hg,1)
        idx=rangesearch(M(:,1),Ndx_g_hg(i,1),r);
        d=[d idx{1}];
end

% size(M)
M(d,:)=[];
% size(M)

% ph_g=EC3_Phase;
% ndx_ph=find(M(:,3)>= ph_g(1,1) & M(:,3)<= ph_g(1,2));
% M=M(ndx_ph,:);


ND_TR=[];
for jj=1:length(TRN)
    nd=[];
    nd=find(M(:,6)==TRN(jj));
    ND_TR(jj,1)=TRN(jj);
    ND_TR(jj,2)=length(nd);
    ND_TR(jj,3)=mean(M(nd,2));
    f=[];[~,f]=meanphase(M(nd,3));
    ND_TR(jj,4)=f;
    nd=[];nd=find(TXVt(:,4)==TRN(jj));
    ND_TR(jj,5)=(TXVt(nd(end),1)-TXVt(nd(1),1));
end

Lowg_tr=ND_TR;











