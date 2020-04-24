
function [Phase_info,In_Spike_C,In_Spike_L]=PF_sectionphase(TXVtPh,TXVt,interv,BinSize,Section)
% X2=Field_edge(1,1)+(Field_edge(1,2)-Field_edge(1,1))/2
% X1=Field_edge(1,1)
% BinSize=(max(behav.TXVt(:,2))./max(xtvts1(:,1)))

%%
b=size(TXVtPh,2);
M=[];M1=[];M2=[];M3=[];X2=[];X3=[];
M1=TXVtPh;
X2=[TXVtPh(:,2)+max(TXVt(:,2))];
M2=[TXVtPh(:,1) X2  TXVtPh(:,3:b)];
X3=[TXVtPh(:,2)+max(TXVt(:,2))*2];
M3=[TXVtPh(:,1) X3  TXVtPh(:,3:b)];
M=[M1;M2;M3];

%%
% close all

% figure
% plot(M(:,2),M(:,1),'.r','markersize',10)
% hold on
% plot([interv(1)*BinSize interv(1)*BinSize],[min(M(:,1)) max(M(:,1))])
% hold on
% plot([interv(2)*BinSize interv(2)*BinSize],[min(M(:,1)) max(M(:,1))])
% 
% figure
% M2=[rows rows rows];
% M3=(1:length(M2)).*BinSize;
% plot(M3,M2)
% hold on
% plot([interv(1)*BinSize interv(1)*BinSize],[0 max(rows)])
% hold on
% plot([interv(2)*BinSize interv(2)*BinSize],[0 max(rows)])
% M100(:,:)=PF(1,:,:);
% %%
% figure
% imagesc(M100(46:50,:))
% figure
% plot(mean(M100(46:50,:)))


%%
B1=interv(1)*BinSize;
B2=interv(2)*BinSize;
nd=find(M(:,2)>=B1 & M(:,2)<=B2);

%%
if isempty(nd)==1 
    [S(:,1) S(:,2)]=sort(M2(:,2),'ascend'); 
    SpNum=round(.2*size(M2,1));
     nd_s=[];
 if strfind(Section,'last')==1
     nd_s=size(M2,1)-SpNum+1:size(M2,1);
      nd=S(nd_s,2);
 elseif strfind(Section,'first')==1
    nd_s=1:SpNum;
    nd=S(nd_s,2);
 end 
end

%%
O=M(nd,b);
In_Spike_C=M(nd,:);
In_Spike_L=TXVtPh(nd,:);
[meanR,meanO]=meanphase(O);

 if isempty(O)==0
        Pval = circ_rtest(O);
    else
        Pval=NaN;
 end
    
Phase_info=[meanO,meanR,Pval];


