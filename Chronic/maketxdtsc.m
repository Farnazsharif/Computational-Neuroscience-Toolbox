%function txdtsc=maketxdtsc(TXDTS,spkt,spkg)
%
%

function txdtsc=maketxdtsc(TXDTS,spkt,spkg,reocrding,filename)
%%
[~,b]=size(TXDTS);
txdtsc=TXDTS;
if txdtsc(1,1)<0
    txdtsc=txdtsc(find(txdtsc(:,1)>0),:);
end
%%
if isempty(strfind(reocrding,'c'))==1
load([filename '.mat'])   
gsub=unique(G_C);
else
gsub=unique(spkg);
end

%%
for ii = 1:length(gsub)
    tt=selectgroup(spkt,spkg,gsub(ii));
	tt=tt(find(tt>=txdtsc(1,1)&tt<=txdtsc(end,1)));
	txdtsc(:,b+ii)=hist(tt,txdtsc(:,1))';
end

