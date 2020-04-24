function TXVtC=make_spike_matrix(Matrix,spkt,spkg)
%%
% spkg=spk.g;
% spkt=spk.ts;
[~,b]=size(Matrix);
TXVtC=Matrix;
if TXVtC(1,1)<0
    TXVtC=TXVtC(find(TXVtC(:,1)>0),:);
end

gsub=unique(spkg);
%%
for ii = 1:length(gsub)
    tt=selectgroup(spkt,spkg,gsub(ii));
	tt=tt(find(tt>=TXVtC(1,1)&tt<=TXVtC(end,1)));
	TXVtC(:,b+ii)=hist(tt,TXVtC(:,1))';
end
