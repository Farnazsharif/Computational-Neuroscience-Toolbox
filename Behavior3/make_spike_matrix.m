function TXVtC=make_spike_matrix(Matrix,spkt,spkg)
%%
% Matrix=TXVt;
% spkt=spk.ts;
% spkg=spk.g;

%%
[~,b]=size(Matrix);
gsub=unique(spkg);

TXVtC=Matrix;
if TXVtC(1,1)<0
    TXVtC=TXVtC(find(TXVtC(:,1)>0),:);
end

for ii = 1:length(gsub)
    tt=selectgroup(spkt,spkg,gsub(ii));
	tt=tt(find(tt>=TXVtC(1,1)&tt<=TXVtC(end,1)));
	TXVtC(:,b+ii)=hist(tt,TXVtC(:,1))';
end
