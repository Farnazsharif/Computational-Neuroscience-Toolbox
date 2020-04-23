function txdtsc = maketxdtsc2(TXDTS,spkt,spkg,ndex)
%txdtsc = maketxdtsc2(TXDTS,spkt,spkg,CellN)
txdtsc=TXDTS;
if txdtsc(1,1)<0
    txdtsc=txdtsc(find(txdtsc(:,1)>0),:);
end

    
%     ndex=find(spkg==CellN);
    spikeTime=(spkt(ndex));
    tt=spikeTime;

	tt=tt(find(tt>=txdtsc(1,1)&tt<=txdtsc(end,1)));
	txdtsc(:,5+1)=hist(tt,txdtsc(:,1))';
