function txdtsc = maketxdtsc3(TXDTS,spkt,ndex)


txdtsc=TXDTS;
if txdtsc(1,1)<0
    txdtsc=txdtsc(find(txdtsc(:,1)>0),:);
end

    spikeTime=(spkt(ndex));
    tt=spikeTime;

	tt=tt(find(tt>=txdtsc(1,1)&tt<=txdtsc(end,1)));
	txdtsc(:,5+1)=hist(tt,txdtsc(:,1))';
