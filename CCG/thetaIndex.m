%function [thetaI,thetaf,nbcount]=thetaIndex(tacg,acg)
%
%compute the "theta index", the ratio of sinusoidal vs non-sinusoidal components
%thetaI: theta Index (0 to 1)
%thetaf: theta frequency
%nbcount: number of spikes in auto-correlogram
%
%tacg: time vector of autocorrelogram
%acg: autocorrelogram


function [thetaI,thetaf,nbcount]=thetaIndex(tacg,acg)

for ii=1:length(acg(1,:))

	nbcount(ii)=sum(acg(:,ii));
	
	if nbcount(ii)>10
		pp=ACGfit(tacg,acg(:,ii),1);
		thetaI(ii)=pp(1)/(pp(1)+pp(3));
		thetaf(ii)=1/pp(2);
	else
		thetaI(ii)=nan;
		thetaf(ii)=nan;
	end
end