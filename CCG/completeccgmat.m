%function Cccgmat=completeccgmat(ccgmat);
%
%complete the symmetric other half of CCG matrix

function Cccgmat=completeccgmat(ccgmat);

Cccgmat=ccgmat;

ncell=length(ccgmat(1,:,1));

nbin=length(ccgmat(:,1,1));

for ii=2:ncell
	for jj=1:ii-1
		Cccgmat(:,ii,jj)=ccgmat(nbin-(0:nbin-1),jj,ii);
	end
end

