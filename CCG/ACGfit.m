%function fitparameter=ACGfit(tACG,ACG,model2use,optfitrange,optplot)
%
%ex optfitrange: [-50 -5 5 50]
% model2use=1: fittype('(a*(cos(x*2*pi()/b)+1)+c)*exp(-abs(x)/d)+e*exp(-x^2/f^2)');
%	fitparameter(1)=tmpfit.a;
%	fitparameter(2)=tmpfit.b;
%	fitparameter(3)=tmpfit.c;
%	fitparameter(4)=tmpfit.d;
%	fitparameter(5)=tmpfit.e;
%	fitparameter(6)=tmpfit.f;
%	fitparameter(7)=tmpfit(tACG(ndx2));
%	fitparameter(8)=tACG(ndx3)-tACG(ndx2);
%	fitparameter(9)=goodness.rsquare;

% model2use=2: fittype('a*exp(-abs(x)/b)+c*exp(-x^2/d^2)');
%	fitparameter(1)=tmpfit.a;
%	fitparameter(2)=tmpfit.b;
%	fitparameter(3)=tmpfit.c;
%	fitparameter(4)=tmpfit.d;
%	fitparameter(5)=tmpfit(tACG(ndx2));
%	fitparameter(6)=tACG(ndx3)-tACG(ndx2);
%	fitparameter(7)=goodness.rsquare;

function fitparameter=ACGfit(tACG,ACG,model2use,optfitrange,optplot)

if nargin>3&length(optfitrange)>0
	ndx1=find(tACG==optfitrange(1));
	ndx2=find(tACG==optfitrange(2));
	ndx3=find(tACG==optfitrange(3));
	ndx4=find(tACG==optfitrange(4));
else

	ndx0=find(tACG==0);
	tmp=find(tACG>=100);
	tmp=tmp(1);
	[mm,mndx]=max(ACG(ndx0:tmp));

	ndx2=ndx0-mndx+1;
	ndx3=ndx0+mndx-1;

	ndx1=1;
	ndx4=length(tACG);

end


if model2use==1

	gg=fittype('(a*(cos(x*2*pi()/b)+1)+c)*exp(-abs(x)/d)+e*exp(-x^2/f^2)');

	gsse=1000000000000;
	for T=100:20:160
		ff=fitoptions('method','nonlin','lower',[0 T-10 0 100 0 1],'upper',[max(ACG) T+10 max(ACG) 10000 max(ACG) 100],'startpoint',[mean(ACG) T mean(ACG) 1000 mean(ACG) 50]);	
		[aa,bb]=fit([tACG(ndx1:ndx2)';tACG(ndx3:ndx4)'],[ACG(ndx1:ndx2);ACG(ndx3:ndx4)],gg,ff);
		if bb.sse<gsse
			tmpfit=aa;
			goodness=bb;
			gsse=goodness.sse;
		end
	end
		
	fitparameter(1)=tmpfit.a;
	fitparameter(2)=tmpfit.b;
	fitparameter(3)=tmpfit.c;
	fitparameter(4)=tmpfit.d;
	fitparameter(5)=tmpfit.e;
	fitparameter(6)=tmpfit.f;
	fitparameter(7)=tmpfit(tACG(ndx2));
	fitparameter(8)=tACG(ndx3)-tACG(ndx2);
	fitparameter(9)=goodness.rsquare;

end

if model2use==2

	gg=fittype('a*exp(-abs(x)/b)+c*exp(-x^2/d^2)');	
	ff=fitoptions('method','nonlin','lower',[0 100 0 1],'upper',[max(ACG) 10000 max(ACG) 100],'startpoint',[mean(ACG) 1000 mean(ACG) 50]);
	[tmpfit,goodness]=fit([tACG(ndx1:ndx2)';tACG(ndx3:ndx4)'],[ACG(ndx1:ndx2);ACG(ndx3:ndx4)],gg,ff);
	
	fitparameter(1)=tmpfit.a;
	fitparameter(2)=tmpfit.b;
	fitparameter(3)=tmpfit.c;
	fitparameter(4)=tmpfit.d;
	fitparameter(5)=tmpfit(tACG(ndx2));
	fitparameter(6)=tACG(ndx3)-tACG(ndx2);
	fitparameter(7)=goodness.rsquare;
	
end
	
if nargin>4
	bar(tACG,ACG)
	hold on
	plot(tACG,tmpfit(tACG),'r')
	hold off
end