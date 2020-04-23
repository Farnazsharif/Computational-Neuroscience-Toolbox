%function refractoryT=ACGrefractoryT(acg)
%
%acg bins should be in ms.

function refractoryT=ACGrefractoryT(acg)

startpoint=ceil(length(acg)/2)+1;

halfacg=acg(startpoint:end);

Shalfacg=smooth1D(halfacg,6);

ndx=LocalMaxima(Shalfacg);

firstmode=ndx(1);

Dhalfacg=diff(halfacg(1:firstmode));

tmp=find(Dhalfacg>std(Dhalfacg));

if length(tmp)>0
	refractoryT=tmp(1)+1;
else
	refractoryT=nan;
end