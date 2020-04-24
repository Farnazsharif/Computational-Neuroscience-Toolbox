 %function [meanR,meanO]=meanphase(O);
%
%compute mean phase vector, amplitude (meanR) and mean phase (meanO).
%O are angles in radians

function [meanR,meanO,a,b]=meanphase(O);

a=sum(cos(O));
b=sum(sin(O));
meanR=sqrt(a^2+b^2)/length(O);

if a>=0&b>=0
	meanO=atan(b/a);
elseif a<=0
	meanO=atan(b/a)+pi;
else
	meanO=atan(b/a)+2*pi;
end

