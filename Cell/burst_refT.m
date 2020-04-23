function [burstIndex,refractoryT]=burst_refT(filename,acgy)
%%
% filename='FM05_1';
% acg=acg.acg1run;

 load([filename '.mat'])

for ii=1:length(G)
    
    if sum(acgy(:,ii))>100
        
        [maxacg,Imaxacg]=max(acgy(41:50,ii));
        
        edge=mean(acgy(1:10,ii));
		
        amp=maxacg-edge;
		
        if amp>0
			burstIndex(ii)=amp/maxacg;
		else
			burstIndex(ii)=amp/edge;
		end
		
		refractoryT(ii)=ACGrefractoryT(acgy(:,ii));
		
	else
	
		burstIndex(ii)=nan;
		refractoryT(ii)=nan;
	end
	
end

acg.burstIndex = burstIndex;
acg.refractoryT = refractoryT;
