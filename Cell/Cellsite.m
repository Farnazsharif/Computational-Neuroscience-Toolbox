function [Cellchn]=Cellsite(ProbN,CellXYZg,probeID)
% load([filename '.mat']);
[nshank,shank_spacing,shank_contour,nsite,site_position,site_contour]=SHKlayout(probeID);

if ProbN==1;
    pn=0;
else
    pn=64;
end

 if probeID == 2
 order=[5,4,6,3,7,2,8,1];
 else  probeID == 3
 order=[10 9 8 7 6 5 4 3 2 1];
 end
     
 
  xysite=[];

for ii = 1:nshank{1,1}
        
        xoffset = (ii-1)*shank_spacing{1,1};
        
        for jj = 1:nsite{1,1}            
            xy=[site_contour{1,1}(:,1)+site_position{1,1}(jj,1)+xoffset,site_contour{1,1}(:,2)+site_position{1,1}(jj,2)] ;
            xy(:,3)=order(jj)+((nsite{1,1}*(ii-1))+pn);
            xysite=cat(1,xysite,xy);
            
        end
end
    xysite;



for celln=1:length(CellXYZg)
[~,ndx]=min(sqrt((xysite(:,1)-CellXYZg(celln,1)).^2+(xysite(:,2)-CellXYZg(celln,2)).^2));
Cellchn(celln,1)=CellXYZg(celln,4);
Cellchn(celln,2)=xysite(ndx,3);
end
Cellchn;
