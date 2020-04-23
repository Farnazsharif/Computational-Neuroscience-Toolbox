


function CellRe(ProbN,filename,Rt)
%%

% ProbN=2;
% filename='FM05_1';
load([filename '.mat']);

Cellchn=Cellinfo.Cellchn_T;
% [208 213]%cellN
for i=1%: length(Cellchn)
    i
    close all
%     celln=Cellchn(i,1);
    ch=Cellchn(i,2);
    cellN=G(i);
    Cellripple(filename,ch,cellN,Rt )

print('-depsc',['RepCell' num2str(cellN)])
end


