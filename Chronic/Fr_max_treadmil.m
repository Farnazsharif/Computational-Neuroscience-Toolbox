function [Fr_max,tr_Select] = Fr_max_treadmil(filename,Celln,smooth)

load([filename '.mat'])
nPFbin=100;
disp([  'initial trial = ' num2str(xttsc(1,3))   '        last trial = ' num2str(xttsc(end,3))  ]);
tr_Select=input('tr_beg ,tr_end = ');
xttsc=xttsc(find(xttsc(:,3)==tr_Select(1) & xttsc(:,1)==1):find(xttsc(:,3)==tr_Select(2) & xttsc(:,1)==100),:);
smoothT=smooth1D(xttsc(:,2),smooth,1);
smoothC=smooth1D(xttsc(:,Celln+4),smooth,1);
rate=smoothC./repmat(smoothT,1,1);
xtsr=[xttsc(:,[1 3 4]) rate]; 
rr=xtsr(:,4);
Fr_max=max(rr);
end

