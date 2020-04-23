function [li,Lt,lc]= licknumber(filename,xbinNumber)

% filename='FM05_1';
% xbinNumber=100;

[xtsl,txrtsl,setyxz]=lickmap(filename,xbinNumber);
matC=reshape(xtsl(:,4),xbinNumber,length(xtsl(:,4))/xbinNumber)';
Rewx=round(setyxz);

[trial]=Trialfinder1(filename);
tr1=find(trial(:,4)==1, 1, 'last')
tr2=find(trial(:,4)==2, 1,'last')

li=[];
i=1;
a1=sum(matC(1:tr1,Rewx(i,3):Rewx(i,4)),2);
i=2;
a2=sum(matC(tr1+1:tr2,Rewx(i,3):Rewx(i,4)),2);
i=3;
a3=sum(matC(tr2+1:length(trial),Rewx(i,3):Rewx(i,4)),2);

li=cat(1,a1,a2,a3);

Lt=sum(matC,2)+0.00001;

lc=(li./Lt)*100;