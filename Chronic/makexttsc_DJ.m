function xttsc=makexttsc_DJ(TXDTS,spk,xbinNumber,SpeedThreshold,G_C)

d=1; %cm for 1 x increment 

txdtsc=maketxdtsc_DJ(TXDTS,spk.ts,spk.g,G_C);
[~,b]=size(TXDTS);
ncell=length(txdtsc(1,:))-b;
binT=min(diff(txdtsc(:,1)));
speed=x2v(txdtsc(:,3)*d,txdtsc(:,1),1);

%find positions 0
ndxzero=find([1;diff(txdtsc(:,2))]<-10);

%for each belt rotation, compute rate in each space bin
nspike=[];
TimeSpent=[];
xbin=[];
trialnb=[];
setnb=[];
for ii = 2:length(ndxzero)
    
    %compute xbin edges
    xstep=txdtsc(ndxzero(ii)-1,2)/xbinNumber;
    binedges=0:xstep:txdtsc(ndxzero(ii)-1,2);
    xbin=[xbin;(1:xbinNumber)'];
    
    %select rotation/speed     
    tmpspeed=speed(ndxzero(ii-1):ndxzero(ii)-1);
    tmp=txdtsc(ndxzero(ii-1):ndxzero(ii)-1,:);
    tmpS=tmp(tmpspeed>SpeedThreshold,:);
    
    for jj=1:length(binedges)-1
        
        %spike count and time spent
        ndx=find(tmpS(:,2)>binedges(jj)&tmpS(:,2)<=binedges(jj+1));
        
        if ~isempty(ndx)
            
            nspike=[nspike;sum(tmpS(ndx,6:end),1)];

            TimeSpent=[TimeSpent;length(ndx)*binT];

        else

            nspike=[nspike;zeros(1,ncell)];

            TimeSpent=[TimeSpent;0.000001];

        end
        
        %trial # and set#
        ndx=find(tmp(:,2)>binedges(jj)&tmp(:,2)<=binedges(jj+1));
        if ~isempty(ndx)
            
            trialnb=[trialnb;tmp(ndx(end),4)];

            setnb=[setnb;tmp(ndx(end),5)];

        else

            trialnb=[trialnb;trialnb(end)];
            setnb=[setnb;setnb(end)];

        end

    end
    
end

xttsc=[xbin TimeSpent trialnb setnb nspike];

