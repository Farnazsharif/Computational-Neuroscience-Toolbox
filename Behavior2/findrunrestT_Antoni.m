%function [restT,runT]=findrunrestT(T,X,speedthreshold,durationthreshold)
%
%speed threshold in cm/s
%duration threshold in s
%find rest (<speed threshold) and run (>speed threshold) periods

function [restT,runT]=findrunrestT_Antoni(T,v,speedthreshold,durationthreshold)


% T=M(:,1);
% v=M(:,3);
% speedthreshold=0.03;

%%
yy=v>speedthreshold;

dyy=diff(yy);

%runT
runstart=find(dyy==1)+1;
runstop=find(dyy==-1);
if isempty(runstart)==1
    runT=[T(1,1) T(end,1)];
    restT=[];
else
    if runstart(1)>runstop(1)&&yy(1)==1
        runstart=[1;runstart];
    end
    if runstop(end)<runstart(end)&&yy(end)==1
        runstop=[runstop;length(yy)];
    end
    
    %%
    runT=[T(runstart) T(runstop)];
    
    %restT
    reststart=find(dyy==-1)+1;
    reststop=find(dyy==1);
    if reststart(1)>reststop(1)&&yy(1)==0
        reststart=[1;reststart];
    end
    if reststop(end)<reststart(end)&&yy(end)==0
        reststop=[reststop;length(yy)];
    end
    
    restT=[T(reststart) T(reststop)];
    
    %duration threshold
    runT=runT(find((runT(:,2)-runT(:,1))>durationthreshold),:);
    restT=restT(find((restT(:,2)-restT(:,1))>durationthreshold),:);
end
%%%%%%%%%%for testing%%%%%%%%%%%%%%%%%
%{
plot(T,v,'k')
hold on
for ii=1:length(runT(:,1))
    plot(runT(ii,:),[20 20],'r')
end
for ii=1:length(restT(:,1))
    plot(restT(ii,:),[15 15],'b')
end
hold off
%}