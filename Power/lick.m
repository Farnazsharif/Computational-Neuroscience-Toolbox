filename='FM05_1';
tstamps=TDTtimestampconv([filename '_shk1TimeStamp.txt']);
lickT=TDTtimestampconv([filename '_TRDlick.txt']);
lickT=lickT-tstamps(1);
L=zeros(length(behav.TXDTS),1);
    for i=1:length(lickT)
        [~,idx]=min(abs(lickT(i)-behav.TXDTS(:,1)));
        L(idx)=L(idx)+1;
    end
behav.TXDTSL=[behav.TXDTS L];

%%
% [restT1,runT1]=findrunrestT(behav.TXDTS(:,1),behav.TXDTS(:,3),5,1);

%% i=8
for i=1:length(restT)
a=find(behav.TXDTS(:,1)==restT(i,1));
b=find(behav.TXDTS(:,1)==restT(i,2));
restT(i,3)=sum(behav.TXDTSL(a:b,6));
end

for i=1:length(runT)
e=find(behav.TXDTS(:,1)==runT(i,1));
f=find(behav.TXDTS(:,1)==runT(i,2));
runT(i,3)=sum(behav.TXDTSL(e:f,6));
end

%%
[trial]=Trialfinder1(filename);
for i=1:length(trial)
trial(i,5)=sum(behav.TXDTSL(trial(i,1):trial(i,2),6));
end
%%
% sum(behav.runT(3))


% L1 = find(behav.TXDTS(a,1)< lickT & lickT <=behav.TXDTS(b,1));
% L2 = find(behav.txrts(c,1)< lickT & lickT <=behav.txrts(d,1));
% length(L1)
% length(L2)
% sum(txrtsl(c:d,6))












