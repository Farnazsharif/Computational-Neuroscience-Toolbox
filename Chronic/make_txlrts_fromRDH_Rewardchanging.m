function txlrts = make_txlrts_fromRDH_Rewardchanging( TX,reward_time )
%%
txlrts = [TX zeros(size(TX,1), 4)];
% Trial making
    ini_trial = 1;
    txlrts(find(diff(TX(:,2))< -10)+1,5)=1;
    txlrts(1,5)=ini_trial ;
    txlrts(:,5) = cumsum(txlrts(:,5));
%    
Ndx=[];
for i=1: length(reward_time)    
[~, Ndx(i,1)]=min(abs(reward_time(i)-TX(:,1)));
end

%%
if mod(length(Ndx),2)==0
    Ndx(1:2,:)=[];
else
    Ndx(1:3,:)=[];
end
%%
Ndx=reshape(Ndx,2,length(Ndx)/2);
RX=[TX(Ndx(1,:),2), TX(Ndx(2,:),2)];

% Reward and set making
if RX(1,2)-RX(1,1)>10
    RX_ndx=find(diff( RX(:,1))< -10 | diff( RX(:,1))> 10 );
    reward_position=[RX(RX_ndx(1),1)  RX(RX_ndx(1)+1,1)]
    txlrts(Ndx(1,:),4)=1;
    set_begtrial=[ 1 ; txlrts(Ndx(1,RX_ndx),5)+1]

    for k=1:length(set_begtrial)
    
        [a,b]=find(txlrts(:,5)==set_begtrial(k));
        txlrts(a(1),6)=1;
    end
     txlrts(:,6) = cumsum(txlrts(:,6));

else
    
    reward_position=1
    if  txlrts(1,2)==0
        txlrts(1,4)=1;
    end
    txlrts(find(diff(TX(:,2))< -10)+1,4)=1;
    txlrts(1,6)=1 ;
    txlrts(:,6) = cumsum(txlrts(:,6));
end



    
     
     


     
     
     
     
     
     
     