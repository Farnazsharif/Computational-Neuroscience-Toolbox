function txlrts = make_txlrts_fromRDH_DJ( TX )

reward_position = input ('reward position?')
set_begtrial=input('starting trial of each set? ex[1 11 31] if no set, enter 1   ' )
ini_trial = 1;
ini_set = 1;

txlrts = [TX zeros(size(TX,1), 4)];
    txlrts(TX(:,2)==reward_position,4)=1;
    txlrts(find(diff(TX(:,2))< -10)+1,5)=1;
    txlrts(1,5)=ini_trial ;
    txlrts(:,5) = cumsum(txlrts(:,5));
    
    for k=1:length(set_begtrial)
        txlrts(1,6)=ini_set ;
        [a,b]=find(txlrts(:,5)==set_begtrial(k));
        txlrts(a(1),6)=1;
    end
     txlrts(:,6) = cumsum(txlrts(:,6));

