function TXYt=TXY(pos,mobPeriods,time_step)

    TXYt=[];
    TPT=[pos(:,1) pos(:,4) pos(:,5)];
    
    if time_step~=0;
    TXYt(:,1)=[TPT(1,1):time_step:TPT(end,1)]';
    TXYt(:,2)=interp1(TPT(:,1),TPT(:,2),TXYt(:,1));
    TXYt(:,3)=interp1(TPT(:,1),TPT(:,3),TXYt(:,1));
    else
        TXYt=TPT;
    end
    
    
    for k=1:length(mobPeriods)
        Ndx=[];
        T=[];
        T=mobPeriods{k};  
        if sum(T)~=0
        for i=1:size(T,1)
            Ndx=[Ndx;find( TXYt(:,1) >= T(i,1) & TXYt(:,1)< T(i,2))];
        end
        TXYt(Ndx,4)=k;
        end
    end