

clearvars;clc;
dirData = 'Y:\OptoMECLEC\';
animal  = {'OML18','OML18','OML18','OML19','OML19'};
sessions  = {'day1','day2','day4','day2','day3'};

for ses=3:5
    ses
    dirses = [dirData animal{ses} '\'  sessions{ses}];
    cd(dirses);
    filename=[animal{ses} sessions{ses}];
    load([filename '.mat' ])
    
    TR=unique(behav.txvtl(:,4));
    TR(TR==0)=[];
    r0=[];TRi=[];TR_all=[];TR_even=[];
    TR_even=find(mod(TR,2)==0);
    TR_On=TR(TR_even);
    
    TR=[];
    TR=unique(behav.txvtl(:,4));
    TR(TR==0)=[];
    r0=[];TRi=[];TR_all=[];TR_Odd=[];
    TR_Off=find(mod(TR,2)==1);
    TR_Off=TR(TR_Off);
    
    %%
    Pos=[];
    Nd_ON=[];
    for i=1:length(TR_On)
        Nd_ON=[Nd_ON;find(behav.txvtl(:,4)==TR_On(i))];
    end
    Pos.ON=behav.txvtl(Nd_ON,:);
    
    Nd_Off=[];
    for i=1:length(TR_Off)
        Nd_Off=[Nd_Off;find(behav.txvtl(:,4)==TR_Off(i))];
    end
    Pos.Off=behav.txvtl(Nd_Off,:);
    save([sessions{ses} '_Pos'],'Pos')
    
end

