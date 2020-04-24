    

midg_C=[];lowg_C=[];Delt=[];
a=2;r=150;
 N_ratio=[];%Mm=[];Ml=[];   
for j=1:8
    j
    Seg=size(Lowg_all{a, j},2);
    TXVt=Trials_Seg{1, j}{1,Seg};
    S=size(Trials_Seg{1, j}{1, 1},2);
    MM=nan(S,20);ML=nan(S,20);Mt=nan(S,20);
    Ndx_g_lg=[];Ndx_g_hg=[];

        for k=1:(Seg-1)
            
        Ndx_g_lg=Lowg_all{a, j}{1, k};
        Ndx_g_hg=Midg_all{a, j}{1, k};
        TRN=Trials_Seg{1, j}{1, k};
        [midg_tr,Lowg_tr]=count_Gamma(Ndx_g_lg,Ndx_g_hg,r,TXVt,TRN);
        MM(:,k)=midg_tr(:,3);%./midg_tr(:,5)
        ML(:,k)=Lowg_tr(:,3);%./midg_tr(:,5)
        Mt(:,k)=midg_tr(:,5);
        end
%         gg=sum(nansum(MM))+sum(nansum(ML));
%         midg_C=[midg_C;MM./gg];
%         lowg_C=[lowg_C;ML./gg];
        midg_C=[midg_C;MM];
        lowg_C=[lowg_C;ML];
        N_ratio=[N_ratio;MM./ML];
        Delt=[Delt;Mt];
        
end

Mm=[Mm;midg_C];Ml=[Ml;lowg_C];
%%

figure('position',[200 200 250 250])
% M=[];
hold on

M=[];M=[];
% H=Ml./max(max(Ml));H(:,9)=H(:,9)+0.02;
% N=Mm./max(max(Mm));
H=Ml./nanmean(nanmean(Ml));%H(:,9)=H(:,9)+.4;
N=Mm./nanmean(nanmean(Mm));
% M=N_ratio;
% M=midg_C;
% M=lowg_C;
% M=lowg_C./Delt;
% M=midg_C./Delt;
% M=M./Delt;
% for kk=1:20 
%    M(:,kk)= (H(:,kk)-N(:,kk))./nansum(H(:,kk)+N(:,kk));
% end
M=H;
% M=H;
Xpos=1:20;


for h=1:20
    O=[];
    O=M(:,h);
    O(find(isnan(O)==1))=[];
    O(find(isinf(O)==1))=[];
    O(O==0)=[];
    a(h)= mean(O);
    s(h) = std(O);
    if h==100
    b1=errorbar(Xpos(h),a(h),(s(h))/sqrt(10*length(O)));    
    else
    b1=errorbar(Xpos(h),a(h),(s(h))/sqrt(length(O)));
    end
    set(b1,'linewidth',2,'color','k')
    hold on
end


hold on
plot(a,'.k','markersize',20)
hold on
plot(a,':b','linewidth',2)
xlim([0 10])
