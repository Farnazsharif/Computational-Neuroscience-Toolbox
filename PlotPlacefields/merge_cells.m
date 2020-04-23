clc
clear
% for TRD2 :out bouners in Path zones > 30 micron : SFA4_S3=30 [123]
% could be double; SFA4_S4=60 [40] common in both zones;  SFA4_S4=21 [33] end of the
% pzone, may be it is in cue zone

% for TRD2 :out bouners in Cue zones <-50 micron : SFA5_S4=243 [-68]
% no problem; SFA5_S4=187 [-56] no problem; SFA5_S3=100 [-60] no problem
%*******************************************************************************************
% for TRD1 :out bouners in Path zones > 30 micron : SFA4_S3=83 [39] could
% be common in both zones 
% could be double; SFA4_S3=10 [47] no problem

% for TRD1 :out bouners in Cue zones <-50 micron : SFA4_S4=90 [-122] very
% wide but no more problem #excluded!!! ,SFA5_S2=68 [-105] no problem #
% excluded
% no problem; SFA5_S3=156 [-57] low spike with 3 trials but no problem;
% SFA3_S2=54 [-68] no problem

%*******************************************************************************************
 %for shank_Avg in TRD1: SFA5_S4=43 should be deletted in Pzone, it is vey small; SFA5_S2=68  out_b 
 %and SFA5_S2=61 should be deletted in Pzone, it is vey small 
 
 %for shank_Avg in TRD2: SFA5_S3=100 [-60] no problem should be deletted in C zone
%*******************************************************************************************
% PVR
% for TRD2 :out bouners in Cue zones:  PVR_V5(76)=0.29 in SFA5_S2 = 208 [-31] excluded because of low stability

%*******************************************************************************************
% Ripple activity
% for TRD1 :out bouners in Pue zones:  Ripple to FR(28)=14.9254 depth=[-54] it is samll spike cell,
% Ripple to FR(10)=8.8 depth=[-26] it is samll spike cell, Ripple to
% FR(9)=9.96 depth=[-27] it is samll spike cell looks like noise
% Ripple to FR(155)=22.7 depth=[+8]in SFA5_S4 (it is not a cell !!!)
% excluded
%*******************************************************************************************
CellD_C=[];
CellD_P=[];
shank_Avg=[];

FN_TRD2={'SFA4_S3','SFA4_S4','SFA5_S4','SFA5_S3','SFA5_S2','SFA3_S3','SFA3_S2'};
FN_TRD1={'SFA4_S3','SFA4_S4','SFA5_S4','SFA5_S2','SFA3_S3','SFA3_S2'};
%%
FN=FN_TRD2;
for h=1:7;
    if h==1
        cd SFA4
        cd (FN{h})
    elseif h==3
        cd SFA5
        cd (FN{h})
    elseif h==6 %***6 in TRD2 5 in TDD1
        cd SFA3
        cd (FN{h})
    else
        cd (FN{h})
    end
    
    cd uncombined
    cd ([FN{h} '_TRD2'])
    
    filename_2=([FN{h} '_TRD2']);
    load([filename_2 '.mat'])
    
    CellD_C=cat(1,CellD_C,Cell_group.C);
    
    if isnan(Cell_group.P)==1
        Cell_group.P=nan(1,5);
        save([filename_2 '.mat'],'-append','Cell_group')
    end
    
    CellD_P=cat(1,CellD_P,Cell_group.P);
    
    if isnan(Cell_group.shank_Avg_P)==1
        Cell_group.shank_Avg_P=nan(1,6);
        save([filename_2 '.mat'],'-append','Cell_group')
    end
    
    shank_A=[Cell_group.shank_Avg_C; Cell_group.shank_Avg_P]';
    shank_Avg=cat(1,shank_Avg,shank_A);
    
    
    if h==2 || h==5 || h==7 %**5 7 in TRD2 4 6 in TRD1
        cd ..
        cd ..
        cd ..
        cd ..
    else
        cd ..
        cd ..
        cd ..
    end
    
end
%%
save(['CellD_C_TRD2.mat'],'CellD_C')

save(['CellD_P_TRD2.mat'],'CellD_P')

save(['shank_Avg_TRD2.mat'],'shank_Avg')



