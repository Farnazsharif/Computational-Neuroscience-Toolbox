%% Generate Place Cells Metrics

% dirData = 'Y:\Data\GrosmarkAD\';
dirData = 'A:\Data\GrosmarkAD\';

sessions = {'Achilles_10252013','Achilles_11012013','Buddy_06272013','Cicero_09012014','Cicero_09102014',...
    'Cicero_09172014','Gatsby_08282013','Gatsby_08022013'};
animal = {'Achilles','Achilles','Buddy','Cicero','Cicero','Cicero','Gatsby','Gatsby'};

for ses =1:length(sessions)%[2 5 7],[1 3 4 6 8]%6:8
    ses
    dirses = [dirData animal{ses} '\' sessions{ses}];
    cd(dirses);
    load([sessions{ses} '.mat'])
    
    %#############################################################################################################
    %Initial conditions
    
    %#1
    trN=5;speed_thr=0.04;
    
    %#2
    delta=0.1;smooth_rate=1;smooth_phase=1;
    folder= 'Metrics_xtvts';
    
    %#3
    Rate_Matrix=xtvts;
    Phase_Matrix=xtvtph;
    BinSize=(max(behav.TXVt(:,2))./max(xtvts2(:,1)));
    
    %##############################################################################################################
    
    filename=sessions{ses};
    mkdir(folder)
    TR=[];
    TR=unique(behav.TXVt(:,4));
    TR(TR==0)=[];
    
    
    if ses==2 | ses==5 | ses==7
        
        r0=[];TRi=[];TR_all=[];TR_odd=[];
        TR_odd=find(mod(TR,2)==1);
        TR_all=TR;Steps=1;% Circular
        
        r0=[1 trN];
        Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
        TRi=TR_all(1:Seg*trN);
        TRi=reshape(TRi,[trN Seg])';
        tr=TRi
        
        tiral=[];
        tic
        for j=1:(Seg)
        [tiral{1,j}]=Linearmaze_stats(filename,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,BinSize,speed_thr,delta,tr(j,:),'Cellinfo3');
        end
        [tiral{1,(Seg+1)}]=Linearmaze_stats(filename,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,BinSize,speed_thr,delta,TR_all,'Cellinfo3');
        toc
        
        cd (folder)
        save(['tiral_1'],'tiral')
        cd ..
        
        
    else
        
        r0=[];TRi=[];TR_all=[];TR_odd=[];
        TR_odd=find(mod(TR,2)==1);
        TR_all=TR(TR_odd);Steps=2;% 1 D linear
        
        r0=[1 trN];
        Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
        TRi=TR_all(1:Seg*trN);
        TRi=reshape(TRi,[trN Seg])';
        tr=TRi
        
        tiral=[];
        tic
        for j=1:(Seg)
            [tiral{1,j}]=Linearmaze_stats(filename,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,BinSize,speed_thr,delta,tr(j,:),'Cellinfo3');
        end
        [tiral{1,(Seg+1)}]=Linearmaze_stats(filename,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,BinSize,speed_thr,delta,TR_all,'Cellinfo3');
        toc
        
        cd (folder)
        save(['tiral_2'],'tiral')
        cd ..
        
        r0=[];TRi=[];TR_all=[];TR_even=[];
        TR_even=find(mod(TR,2)==0);
        TR_all=TR(TR_even);Steps=2;% 1 D linear
        
        r0=[1 trN];
        Seg=((length(TR_all)-mod(length(TR_all),trN))./trN);
        TRi=TR_all(1:Seg*trN);
        TRi=reshape(TRi,[trN Seg])';
        tr=TRi
        
        tiral=[];
        tic
        for j=1:(Seg)
            [tiral{1,j}]=Linearmaze_stats(filename,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,BinSize,speed_thr,delta,tr(j,:),'Cellinfo3');
        end
        [tiral{1,(Seg+1)}]=Linearmaze_stats(filename,Rate_Matrix,Phase_Matrix,smooth_rate,smooth_phase,BinSize,speed_thr,delta,TR_all,'Cellinfo3');
        toc
        
        cd (folder)
        save(['tiral_3'],'tiral')
        cd ..
        
    end
    
end

%%

%'Metrics_xtvts',     trials, xtvts,  delta=0.1; smooth_rate=1; smooth_phase=1; Binsize=2

%'Metrics_smooth_1',  trials, xtvts1, delta=0.1; smooth_rate=1; smooth_phase=1; Binsize=2

% xtvts2_3 binsize=3 has no metrics

%'CFC3' low and mid gamma ndx nad CFC for 10 trials

%'CFC2' CFC for 5 trials

% Cellinfo3 TXVtPh
% Cellinfo2 trash
% Cellinfo trash

%%
j=0;
% Trial_all=[];
Com_all=[];
for ses =1:length(sessions)
    ses
    dirses = [dirData animal{ses} '\' sessions{ses}];
    cd(dirses);
    load([sessions{ses} '_sessInfo.mat'])
    load([sessions{ses} '.mat'])
    filename=sessions{ses};
    
    cd Metrics_smooth_1
    if  ses == 1 | ses == 3 | ses == 4 | ses == 6 | ses == 8
        tiral=[];
        load('tiral_2.mat')
        j=j+1;
        Trial_all{1,j}=tiral;
        tiral=[];
        load('tiral_3.mat')
        j=j+1;
        Trial_all{1,j}=tiral;
        %         Com_all{1,j}=Com;
    else
        tiral=[];
        load('tiral_1.mat')
        j=j+1;
        Trial_all{1,j}=tiral;
        
    end
end
%%

%%

ses=1;
sessions  = {'Achilles_10252013_even'};
Trial_all{2,1}=[1319 1315 1317 1314 1310 1304 1305 1212 1302 1115 1117 1040 1102 1106 1113 1037 1038 1032 1033 1020 1021 1022 934 1004 1019 931 925 928 830 913 915 ...
    827 828 819 807 811 817 620 805 615 616 618 619 612 538 603 605 534 535 529 532 523 526 516 518 519 522 419 420 510 413 414 415 416 410 412 205 305 117 ...
    122 128 110 112 113 108 ];
sessions  = {'Achilles_10252013_odd'};
Trial_all{2,2}=[1319 1315 1317 1310 1313 1119 1302 1114 1117 1037 1039 1040 1102 1106 1033 1020 931 927 928 929 915 615 618 807 538 535 532 533 523 516 519 417 418 419 ...
    413 414 415 416 411 412 205 305 124 128 117 121 110 112 103 ];

ses=2;
sessions  = {'Achilles_11012013'};
Trial_all{2,3}=[1311 1228 1211 1225 1127 1209 1121 1122 1114 1117 1118 1031 1038 1112 1022 1029 1121 1122 1114 1117 1118 1119 1031 1038 1112 1022 1029 1009 1011 ...
    1006 1007 924 925 927 918 905 907 914 817 812 804 608 617 517 519 522 605 512  413 415 407 409 403 406 302 308 124 202 211 116 119 123 103 105 107 ];

ses=3;
sessions  = {'Buddy_06272013_even'};
Trial_all{2,4}=[610 611 1303 517 518 509 415 407 220 219 213 203];
sessions = {'Buddy_06272013_odd'};
Trial_all{2,5}=[1304 802 809 709 609 604 516 517 506 418 414 403 404 219 215 203];

ses=4;
sessions  = {'Cicero_09012014_even'};
Trial_all{2,6}=[1117 1207 1208 1013 915 1003 911 913 914 906 909 815 905 810 814 804 807 620 621 622 623 616 618 612 613 614 615 610 611 413 606 316 317 402 311 313 115 207 104 106];
sessions  = {'Cicero_09012014_odd'};
Trial_all{2,7}=[1208 1011 1013 915 1003 911 912 913 914 909 910 905 804 807 810 620 621 622 617 618 619 612 607 610 413 606 317 402 311 312 313 106 115];

ses=5;
sessions  = {'Cicero_09012014'};
Trial_all{2,8}=[1302 1219 1205 1206 1208 1202 1107 1110 1102 1105 1010 1011 1013 931 1005 928 929 930 926 920 922 911 916 918 906 908 910 902 903 825 828 820 822 824 ...
    817 615 803 609 611 612 607 602 308 309 404 109 108]

ses=6;
sessions  = {'Cicero_09172014_even'};
Trial_all{2,9}=[1215 1116 1118 1108 1111 1007 1106 934 1004 921 922 923 905 810 812 902 612 613 608 610 611 604 605 606 406 118 404 112 114 115 102 105 106];
sessions  = {'Cicero_09172014_odd'};
Trial_all{2,10}=[1216 1116 1118 1114 1108 1109 1007 1106 1004 1006 921 922 923 908 913 812 611 808 117 118 602 106];

ses=7;
sessions  = {'Gatsby_08282013'};
Trial_all{2,11}=[1404 724 725 1402 720 721 717 718 715 703 705 706  614 702 607 608 502 503 203];

ses=8;
sessions  = {'Gatsby_08022013_even'}
Trial_all{2,12}=[1521 1522 1502 1404 1412 1313 1314 1302 1304 1217 405 1207 1208 524 525 518 521 507 509 ];
sessions= {'Gatsby_08022013_odd'}
Trial_all{2,13}=[1521 1504 1505 1410 1313 1302 1304 1304 1214 1210 1205 1207 1208 522 517 507 509 504 302 ];

%%
save(['Trial_all_s1.mat'],'Trial_all')


