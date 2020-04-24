

clearvars;clc;

sessions  = {'Achilles_10252013_even'}
CellG_even{1,1}=[1319 1315 1317 1314 1310 1304 1305 1212 1302 1115 1117 1040 1102 1106 1113 1037 1038 1032 1033 1020 1021 1022 934 1004 1019 931 925 928 830 913 915 ...
    827 828 819 807 811 817 620 805 615 616 618 619 612 538 603 605 534 535 529 532 523 526 516 518 519 522 419 420 510 413 414 415 416 410 412 205 305 117 ...
    122 128 110 112 113 108 ];
sessions  = {'Achilles_10252013_odd'}
CellG_odd{1,1}=[1319 1315 1317 1310 1313 1119 1302 1114 1117 1037 1039 1040 1102 1106 1033 1020 931 927 928 929 915 615 618 807 538 535 532 533 523 516 519 417 418 419 ...
    413 414 415 416 411 412 205 305 124 128 117 121 110 112 103 ];

sessions  = {'Buddy_06272013_even'}
CellG_even{3,1}=[610 611 1303 517 518 509 415 407 220 219 213 203];
sessions = {'Buddy_06272013_odd'}
CellG_odd {3,1}=[1304 802 809 709 609 604 516 517 506 418 414 403 404 219 215 203];

sessions  = {'Cicero_09012014_even'}
CellG_even{4,1}=[1117 1207 1208 1013 915 1003 911 913 914 906 909 815 905 810 814 804 807 620 621 622 623 616 618 612 613 614 615 610 611 413 606 316 317 402 311 313 115 207 104 106];
sessions  = {'Cicero_09012014_odd'}
CellG_odd{4,1}=[1208 1011 1013 915 1003 911 912 913 914 909 910 905 804 807 810 620 621 622 617 618 619 612 607 610 413 606 317 402 311 312 313 106 115];

sessions  = {'Cicero_09172014_even'};
CellG_even{6,1}=[1215 1116 1118 1108 1111 1007 1106 934 1004 921 922 923 905 810 812 902 612 613 608 610 611 604 605 606 406 118 404 112 114 115 102 105 106];
sessions  = {'Cicero_09172014_odd'};
CellG_odd{6,1}=[1216 1116 1118 1114 1108 1109 1007 1106 1004 1006 921 922 923 908 913 812 611 808 117 118 602 106];

sessions  = {'Gatsby_08022013_even'}
CellG_even{8,1}=[1521 1522 1502 1404 1412 1313 1314 1302 1304 1217 405 1207 1208 524 525 518 521 507 509 ]
sessions= {'Gatsby_08022013_odd'}
CellG_odd{8,1}=[1521 1504 1505 1410 1313 1302 1304 1304 1214 1210 1205 1207 1208 522 517 507 509 504 302 ]; 


sessions  = {'Achilles_11012013'}
CellG{2,1}=[1311 1228 1211 1225 1127 1209 1121 1122 1114 1117 1118 1031 1038 1112 1022 1029 1121 1122 1114 1117 1118 1119 1031 1038 1112 1022 1029 1009 1011 ...
    1006 1007 924 925 927 918 905 907 914 817 812 804 608 617 517 519 522 605 512  413 415 407 409 403 406 302 308 124 202 211 116 119 123 103 105 107 ];
sessions  = {'Cicero_09012014'}
CellG{5,1}=[1302 1219 1205 1206 1208 1202 1107 1110 1102 1105 1010 1011 1013 931 1005 928 929 930 926 920 922 911 916 918 906 908 910 902 903 825 828 820 822 824 ...
    817 615 803 609 611 612 607 602 308 309 404 109 108]
sessions  = {'Gatsby_08282013'}
CellG{7,1}=[1404 724 725 1402 720 721 717 718 715 703 705 706  614 702 607 608 502 503 203]


dirData = 'Y:\Data\GrosmarkAD\';
% dirData = 'A:\Data\GrosmarkAD\';

sessions  = {'Achilles_10252013','Achilles_11012013','Buddy_06272013','Cicero_09012014','Cicero_09102014',...
    'Cicero_09172014','Gatsby_08282013','Gatsby_08022013'};
animal = {'Achilles','Achilles','Buddy','Cicero','Cicero','Cicero','Gatsby','Gatsby'};
j=0;

 PPhase_odd_1=[];PPhase_odd_4=[];PPhase_odd_all=[];
 PPhase_even_1=[];PPhase_even_4=[];PPhase_even_all=[];
 PPhase_1=[];PPhase_4=[];PPhase_all=[];
 
for ses =[1 3 4 6 8]%1:length(sessions)%[2 5 7],[1 3 4 6 8]%6:8
    
    dirses = [dirData animal{ses} '\' sessions{ses}];
    cd(dirses);
    load([sessions{ses} '_sessInfo.mat'])
    load([sessions{ses} '.mat'])  
    filename=sessions{ses};
    
    %%
    cd Metrics
    if  ses == 1 | ses == 3 | ses == 4 | ses == 6 | ses == 8
     %%   
        goodCells= CellG_odd{ses, 1};
        nd=[];
        for i=1:length(goodCells)
            nd=[nd;find(G==goodCells(i))];
        end
         
        load('tiral_odd.mat') 
        meanO1_C=[];meanO4_C=[];meanO1_all=[];Linear_C=[];
        for k=1:length(nd)
            
            CellN=[];
            CellN=nd(k);
            meanO_all=[];meanO1=[];meanO4=[];Range=[];Linear=[];
            
           for j=1:size(tiral,2)
                
                meanO1=[meanO1; tiral{1,j}.Phase_info(1,1,CellN)];
                meanO4=[meanO4; tiral{1,j}.Phase_info(4,1,CellN)];
                meanO_all=[meanO_all; tiral{1,j}.Phase_info(4,1,CellN)];
                     
            end 
            
            meanO1_C=[meanO1_C meanO1];
            meanO4_C=[meanO4_C meanO4];
            meanO1_all=[meanO1_all meanO_all];
        end
        
        
       for j=1:size(tiral,2)
       Linear(j,:,:)=[tiral{1,j}.Linear(nd,:)];
       end
       size(Linear)
        
            PPhase_odd_1{ses,1}=[meanO1_C'];
            PPhase_odd_4{ses,1}=[meanO4_C'];
            PPhase_odd_all{ses,1}=[meanO1_all'];
            Linear_odd_all{ses,1}=[Linear];
        
       goodCells=[];
        goodCells= CellG_even{ses, 1};
        nd=[];
        for i=1:length(goodCells)
            nd=[nd;find(G==goodCells(i))];
        end
       

         
        load('tiral_even.mat') 
           meanO1_C=[];meanO4_C=[];meanO1_all=[];Linear_C=[];
        for k=1:length(nd)
            
            CellN=[];
            CellN=nd(k);
            meanO_all=[];meanO1=[];meanO4=[];Range=[];Linear=[];
            
            for j=1:size(tiral,2)
                
                meanO1=[meanO1; tiral{1,j}.Phase_info(1,1,CellN)];
                meanO4=[meanO4; tiral{1,j}.Phase_info(4,1,CellN)];
                meanO_all=[meanO_all; tiral{1,j}.Phase_info(4,1,CellN)];
                     
            end
            
            meanO1_C=[meanO1_C meanO1];
            meanO4_C=[meanO4_C meanO4];
            meanO1_all=[meanO1_all meanO_all];
            
        end
        
       for j=1:size(tiral,2)
       Linear(j,:,:)=[tiral{1,j}.Linear(nd,:)];
       end
        size(Linear)
        
        
            PPhase_even_1{ses,1}=[meanO1_C'];
            PPhase_even_4{ses,1}=[meanO4_C'];
            PPhase_even_all{ses,1}=[meanO1_all'];
            Linear_even_all{ses,1}=[Linear];
  
          cd ..
 
          
    else
       
        
        goodCells=[];
        goodCells= CellG{ses, 1};
        nd=[];
        for i=1:length(goodCells)
            nd=[nd;find(G==goodCells(i))];
        end
         
        load('tiral.mat') 
        for k=1:length(nd)
            
            CellN=[];
            CellN=nd(k);
            meanO_all=[];meanO1=[];meanO4=[];Range=[];Linear=[];
            
            for j=1:size(tiral,2)
                
                meanO1=[meanO1; tiral{1,j}.Phase_info(1,1,CellN)];
                meanO4=[meanO4; tiral{1,j}.Phase_info(4,1,CellN)];
                meanO_all=[meanO_all; tiral{1,j}.Phase_info(4,1,CellN)];
                Linear=[Linear;tiral{1,j}.Linear(CellN,:)];
                
            end
            
            PPhase_1{ses,1}=[meanO1'];
            PPhase_4{ses,1}=[meanO4'];
            PPhase_all{ses,1}=[meanO_all'];
            
        end
        
    end
    
    

end


%%
%Linear={SPI per spike,SPI per second,Sparsity,Coefficient,Selectivity,Olpher1,Olpher2,FR}

N=[];
M=[];
L=[];
for j=1:size(PPhase_even_1,1)
    
    if isempty(PPhase_even_1{j, 1})==0
    M=[M; PPhase_even_1{j, 1}(:,1:5)];
    M=[M; PPhase_odd_1{j, 1}(:,1:5)];
    size(M)
    
    N=[N; PPhase_even_4{j, 1}(:,1:5)];
    N=[N; PPhase_odd_4{j, 1}(:,1:5)];
    size(M)
    
    L=[L; Linear_odd_all{j, 1}(:,1:5)];
    L=[L; Linear_even_all{j, 1}(:,1:5)];
    size(L)
    end  
end
%%
clc
%Linear={SPI per spike,SPI per second,Sparsity,Coefficient,Selectivity,Olpher1,Olpher2,FR}
L=[];
C=8;
for j=1:size(Linear_odd_all,1)
         L1=[];

    for h=1:5     
    if isempty(Linear_odd_all{j, 1})==0  
    L1(:,h)=Linear_odd_all{j,1}(h,:,C);
    size(L1)
    end     
    end
            L=[L;L1];
     
    L1=[];
    for h=1:5     
    if isempty(Linear_even_all{j, 1})==0  
    L1(:,h)=Linear_even_all{j,1}(h,:,C);
    size(L1)
    end     
    end

            L=[L;L1];

    
end

%%
Q={'Olypher'};

U2=L(:,1);
U3=L(:,2);
U5=L(:,3);
U6=L(:,4);
U2(find(isnan(U2)==1))=[];
U3(find(isnan(U3)==1))=[];
U5(find(isnan(U5)==1))=[];
U6(find(isnan(U6)==1))=[];
% frist
CA1.C_1=U2;
CA1.P_1=U3;%[U3;U3];%U3;
CA1.C_2=U5;
CA1.P_2=U6;%[U6;U6];%;

U2=L(:,5);
U3=0;
U5=0;
U6=0;
CA3.C_1=U2;
CA3.P_1=U3;%[U3;U3];%U3;%
CA3.C_2=U5;
CA3.P_2=U6;%[U6;U6];%U6;

% Xch=CA1.C_1;
% nd_slop=find(Xch<mean(Xch)-1.4.*std(Xch))
% CA1.C_1(nd_slop)=[];

Boxplot_f(CA3,CA1,Q{1})
barplot_f(CA3,CA1,Q{1})
%%
figure;
Xch=CA1.C_1;
A=(mean(Xch)-1.5.*std(Xch)) 
plot(Xch,'.')
hold on
plot ([0 length(Xch) ], [A A])


%%
Q={'Last'};

U2=N(:,1);
U3=N(:,2);
U5=N(:,3);
U6=N(:,4);
U2(find(isnan(U2)==1))=[];
U3(find(isnan(U3)==1))=[];
U5(find(isnan(U5)==1))=[];
U6(find(isnan(U6)==1))=[];
% frist
CA1.C_1=U2;
CA1.P_1=U3;%[U3;U3];%U3;
CA1.C_2=U5;
CA1.P_2=U6;%[U6;U6];%;

U2=N(:,5);
U3=0;
U5=0;
U6=0;
CA3.C_1=U2;
CA3.P_1=U3;%[U3;U3];%U3;%
CA3.C_2=U5;
CA3.P_2=U6;%[U6;U6];%U6;

Boxplot_f(CA3,CA1,Q{1})
barplot_f(CA3,CA1,Q{1})






%%