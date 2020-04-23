function Cluseperator(filename,shankN,chN)
% Cluseperator('SFA3_S1',6,10)

% filename='SFA3_S1';
% shankN=6;
% chN=10;

[G_C]=SpikeTimeWave3_new(filename,1:shankN,chN);
save([filename '.mat'],'G_C')



name='uncombined';
dirinfo1=dir(name);
[ds,~]=size(dirinfo1); %%check the number of folders

fln1=0;
for hh=1:ds
    if isempty(strfind(dirinfo1(hh).name,'SF'))==1 %% detect folder name 'DJ'
        fln1=fln1+1;
    end
end


    for i=1:shankN;
        
        load ([filename '_clu_' num2str(i) '.mat']) %% load combined clu
        
            SZE=[];  %% check the length of each session and each shank
        for k=1:(ds-fln1)

            cd uncombined
            cd(dirinfo1(k+fln1).name)
           
            a1=Loadres([dirinfo1(k+fln1).name '.res.' num2str(i)]);
           
            SZE(k,1)=length(a1)
            
            cd ..
            cd ..
            
        end
   
       SIZE=[]; %% make a matrix for each shanks clu
        for jj=1:length(SZE)
            SIZE(jj,2)=sum(SZE(1:jj));
            if jj==1
                SIZE(jj,1)=1;
            else
                SIZE(jj,1)=SIZE(jj-1,2)+1;
            end
        end
        
      g=[];%% make seperated clu
        for jj=1:length(SZE)
            g=clu.g(SIZE(jj,1):SIZE(jj,2));
            cd 'uncombined'
            cd(dirinfo1(jj+fln1).name)
            save ([dirinfo1(jj+fln1).name '_clu_' num2str(i) '.mat'] , 'g')
            save( [dirinfo1(jj+fln1).name '.mat'],'G_C')
           
  
            cd ..
            cd ..
        end
    
    end
