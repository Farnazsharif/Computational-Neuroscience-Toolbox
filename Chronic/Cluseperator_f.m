function Cluseperator_f(filename,shankN,chN,apparatus_order)

% Cluseperator('SFA3_S1',6,10,{'SFA3_S1_TRD1','SFA3_S1_TRD2','SFA3_S1_maze'})
% apparatus_order={'SFA3_S1_TRD1','SFA3_S1_TRD2','SFA3_S1_maze'}
% filename='SFA3_S1';
% shankN=6;
% chN=10;

[G_C,spk,waveform]=SpikeTimeWave3_new(filename,1:shankN,chN);


    for i=1:shankN;

        
        load ([filename '_clu_' num2str(i) '.mat']) %% load combined clu
        
            SZE=[];  %% check the length of each session and each shank
        for k=1:length(apparatus_order)
            cd uncombined
            cd(apparatus_order{k})
           
            a1=Loadres([apparatus_order{k} '.res.' num2str(i)]);
           
            SZE(k,1)=length(a1);
            
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
          
      g=[];   %% make seperated clu
        for jj=1:length(apparatus_order)
            g=clu.g(SIZE(jj,1):SIZE(jj,2));
            cd 'uncombined'
            cd(apparatus_order{jj})
            save ([apparatus_order{jj} '_clu_' num2str(i) '.mat'] , 'g')
            save( [apparatus_order{jj} '.mat'],'G_C')
            cd ..
            cd ..
        end
    
    end
save([filename '.mat'],'G_C')
save([filename '.mat'],'-append','apparatus_order','spk','waveform','SIZE')