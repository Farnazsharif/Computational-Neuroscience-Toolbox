function G_maker(filename_C,shankN,chN)
% filename_C='DJ57_S1C';
% shankN=6;
% chN=10;
% name='uncombined';
[G_C,spk,spkinfo.waveform]=SpikeTimeWave3(filename_C,1:shankN,chN);
save('G_C.mat','G_C')

% dirinfo1=dir(name);
% [ds,~]=size(dirinfo1);
% 
% fln1=0;
% for hh=1:ds
%     if isempty(strfind(dirinfo1(hh).name,'DJ'))==1 %% detect folder name 'DJ'
%         fln1=fln1+1;
%     end
% end
% 
% cd(name)
% 
% for j=1:ds-fln1
% 
% cd(dirinfo1(j+fln1).name)
% filename=dirinfo1(j+fln1).name;
% 
% load([filename '.mat'])
% ndx=find(ismember(G_C,G)==0);
% 
% if length(G_C)~=length(G)
% G=G_C;
% G(ndx)=NaN;
% save([filename '.mat'],'-append','G')
%   
% elseif isempty(ndx)==0;
% disp( 'Warning ! Cell oreders are mismatched')
% ndx
% end
% 
% cd ..
% end
% cd ..

