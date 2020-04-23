function [Fr_max,tr_Select] = max_finder(Celln,smooth,Speed_threshold_maze)


%[Fr_max,tr_Select] = max_finder(17,10,50,0)
% Celln=17;
% smooth=10;
% Speed_noise_maze=50;
% Speed_threshold_maze=0;


name='uncombined';
dirinfo1=dir(name);
[ds,~]=size(dirinfo1); %%check the number of folders
fln1=0;

for hh=1:ds
    if isempty(strfind(dirinfo1(hh).name,'DJ'))==1 %% detect folder name 'DJ'
        fln1=fln1+1;
    end
end

Fr_max=[];
cd(name);

for k=1:(ds-fln1)

    cd(dirinfo1(k+fln1).name)
    filename=dirinfo1(k+fln1).name
    
    if isempty(strfind(dirinfo1(k+fln1).name,'T'))==0
        
        
        if exist('Empty.mat')==0
            
            [Fr_max(k),tr_Select(k,:)] = Fr_max_treadmil(filename,Celln,smooth);
            
        end
        
    else
        
        
        if exist('Empty.mat')==0
            
            [~,~, Fr_max(k)]=Fr_max_maze(filename,Celln,Speed_threshold_maze,smooth);
            
        end
    end
    
    cd ..
    Fr_max
    
end
cd ..

% Fr_max=max(Fr_max);
