function plot_PF_rs_W(Celln,smooth,Speed_threshold_maze)

% plot_PF_rs_W(16,10,0)
% Celln=10;
% smooth=10;
% Speed_noise_maze=50;
% Speed_threshold_maze=0;
%%

% [Fr_max,tr_Select] = max_finder(Celln,smooth,Speed_threshold_maze);

%%

name='uncombined';
dirinfo1=dir(name);
[ds,~]=size(dirinfo1);

fln1=0;
for hh=1:ds
    if isempty(strfind(dirinfo1(hh).name,'SF'))==1 %% detect folder name 'DJ'
        fln1=fln1+1;
    end
end

cd(name)


pn=reshape(1:(ds-fln1)*3,3,ds-fln1);
%figure('Position',[300 0 1200 1300]) %figure('Position',[left bottom width hight])
figure('Position',[300 0 1000 1000])

for j=1:ds-fln1
    
    
    cd(dirinfo1(j+fln1).name)
    filename=dirinfo1(j+fln1).name;
    load([filename '.mat'])
    
%     if exist('Empty.mat')==0
        
        if isempty(strfind(dirinfo1(j+fln1).name,'T'))==0
            
            
            if exist('Empty.mat')==0
                
            subplot(ds-fln1,3,pn(1,j))
%              (filename,Celln,smooth) 
            plot_selectT(filename,Celln,smooth,max(Fr_max),tr_Select(j,:))
            
            subplot(ds-fln1,3,pn(2,j))
            raster_plot(filename,Celln,tr_Select(j,:))
  
            
            end
            subplot(ds-fln1,3,pn(3,j))
            waveform_plot(filename,Celln)
            
            
        else
            
            
                    if exist('Empty.mat')==0
            subplot(ds-fln1,3,pn(1,j))
            plot_maze(filename,Celln,Speed_threshold_maze,smooth,max(Fr_max))
            
            subplot(ds-fln1,3,pn(2,j))
            maze_track(filename,Celln,Speed_threshold_maze,smooth)
                    end
            subplot(ds-fln1,3,pn(3,j))
            waveform_plot(filename,Celln)
            
            
            
        end
        hold on
        
%     end
    
    cd ..
end
cd ..
subplot(ds-fln1,3,2)
title(['Cell# =' num2str(Celln)] , 'FontSize', 20)


