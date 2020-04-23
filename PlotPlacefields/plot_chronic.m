function plot_chronic(Celln,smooth,Speed_threshold_maze,apparatus_order,Mouse_filename)
% ***** Notice: to plot all the variables in ACG,Frate,spkinfo matrixes,
% check the mistmatch between G and G_C
% 
% Celln=115;
% smooth=3;
% Speed_threshold_maze=0;
% apparatus_order={'SFA5_S1_TRD1','SFA5_S1_TRD2','SFA5_S1_maze'};
% Mouse_filename=apparatus_order;
% [Fr_max,tr_Select] = max_finder(Celln,smooth,Speed_threshold_maze);
% apparatus_order={'SFA5_S1_maze','SFA5_S1_TRD1','SFA5_S1_TRD2'};

RowN=length(apparatus_order)+1;

name='uncombined';
cd(name)

column=4;
pn=reshape(1:RowN*column,column,RowN);

figure('Position',[300 0 1000 1000])

 for j=1:length(apparatus_order)

    cd(apparatus_order{j})
    filename=apparatus_order{j};
    load([filename '.mat'])
    
        
  if isempty(strfind(filename,'m'))==0

        if exist('Empty.mat')==0
%               myc=jet(256);
%               myc(1, :)=[1,1,1];
            subplot(RowN,column,pn(1,j));
            plot_maze(filename,Celln,Speed_threshold_maze,smooth,max(1));
%           colormap(ax1,myc);

            subplot(RowN,column,pn(2,j));
            maze_track(filename,Celln,Speed_threshold_maze)
            
            subplot(RowN,column,pn(3,j));
            XN=find(G==G_C(Celln));
            if isempty(XN)==0
            bar(acg.acg1t,acg.acg1run(:,XN),'FaceColor',[0 0.3 0.9],'EdgeColor',[0 0.3 0.9])
            xlim([-50 50])
            set(gca,'YTickLabel',[])
            set(gca,'XTickLabel',[])
            end
            subplot(RowN,column,pn(1,j+1));
            plot_maze(filename,Celln,Speed_threshold_maze,smooth,max(1),1)
                        
            subplot(RowN,column,pn(2,j+1));
            maze_track(filename,Celln,Speed_threshold_maze,1)
            
            subplot(RowN,column,pn(3,j+1));
            plot_maze(filename,Celln,Speed_threshold_maze,smooth,max(1),2)
                        
            subplot(RowN,column,pn(4,j+1));
            maze_track(filename,Celln,Speed_threshold_maze,2)    
            end
            subplot(RowN,column,pn(4,j));
            waveform_plot(filename,Celln) 
    
        else 
        
        
        if exist('Empty.mat')==0
            if j==3
                j=4;
            end
            ax9=subplot(RowN,column,pn(1,j));
            filename2=([filename 'PlaceField_Total']);
            plotPFSetObj_chr(filename,filename2,Celln,10)
            colormap(ax9,jet(256))
            
            subplot(RowN,column,pn(2,j));
            raster_plot(filename,filename2,Celln)
            
            subplot(RowN,column,pn(3,j));
            XN=find(G==G_C(Celln));
            if isempty(XN)==0
            bar(acg.acg1t,acg.acg1run(:,XN),'FaceColor',[0 0.3 0.9],'EdgeColor',[0 0.3 0.9])
            xlim([-50 50])
            set(gca,'YTickLabel',[])
            set(gca,'XTickLabel',[])
            end
          
        end
        subplot(RowN,column,pn(4,j));
        waveform_plot(filename,Celln)

        
      
    end
     hold on
    
    
    
    cd ..
 end
cd ..

subplot(RowN,column,pn(4,1))
title(['Cell# =' num2str(Celln)] , 'FontSize', 20)

    cd ([Mouse_filename '_Pics'])
    print('-djpeg',['Cell#' num2str(Celln)])
    cd ..
   