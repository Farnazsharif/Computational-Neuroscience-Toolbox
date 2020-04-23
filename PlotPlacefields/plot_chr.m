function plot_chr(Celln,smooth,Speed_threshold_maze,apparatus_order,Mouse_filename)
% ***** Notice: to plot all the variables in ACG,Frate,spkinfo matrixes,
% check the mistmatch between G and G_C
% Celln=202;
% smooth=10;
% Speed_threshold_maze=0;
% function plot_chr(Celln,smooth,Speed_threshold_maze,apparatus_order,Mouse_filename)
% ***** Notice: to plot all the variables in ACG,Frate,spkinfo matrixes,
% check the mistmatch between G and G_C

% [Fr_max,tr_Select] = max_finder(Celln,smooth,Speed_threshold_maze);
Fr_max=1;
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
    
        
  if isempty(strfind(filename,'T'))==1

        if exist('Empty.mat')==0
          
            ax5=subplot(RowN,column,pn(1,j));
            plot_maze(filename,Celln,Speed_threshold_maze,smooth,max(Fr_max));


            ax6=subplot(RowN,column,pn(2,j));
            maze_track(filename,Celln,Speed_threshold_maze)
            
            ax7=subplot(RowN,column,pn(3,j));
            XN=find(G==G_C(Celln));
            if isempty(XN)==0
            bar(acg.acg1t,acg.acg1run(:,XN),'FaceColor',[0 0.3 0.9],'EdgeColor',[0 0.3 0.9])
            xlim([-50 50])
            set(gca,'YTickLabel',[])
            set(gca,'XTickLabel',[])
            end
            ax8=subplot(RowN,column,pn(1,j+1));
            plot_maze(filename,Celln,Speed_threshold_maze,smooth,max(Fr_max),1)
            colormap(ax8,myc)
            
            ax9=subplot(RowN,column,pn(2,j+1));
            maze_track(filename,Celln,Speed_threshold_maze,1)
            
            ax10=subplot(RowN,column,pn(3,j+1));
            plot_maze(filename,Celln,Speed_threshold_maze,smooth,max(Fr_max),2)
            colormap(ax10,myc)
            
            ax11=subplot(RowN,column,pn(4,j+1));
            maze_track(filename,Celln,Speed_threshold_maze,2)    
             end
            ax12=subplot(RowN,column,pn(4,j));
            waveform_plot(filename,Celln) 
    
        else 
        
        
        if exist('Empty.mat')==0
            
            ax1=subplot(RowN,column,pn(1,j));
            
            filename2=([filename 'PlaceField_Total']);
            plotPFSetObj_chr(filename,filename2,Celln,10)
            colormap(ax1,jet(256))
            
            ax2=subplot(RowN,column,pn(2,j));
            raster_plot(filename,filename2,Celln)
            
            ax3=subplot(RowN,column,pn(3,j));
            XN=find(G==G_C(Celln));
            if isempty(XN)==0
            bar(acg.acg1t,acg.acg1run(:,XN),'FaceColor',[0 0.3 0.9],'EdgeColor',[0 0.3 0.9])
            xlim([-50 50])
            set(gca,'YTickLabel',[])
            set(gca,'XTickLabel',[])
            end
          
        end
        ax4=subplot(RowN,column,pn(4,j));
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