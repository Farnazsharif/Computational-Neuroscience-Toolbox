
    function plot_maze(filename,celln,Speed_threshold,smooth,Fr_max,Set)
    %%  plot_maze(filename,1,0,3,1,2)
%     filename='SFA5_S1_maze';
%     Speed_threshold=0 ;
%     celln=33;
%     smooth=3;
%     nargin=5;
%     Set=2;
%% 
    b=4;
    load([filename '.mat'])
    load([filename 'Reward_info.mat'])
    
    Speedndx=find(TXYVC(:,4) <= Speed_threshold);
    TXYVC(Speedndx,:)=[];
    if nargin<6
            [fmap,TimeSpent]=PF2D(TXYVC(:,1:3),TXYVC(:,celln+b),100,smooth);
             % imagesc((fmap./Fr_max).^(1)) 
    % selecting the Set to plot
    else
    [~,TimeSpent]=PF2D(TXYVC(:,1:3),TXYVC(:,celln+b),100,smooth);
    [~,Trans_ndx]=min(abs(TXYVC(:,1)-sum(Delay_Trans)));
    if Set==1
    TXYVC=TXYVC(1:Trans_ndx-1,:);
    else
    TXYVC=TXYVC(Trans_ndx:length(TXYVC),:);
    end
    [fmap,~]=PF2D(TXYVC(:,1:3),TXYVC(:,celln+b),100,smooth);
   
    end

    %% Scalling and smoothing the map
    
    factor=10;
    AA=sort(reshape(fmap,1,100*100));
    if max(AA)/10>prctile(AA,[97])
    fr_thr=prctile(AA,[97]);
    else
    fr_thr=prctile(AA,[94]);
    end
    myc=jet(256);
    myc(1, :)=[1,1,1];
    TimeSpent_thr=1/(mean((TXYVC(:,4)))*7);
%     TimeSpent_thr=0.077;
    fmap_total=fmap.*(TimeSpent>TimeSpent_thr)+((TimeSpent>TimeSpent_thr)-1);
    imagesc(fmap_total)
    colormap(myc);
    C_step=length(0:fr_thr)/length(myc);
    caxis([-length(0:fr_thr+C_step)/length(myc)  fr_thr])

    axis xy 
    title(['fr95= ' num2str(prctile(AA,[95])) '   ' 'fr/thr= ' num2str(fr_thr) ])
    
%%
% m = size(jet(256),1);
% cmin = min(fmap_total(:));
% cmax = max(fmap_total(:));
% C1=min(m,round((m-1)*(fmap_total-cmin)/(cmax-cmin))+1); 

    %% plotting the reward positions
    if nargin>5
    hold on
   
    [~,by]=find(fmap_total>-1);
    r=15;
    nsegments=100;
    ndx=find(Rew_s(:,3)==Set);
    x1=Rew_s(ndx(1),1)*100;
    y1=Rew_s(ndx(1),2)*100;
    x2=Rew_s(ndx(2),1)*100;
    y2=Rew_s(ndx(2),2)*100;
 

    if x1==0 && y1==0
    x1=x1+by(1)/2;
    y1=y1+by(1)/2;
    x2=x2-by(1)/2;
    y2=y2-by(1)/2;
    s1=0;s2=pi/2;s3=pi;s4=3*pi/2;
    elseif x1==100 && y1==100
    x2=x2+by(1)/2;
    y2=y2+by(1)/2;
    x1=x1-by(1)/2;
    y1=y1-by(1)/2;
    s3=0;s4=pi/2;s1=pi;s2=3*pi/2;
    elseif x1==0 && y1==100 
    x2=x2-by(1)/2;
    y2=y2+by(1)/2;
    x1=x1+by(1)/2;
    y1=y1-by(1)/2;  
    s1=3*pi/2;s2=2*pi;s3=pi/2;s4=pi;
    else
    x1=x1-by(1)/2;
    y1=y1+by(1)/2;
    x2=x2+by(1)/2;
    y2=y2-by(1)/2;
    s3=3*pi/2;s4=2*pi;s1=pi/2;s2=pi;
    end
   
    
    th = s1:2*pi/nsegments:s2;
    xunit = r * cos(th) + x1;
    yunit = r * sin(th) + y1;
    plot(xunit, yunit,'k','linewidth',3);
    hold on
    th = s3:2*pi/nsegments:s4;
    xunit = r * cos(th) + x2;
    yunit = r * sin(th) + y2;
    plot(xunit, yunit,'k','linewidth',3);
%     hold off
    end
    %%

    % subplot(1,2,2)
    % plot(TXYVC(:,3),TXYVC(:,2))
    % hold on
    % 
    % spkndx=find(TXYVC(:,celln+b)>0);
    % plot(TXYVC(spkndx,3),TXYVC(spkndx,2),'r.')
    % 
    % set(gca,'YTickLabel',[])
    % set(gca,'YTickLabel',[0:10:100])
    % set(gca,'Ytick',[0:0.1:1])
    % 
    % set(gca,'XTickLabel',[])
    % set(gca,'XTickLabel',[0:10:100])
    % set(gca,'Xtick',[0:0.1:1])

    % hold off
