%function LayerInfo = Clayout(filename)
%
%plot shank layout with cell XY position
%plot different cell groups with the differents symbols you specify
%if Layerinfo specified, ask to draw some limits and store them in a structure (Layerinfo .name .xylimits)


function  Clayout_2(subG,Rf_ch_shank,CA,probeID ,channel_order,CellXYZg)
%%
% filename='SFA4_S3_TRD2';
% Rf_ch_shank=Rf_ch_shank_pow;
% channel_order =[10 9 8 7 6 5 4 3 2 1];
% CA=1;
% probeID=3;
% tic
% load([filename,'.mat'])
% toc
% LayerInfo = Clayout_2(PrePostG,filename,CellN,Rf_ch_shank,CA,probeID ,channel_order,symbols)

%%

% PrePostG=Cellinfo.PrePostG_T;
 Starting_Shank = 1;%9;                 %ex: 1 for first probe, 9 for second probe
% probeID = 3;%2;                        %1 for Buzsaki32 (8 sites x 4 shanks), 2 for Buzsaki64 (8 sites x 8 shanks), 3 for Buzsaki64 (10 sites x 6 shanks)
% [5 4 6 3 7 2 8 1];  %deepest to most superficial site




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% info to set %%%%%%%%%%%%%%%%%%%%%%%%%%




% subG={G(CellN)};
% subG={G(CellN1) G(CellN2) G(CellN3)};


Msize=[8,8,8,8,8];
LW=[2 2 2 2 2];%MarkerSize
%pairs for connectome here
Gpair = [];%PrePostG;

Csymbol = {'r'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%shank layout
[nshank,shank_spacing,shank_contour,nsite,site_position,site_contour]=SHKlayout(probeID);

%Compute cell position
% CellXYZg=CellPosition(G,spkinfo.waveform.average,channel_order,spkinfo.samplerate,probeID,Starting_Shank);  %if not in filename.mat already, else --> CellXYZ=spkinfo.XYZposition;
% CellXYZg=Cellinfo.CellXYZg;

%% plot shank and site position and ripple postion

if CA==1
% [~,ndx2]=sort(Rf_ch_shank(:,5)); 
M=cool(6);
    for ii = 1:nshank{1,1}
        
%         ndx3=find(ii==ndx2);
        Re=find(Rf_ch_shank(ii,6)==channel_order);
        xoffset = (ii-1)*shank_spacing{1, 1}(1,1);
        
        fill(shank_contour{1, 1}(:,1)+xoffset,shank_contour{1, 1}(:,2),[56/256 61/256 150/256],'EdgeColor','none')
        %     hold on
        %      plot(xoffset,Rf_ch_shank(ii,3)*200,'d','color',[0 102/256 0],'linewidth',5)
        hold on
        for jj =1:nsite{1,1}
            
            fill(site_contour{1, 1}(:,1)+site_position{1, 1}(jj,1)+xoffset,site_contour{1, 1}(:,2)+site_position{1, 1}(jj,2),[0.85 0.85 0],'EdgeColor','none')
            hold on
%             fill(site_contour(:,1)+site_position(Re,1)+xoffset,site_contour(:,2)+site_position(Re,2),M(ndx3,:),'EdgeColor','none')
%             hold on
            plot([mean(site_contour{1, 1}(:,1)+site_position{1, 1}(Re,1)+xoffset)-80 mean(site_contour{1, 1}(:,1)+site_position{1, 1}(Re,1)+xoffset)+80],[mean(site_contour{1, 1}(:,2)+site_position{1, 1}(Re,2)) mean(site_contour{1, 1}(:,2)+site_position{1, 1}(Re,2))],'k--','LineWidth',1)
            hold on
            
        end
    end
    
    %% plot shank and site position
else
    for ii = 1:nshank
        
        xoffset = (ii-1)*shank_spacing;
        
        fill(shank_contour(:,1)+xoffset,shank_contour(:,2),[56/250 61/250 150/250],'EdgeColor','none')
        %     hold on
        %      plot(xoffset,Rf_ch_shank(ii,3)*200,'d','color',[0 102/256 0],'linewidth',5)
        hold on
        for jj =1:nsite
            fill(site_contour(:,1)+site_position(jj,1)+xoffset,site_contour(:,2)+site_position(jj,2),[0.85 0.85 0],'EdgeColor','none')
            
        end
    end
    
end

%% plot cell position

for gg = 1:length(subG)
    
    ndxG=find(ismember(CellXYZg(:,4),subG{gg}));
%     plot(CellXYZg(ndxG,1),CellXYZg(ndxG,2),'o','MarkerSize',Msize(gg),'LineWidth',LW(gg),'color',[3/256 192/256 60/256])
    plot(CellXYZg(ndxG,1),CellXYZg(ndxG,2),'o','MarkerSize',3,'LineWidth',2,'color','r')
%     plot(CellXYZg(ndxG,1),CellXYZg(ndxG,2),'o','MarkerSize',Msize(gg),'LineWidth',LW(gg),'color','k')
    
end

%% draw local connectome
if length(Gpair)>0
    for ii = 1:length(Gpair(:,1))
        ndx1=find(CellXYZg(:,4)==Gpair(ii,1));
        ndx2=find(CellXYZg(:,4)==Gpair(ii,2));
        if length(ndx1)>0 & length(ndx2)>0
            xy=CellXYZg([ndx1 ndx2],1:2);
            plot(xy(:,1),xy(:,2),Csymbol{1})
        end
    end
end


axis image
axis off

%draw layers
% if nargout > 0
%     LayerInfo=[];
%     ndx=1;
%     RR=1;
%     while RR==1
%         'draw boundary + "return"'
%         [xx,yy]=ginput;
%         h(ndx)=plot(xx,yy,'k--','LineWidth',1);
%
%         R=input('Type boundary name + "return". "return" to skip  ');
%         if length(R)>0
%             LayerInfo(ndx).name = R;
%             LayerInfo(ndx).xylimits = [xx yy];
%             ndx = ndx + 1;
%         else
%             set(h(ndx),'Visible','off')
%         end
%
%         RR=input('Draw another? (1 --> yes, other --> no)');
%     end
%
% end


hold off

%%



