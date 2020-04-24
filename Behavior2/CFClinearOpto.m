%% CFC MEC silencing linear track
clearvars;
dir = 'A:\OptoMECLEC\OML18\day2\';
cd(dir);
load([dir 'day1_Pos.mat']);
load([dir 'LasreInt.mat']);
radCh = 21; % 0 base
lmCh = 47;
pyrCh = 1;
dgCh = 19;

%% Get laser ON / OFF intervals
% onT = []; offT = [];
% for i = 1:length(Pos.ON(:,1))
%     if Pos.ON(i,3) > 0.0
%        onT = cat(1,onT,Pos.ON(i,1));
%     end
% end
% for i = 1:length(Pos.Off(:,1))
%     if Pos.Off(i,3) > 0.0
%        offT = cat(1,offT,Pos.Off(i,1));
%     end
% end
% onInt = inds_to_intvl(onT);
% offInt = inds_to_intvl(offT);
% onInt2 = []; offInt2 = [];
% for i = 1:length(onInt)
%     if onInt(i,2)-onInt(i,1) > 2
%        onInt2 = cat(1,onInt2,onInt(i,:));
%     end
% end
% for i = 1:length(offInt)
%     if offInt(i,2)-offInt(i,1) > 2
%        offInt2 = cat(1,offInt2,offInt(i,:));
%     end
% end

%% 
layerCh = [lmCh radCh dgCh];
%ints{1} = offInt2;  ints{2} = onInt2;
ints{1} = LasreInt.Off_Intv;ints{2} = LasreInt.On_Intv;

lfp = bz_GetLFP(layerCh,'noPrompts',true);

Fs = 1250;
phaserange = 5:0.5:10;
amprange = 25:5:150;
Nbins=18; binedges=linspace(-pi,pi,Nbins+1); binedges(1)=-pi-0.1; binedges(Nbins+1)=pi+0.1;

        % get cwt in intervals
        for i = 1:length(ints)
            wave1 = bz_WaveSpec(lfp,'fvector',phaserange,'intervals',ints{i},'chanID',layerCh(1));
            phase{i} = angle(wave1.data'); 
            for ch = 1:length(layerCh)
                wave2 = bz_WaveSpec(lfp,'fvector',amprange,'intervals',ints{i},'chanID',layerCh(ch));
                amp{i}(:,:,ch) = abs(wave2.data'); 
            end
        end
        
        % tort
        for i = 1:length(ints)
            tort{i} = zeros(length(phaserange),length(amprange),length(layerCh));   
            for if1=1:length(phaserange)  
                idbins=false(size(phase{i},2),Nbins);
                for ibin=1:Nbins
                    idbins(:,ibin)=and(phase{i}(if1,:)>binedges(ibin),phase{i}(if1,:)<(binedges(ibin+1)));
                end
                for if2 = 1:length(amprange) 
                    for isen = 1:length(layerCh)
                        pamp=zeros(1,Nbins);
                        for ibin = 1:Nbins
                            pamp(ibin) = mean(amp{i}(if2,idbins(:,ibin),isen));
                        end
                        pamp=pamp/sum(pamp);
                        tort{i}(if1,if2,isen)=(log(Nbins)+sum(pamp.*log(pamp)))/log(Nbins);
                     end
                end
            end
        end
        
%% plot
layerName = {'LM','rad','DG'};
        figure;
        for ch=1:length(layerCh)
            dibujar=abs(tort{1}(:,:,ch));
            subplot(2,3,ch); hold on;
            contourf(phaserange,amprange,dibujar',20,'LineColor','none');
            maxC = max([max(max(tort{1}(:,:,ch))),max(max(tort{2}(:,:,ch)))]);clim([0 maxC]);
            colorbar ('SouthOutside'); colormap jet; 
            xlim([phaserange(1) phaserange(end)]);
            ylim([amprange(1) amprange(end)]);
            title(layerName(ch));
            if ch>1
                set(gca,'YTick',[]);
            end
        end
        for ch=1:length(layerCh)
            dibujar=abs(tort{2}(:,:,ch));
            subplot(2,3,ch+3); hold on;
            contourf(phaserange,amprange,dibujar',20,'LineColor','none');
            maxC = max([max(max(tort{1}(:,:,ch))),max(max(tort{2}(:,:,ch)))]);clim([0 maxC]);
            colorbar ('SouthOutside'); colormap jet;
            xlim([phaserange(1) phaserange(end)])
            ylim([amprange(1) amprange(end)])
            if ch>1
                set(gca,'YTick',[]);
            end
        end         
        
%% 

% lfpONs = bz_GetLFP(layerCh,'noPrompts',true,'intervals',LasreInt.On_Intv);
% lfpOFFs = bz_GetLFP(layerCh,'noPrompts',true,'intervals',LasreInt.Off_Intv);
% pyrONs = bz_GetLFP(pyrCh,'noPrompts',true,'intervals',LasreInt.On_Intv);
% pyrOFFs = bz_GetLFP(pyrCh,'noPrompts',true,'intervals',LasreInt.Off_Intv);
% lfpON=[];lfpOFF=[];pyrON=[];pyrOFF=[];
% for i = 1:length(lfpONs)
%     lfpON = cat(1,lfpON,lfpONs(i).data);
%     pyrON = cat(1,pyrON,pyrONs(i).data); 
% end
% for i = 1:length(lfpOFFs)
%     lfpOFF = cat(1,lfpOFF,lfpOFFs(i).data);
%     pyrOFF = cat(1,pyrOFF,pyrOFFs(i).data); 
% end  
% 
% layerCh = [pyrCh radCh lmCh dgCh];
% lfpON =[];
% for i = 1:length(LasreInt.On_Intv)
%     temp = LoadBinary('day1.lfp','frequency',1250,'nChannels',64,'channels',layerCh+1,...
%             'start',LasreInt.On_Intv(i,1),'duration',(LasreInt.On_Intv(i,2)-LasreInt.On_Intv(i,1)));
%     lfpON = cat(1,lfpON,temp); clear temp;
% end
% lfpOFF =[];
% for i = 1:length(LasreInt.Off_Intv)
%     temp = LoadBinary('day1.lfp','frequency',1250,'nChannels',64,'channels',layerCh+1,...
%             'start',LasreInt.Off_Intv(i,1),'duration',(LasreInt.Off_Intv(i,2)-LasreInt.Off_Intv(i,1)));
%     lfpOFF = cat(1,lfpOFF,temp); clear temp;
% end
% 
%     % theta spectrum
%     params = struct ('tapers',{[3 5]},'pad',{0},'Fs',{1250},'fpass',{[1 60]},'err',{[1 0.95]},'trialave',{1});
%     [spect1,fbins1,Serr1] = mtspectrumc(double(lfpOFF(:,1)),params);
%     [spect2,fbins2,Serr2] = mtspectrumc(double(lfpON(:,1)),params);
%     figure;
%     plot_vector(spect1,fbins1,'n',Serr2,'b',2);hold on;
%     plot_vector(spect2,fbins2,'n',Serr1,'r',2);hold on;
%     set(gca,'XLim',[params.fpass(1) params.fpass(2)]); title('OFF/ ON');
% 

%%     
load([dir 'LasreInt.mat']);
lfpONs = bz_GetLFP(layerCh,'noPrompts',true,'intervals',LasreInt.On_Intv);
lfpON = [];    
for i = 1:length(lfpONs)
    lfpON = cat(1,lfpON,lfpONs(i).data);
end
lfpOFFs = bz_GetLFP(layerCh,'noPrompts',true,'intervals',LasreInt.Off_Intv);
lfpOFF = [];    
for i = 1:length(lfpONs)
    lfpOFF = cat(1,lfpOFF,lfpOFFs(i).data);
end

    % theta spectrum
    params = struct ('tapers',{[3 5]},'pad',{0},'Fs',{1250},'fpass',{[1 60]},'err',{[1 0.95]},'trialave',{1});
    [spect1,fbins1,Serr1] = mtspectrumc(double(lfpOFF(:,3)),params);
    [spect2,fbins2,Serr2] = mtspectrumc(double(lfpON(:,3)),params);
    figure;
    subplot(1,2,1);plot_vector(spect1,fbins1,'n',Serr1,'b',2);hold on;
    set(gca,'XLim',[params.fpass(1) params.fpass(2)]);title('OFF'); 
    subplot(1,2,2);plot_vector(spect2,fbins2,'n',Serr2,'r',2);hold on;
    set(gca,'XLim',[params.fpass(1) params.fpass(2)]);title('ON'); 
    

tortON = CFCtort2(double(lfpON(1:750000,3)),double(lfpON(1:750000,:)),5:0.5:10,25:5:150,1250,1);
tortOFF = CFCtort2(double(lfpOFF(1:750000,3)),double(lfpOFF(1:750000,:)),5:0.5:10,25:5:150,1250,1);    




    