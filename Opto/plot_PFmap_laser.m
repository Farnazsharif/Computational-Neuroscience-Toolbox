function plot_PFmap_laser(filename, fastON)

% filename='DJ_m8_1';
% fastON=1;

%%
figure('Position',[1 2560  500 600])
%object add trial
if fastON==1
    timeThreshold = 10;
elseif fastON==0
    timeThreshold = 100000;
end
%%
load([filename '.mat'])

switch (task_version)
    case 11
        ntrials = unique(behav.txlrt(:,5));
        object_none = ntrials(2:end);
        object_add = -1;
        trial_add = -1;
    case 12
        ntrials = unique(behav.txlrt(:,5));
        object_none = ntrials(1):BeltInfo.New-2;
        object_add = -1;
        trial_add = -1;
        %object_add = BeltInfo.New: ntrials(end);
        %trial_add = BeltInfo.New;
        
end

%%
%find laser positions
laserON_pos =min (interp1(behav.txlrt(:,1),behav.txlrt(:,2),opto.laserON_OFF(:,1)));
laserOFF_pos =max (interp1(behav.txlrt(:,1),behav.txlrt(:,2),opto.laserON_OFF(:,2)));


% if   timeThreshold == 100000;
%     laserOFF_pos =(interp1(behav.txlrt(:,1),behav.txlrt(:,2),opto.laserON_OFF(:,2)));
% else 
%     laserOFF_pos =max (interp1(behav.txlrt(:,1),behav.txlrt(:,2),opto.laserON_OFF(:,2)));
% 
% end
    
%%
%belt specification
ZeroPosition = BeltInfo.Zero;
BeltLength = BeltInfo.Length;
objectX_bgn= BeltInfo.object_bgn;
objectX_end= BeltInfo.object_end;
objectC= BeltInfo.objectC;

objectX_bgn=(objectX_bgn-ZeroPosition)+((objectX_bgn-ZeroPosition)<0)*BeltLength;
objectX_end=(objectX_end-ZeroPosition)+((objectX_end-ZeroPosition)<0)*BeltLength;
%%

%fastcontroltrial and fastlasertrial
fasttrial = fastTrialseg(behav.txlrt(:,[1 2 5]),[laserON_pos laserOFF_pos],timeThreshold);
%alltrials = unique(behav.txlrt(:,5));
%%
lasertrials = findlasertrial(behav.txlrt(:,[1 2 5]),opto.laserON_OFF(:,1));
%%
fastlasertrial = intersect(fasttrial, lasertrials);
%%
fastcontroltrial = setdiff(fasttrial, lasertrials);

nfastlasertrial = length(fastlasertrial);

nfastcontroltrial = length(fastcontroltrial);

%object added trials
Objctrtrial = sum(fastcontroltrial < trial_add);
fast_add = intersect(fastcontroltrial, object_add);
fast_none = intersect(fastcontroltrial, object_none);
Objlasertrial = sum(fastlasertrial < trial_add);
%%
load([filename '_PlaceField.mat'])


%Total_pos = max(behav.txlrt(:,2));
ndx_beltMax = find(diff(behav.txlrt(:,2))<-10);
Total_pos = median(behav.txlrt(ndx_beltMax,2));
%%%% wheel increment
d = BeltLength/Total_pos;
Total_pos = median(behav.txlrt(ndx_beltMax,2))*d;

%laser/reward position
ndx_trialBGN = find(diff(behav.txlrt(:,5))==1)+1;
trialBGN_pos = behav.txlrt(ndx_trialBGN(1),2)*d;

ndx_Reward = find(behav.txlrt(:,4)==1);
reward_pos = behav.txlrt(ndx_Reward(1),2)*d;

nbin = max(xttsc(:,1));

Reward_bin = ((reward_pos-trialBGN_pos)+((reward_pos-trialBGN_pos)<0)*Total_pos)/Total_pos*nbin;
%%%%%%%laserON_bin = (laserON_pos-trialBGN_pos)/Total_pos*nbin;
%%%%%%%laserOFF_bin = (laserOFF_pos-trialBGN_pos)/Total_pos*nbin;
laserON_bin = (laserON_pos-trialBGN_pos)/Total_pos*nbin;
laserOFF_bin = (laserOFF_pos-trialBGN_pos)/Total_pos*nbin;
%%
%object position
objectX_bgn = ((objectX_bgn-trialBGN_pos)+((objectX_bgn-trialBGN_pos)<0)*Total_pos)/Total_pos*nbin;
objectX_end = ((objectX_end-trialBGN_pos)+((objectX_end-trialBGN_pos)<0)*Total_pos)/Total_pos*nbin;


if  size(object_add,2) > 1
    NobjectX_bgn = BeltInfo.Newobject_bgn;
    NobjectX_end = BeltInfo.Newobject_end;
    NobjectC= BeltInfo.NewobjectC;
    
    NobjectX_bgn=(NobjectX_bgn-ZeroPosition)+((NobjectX_bgn-ZeroPosition)<0)*BeltLength;
    NobjectX_end=(NobjectX_end-ZeroPosition)+((NobjectX_end-ZeroPosition)<0)*BeltLength;
    NobjectX_bgn = ((NobjectX_bgn-trialBGN_pos)+((NobjectX_bgn-trialBGN_pos)<0).*Total_pos)./Total_pos*nbin;
    NobjectX_end = ((NobjectX_end-trialBGN_pos)+((NobjectX_end-trialBGN_pos)<0).*Total_pos)./Total_pos*nbin;
    
end


%compute place fields
smoothT=smooth1D(xttsc(:,2),10,1);
smoothC=smooth1D(xttsc(:,5:end),10,1);
ncell=length(smoothC(1,:));
rate=smoothC./repmat(smoothT,1,ncell);
xtr=[xttsc(:,[1 3]) rate];
%%

for ii = 1:ncell
    
    ctrmat=[];
    for jj = 1:length(fast_none)
        
        ndx = find(xtr(:,2)==fast_none(jj));
        
        if length(ndx) == nbin
            
            ctrmat = [ctrmat;xtr(ndx,2+ii)'];
            
        end
        
    end
    
    
    PFmat_add=[];
    if length(fast_add)>1
        for jj = 1:length(fast_add)
            
            ndx = find(xtr(:,2)==fast_add(jj));
            
            if length(ndx) == nbin
                
                PFmat_add = [PFmat_add;xtr(ndx,2+ii)'];
                
            end
            
        end
        
    end
    
 %%   

    lasermat=[];
    for jj = 1:length(fastlasertrial)
        
        ndx = find(xtr(:,2)==fastlasertrial(jj));
        
        if length(ndx) == nbin
            
            lasermat = [lasermat;xtr(ndx,2+ii)'];
            
        end
        
    end
   
    Pvalues=[];
    for jj = 1:length(ctrmat(1,:))
        ctrdata=ctrmat(:,jj);
        laserdata=lasermat(:,jj);
        [hh,PP]=ttest2(ctrdata,laserdata);
        Pvalues=[Pvalues,PP];
    end
    
    ndx=find(Pvalues<0.05);
    signif=Pvalues(ndx);

    
 %%   
    if size(PFmat_add,1)>0
        
        ctrmat = ctrmat(sum(ctrmat,2)>0,:);
        PFmat_add = PFmat_add(sum(PFmat_add,2)>0,:);
        PFmat = [ctrmat; PFmat_add];
    else
        if  isempty(strfind(filename,'DJ_m11_1'))==0
            if (ii==34) ctrmat = ctrmat(sum(ctrmat,2)>10,:); end
        end
        if isempty(strfind(filename,'DJ_m12_1'))==0
            if (ii==3) ctrmat = ctrmat(sum(ctrmat,2)>500,:); end
        end
        ctrmat = ctrmat(sum(ctrmat,2)>0,:);
        PFmat = ctrmat;
    end
    %%
    subplot(2,2,1);imagesc(PFmat)
    hold on
%      plot([laserON_bin laserON_bin],[0 length(fastcontroltrial)],'y--')
%      plot([laserOFF_bin laserOFF_bin],[0 length(fastcontroltrial)],'y--')

    for jj = 1:length(objectX_bgn)
        plot([objectX_bgn(jj) objectX_bgn(jj)],[0 length(fastcontroltrial)],objectC{jj})
        plot([objectX_end(jj) objectX_end(jj)],[0 length(fastcontroltrial)],objectC{jj})
    end
    
    if size(PFmat_add>0)
        hold on
        plot([NobjectX_bgn(1)-1 NobjectX_bgn(1)-1],[Objctrtrial-1 length(fastcontroltrial)],NobjectC{1})
        plot([NobjectX_end(1)   NobjectX_end(1)],  [Objctrtrial-1 length(fastcontroltrial)],NobjectC{1})
        plot([0 nbin],[Objctrtrial-1 Objctrtrial-1],'w--')
        Pvalues=[];
        for jj = 1:nbin
            [hh,pp]=ttest2(ctrmat(:,jj),PFmat_add(:,jj));
            Pvalues=[Pvalues;pp];
        end
    end
    title(G(ii));
%     hold off
    %% laser
    
    subplot(2,2,2);imagesc(lasermat)
    hold on
    plot([laserON_bin laserON_bin],[0 length(fastlasertrial)],'y--')
    plot([laserOFF_bin laserOFF_bin],[0 length(fastlasertrial)],'y--')
    if Objlasertrial>0
        plot([0 nbin],[Objlasertrial Objlasertrial],'y--')
    end
%     hold off
%     hold on 
 %%  plot mean and std of laser 
    x= 1:nbin;
    y = mean(lasermat);
    xi= [x x(end:-1:1) x(1) x(2)];
    yi= y-std(lasermat,1)/sqrt(size(lasermat,1));
    yj= y+std(lasermat,1)/sqrt(size(lasermat,1));
    yij= [yi, yj(end:-1:1) yi(1) yi(2)];
    
    subplot(2,2,4);plot(1:nbin,mean(lasermat),'r')
    hold on
    plot(xi, yij, 'k') 
   
   
    plot([laserON_bin laserON_bin],[0 max([mean(ctrmat) mean(lasermat)])*1.1],'g--')
    plot([laserOFF_bin laserOFF_bin],[0 max([mean(ctrmat) mean(lasermat)])*1.1],'g--')
    hold off
 %%  plot mean and std of contol 
    x= 1:nbin;
    y = mean(ctrmat);
    xi= [x x(end:-1:1) x(1) x(2)];
    yi= y-std(ctrmat,1)/sqrt(size(ctrmat,1));
    yj= y+std(ctrmat,1)/sqrt(size(ctrmat,1));
    yij= [yi, yj(end:-1:1) yi(1) yi(2)];
    
    
    if Objctrtrial>0
%         subplot(3,1,3);plot(1:nbin,mean(ctrmat(1:Objctrtrial,:)),'b',1:nbin,mean(ctrmat(Objctrtrial+1:end,:)),'b--')
        subplot(2,2,3);plot(1:nbin,mean(ctrmat),'b',1:nbin,mean(PFmat_add),'r')
        hold on
%         plot(1:nbin,mean(lasermat(1:Objlasertrial,:)),'r',1:nbin,mean(lasermat(Objlasertrial+1:end,:)),'r--')
    else
        subplot(2,2,3);plot(1:nbin,mean(ctrmat),'b')
        
%         subplot(2,2,4);plot(1:nbin,mean(ctrmat),'b',1:nbin,mean(lasermat),'r')
    end
    hold on
     
    plot(xi, yij) 
    



%%    
    if size(PFmat_add,1)>0
        ymax = max([max(mean(ctrmat)), max(mean(PFmat_add))]);
        ymin = min([min(mean(ctrmat)), min(mean(PFmat_add))]);
        
        
        
    else
        ymax= max(yij); %ymax = max(mean(ctrmat));
        ymin = min(yij); %ymin = min(mean(ctrmat))*0.8;
    end
    
    for jj = 1:length(objectX_bgn)
        plot([objectX_bgn(jj) objectX_bgn(jj)],[0 ymax*1.1],objectC{jj})
        plot([objectX_end(jj) objectX_end(jj)],[0 ymax*1.1],objectC{jj})
    end
    fill([Reward_bin 100 100 Reward_bin],[ymax*1.05 ymax*1.05 ymax*1.07 ymax*1.07],'b')
    if isnan(max(mean(PFmat)))==0 && max(mean(PFmat))>0
        ylim([ymin*0.8 ymax*1.1])
    end
    
    if size(PFmat_add)>0
        plot([NobjectX_bgn(1) NobjectX_bgn(1)],[0 ymax*1.1],'k')
        plot([NobjectX_end(1) NobjectX_end(1)],[0 ymax*1.1],'k')
        
        ndx=find(Pvalues(:,1)<0.05);
        plot(ndx,ones(size(ndx))*ymax*1.05,'r*')
        title (['blue: none, red:'  BeltInfo.Newobject_ID])
    end
    
  %%  
    hold off
    
    R=input('print? 1 for yes')
    
    if R==1
        cd ..
        cd ('figure')
        if fastON ==0
            if task_version ==12
                print('-depsc',['o_' filename(4:end) '_' num2str(G(ii))])
            else
                print('-depsc',[filename(4:end) '_' num2str(G(ii))])
            end
        elseif fastON==1
            print('-depsc',['f' filename(4:end) '_' num2str(G(ii))])
        end
        cd ..
        cd (filename)
    end
    
end