function plotPFSetObj_whole(filename,a,xttsc,gu1, gu2, gu3, gu4)

% set plotting window size

% filename=filename1;
% a=Cndex;
% xttsc=xttscT; 
load([filename '.mat'])

SetID=1:10; %indicate range of set here
incr=15; %position increments between reward and vaccum 


% load([filename 'PlaceField2S.mat'])

nPFbin=100;

%smooth
smoothT=smooth1D(xttsc(:,2),10,1);
smoothC=smooth1D(xttsc(:,a+4),10,1);
ncell=length(smoothC(1,:));
rate=smoothC./repmat(smoothT,1,ncell);
xtsr=[xttsc(:,[1 3 4]) rate]; 
size(xtsr);

%select sets
xtsr=xtsr(find(ismember(xtsr(:,3),SetID)),:);

%trim beginning and end to x = 1 and x = 100
bgn=find(xtsr(:,1)==1);
ending=find(xtsr(:,1)==nPFbin);
xtsr=xtsr(bgn(1):ending(end),:);

%set bgn end indexes and reward x
sets=unique(xtsr(:,3));
setyxz=nan(length(sets),4);

for ii = 1:length(sets)
  
    ndx=find(xtsr(:,3)==sets(ii));
    setyxz(ii,1)=floor(ndx(1)/100);
    setyxz(ii,2)=floor(ndx(end)/100);
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&    
    %Reward
    ndx=find(behav.txrts(:,5)==sets(ii));
    tmp2=behav.txrts(ndx,:);    
    setyxz(ii,3)=mode(tmp2(find(tmp2(:,3)==1),2)); 
    

    %Vacuum
    setyxz(ii,4)=setyxz(ii,3)+incr;
    
end
%%
TS=unique(behav.txrts(:,4:5),'rows');
maxpos=[];
for ii=1:length(TS(:,1))
    ndx=find(behav.txrts(:,4)==TS(ii,1)&behav.txrts(:,5)==TS(ii,2));
    maxpos=[maxpos;max(behav.txrts(ndx,2))];
end
maxpos=mode(maxpos);
setyxz(:,3:4)=setyxz(:,3:4)/maxpos*100;

ii=1;
% while ii < ncell becuse we have only one cell
    
    %display belt info
%      subplot('Position', [0.25, 0.96, 0.2, 0.04]);
 subplot('Position', [gu1, gu2+gu4+0.005, gu3, 0.04]);
 
    w=Beltinfo.object_end-Beltinfo.object_bgn;
    y = zeros(length(Beltinfo.object_bgn),1);
    dy = ones(length(Beltinfo.object_bgn),1).*10;

    hold on

    for jj=1:length(Beltinfo.object_bgn) 
        u=rectangle('position',[Beltinfo.object_bgn(jj) y(jj) w(jj) dy(jj)],'FaceColor',Beltinfo.objectC{jj});
    end

    set(gca,'YTickLabel',[])
%     set(gca,'Color',[0.8,0.8,0.8])
    set(gca,'XTickLabel',[])
    
    hold off
    
    %image firing rate
    rr=xtsr(:,ii+3);
    matC=reshape(rr,nPFbin,length(rr)/nPFbin)';
    matC=matnorm(matC,2);
%      subplot('Position',[0.25, 0.5, 0.2, 0.45])
  subplot('Position',[gu1, gu2, gu3, gu4])
    imagesc(matC)
    colormap jet
%     xlabel(['Cell #' num2str(G(ii))]);
    hold on
    set(gca,'XTickLabel',[])
set(gca,'XTickLabel',[40 80 120 160 200])
set(gca,'Xtick',[20 40 60 80 100])
    %plot reward
    for jj=1:length(setyxz(:,1))
        plot([setyxz(jj,3) setyxz(jj,3)],[setyxz(jj,1) setyxz(jj,2)],'w')
        plot([setyxz(jj,4) setyxz(jj,4)],[setyxz(jj,1) setyxz(jj,2)],'w')
        plot([setyxz(jj,3) setyxz(jj,4)],[setyxz(jj,1) setyxz(jj,1)],'w')
        plot([setyxz(jj,3) setyxz(jj,4)],[setyxz(jj,2) setyxz(jj,2)],'w')
    end
    
    %draw objects
    yspan=length(matC(:,1));
    for jj=1:length(Beltinfo.object_bgn)
        plot([Beltinfo.object_bgn(jj) Beltinfo.object_bgn(jj) Beltinfo.object_end(jj) Beltinfo.object_end(jj) Beltinfo.object_bgn(jj)]/Beltinfo.length*nPFbin,[yspan+1 -1 -1 yspan+1 yspan+1],Beltinfo.objectC{jj})
    end  

    hold off
% set(gcf,'color','w');
%     R=input('print? 1 for yes')
%     if R==1
%         print('-depsc',['cell' num2str(ii)]) 
%         ii = ii+1;
%     elseif R == 3
%         ii = ii - 1;
%     else
%         ii = ii + 1;
%     end

end

