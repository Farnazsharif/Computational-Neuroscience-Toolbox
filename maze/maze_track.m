
function maze_track(filename,celln,Speed_threshold,Set)
% %%
% filename='SFA5_S2_maze';
% celln=26
% %%



%%

b=4;
load([filename '.mat'])
load([filename 'Reward_info.mat'])
Speedndx=find(TXYVC(:,4) <= Speed_threshold);
TXYVC(Speedndx,:)=[];

plot(TXYVC(:,3),TXYVC(:,2),'.b')
hold on
if nargin<4
spkndx=find(TXYVC(:,celln+b)>0);
plot(TXYVC(spkndx,3),TXYVC(spkndx,2),'r.')

else
    [~,Trans_ndx]=min(abs(TXYVC(:,1)-sum(Delay_Trans)));
    if Set==1
    TXYVC=TXYVC(1:Trans_ndx-1,:);
    else
    TXYVC=TXYVC(Trans_ndx:length(TXYVC),:);
    end
   
spkndx=find(TXYVC(:,celln+b)>0);
plot(TXYVC(spkndx,3),TXYVC(spkndx,2),'r.')
end


set(gca,'YTickLabel',[])
% set(gca,'YTickLabel',[0:20:100])
% set(gca,'Ytick',[0:0.2:1])
ylim([0 1])
set(gca,'XTickLabel',[]);
% set(gca,'XTickLabel',[0:20:100]);
% set(gca,'Xtick',[0:0.2:1]);
xlim([0 1])
