
function TXVtPh=Cell_TXVtPh(filename,CellN,Fs,thetaP)

load([filename '_sessInfo.mat'])
load([filename '.mat'])
%%
Spike_time=spk.ts(find(spk.g==CellN));
nd_maze=find(Spike_time>=sessInfo.Epochs.MazeEpoch(1,1) & Spike_time<=sessInfo.Epochs.MazeEpoch(1,2));
Spike_maze=Spike_time(nd_maze);
nd=[];
for i=1:length(Spike_maze)
    [~,nd(i)]= min(abs(behav.TXVt(:,1)-Spike_maze(i)));
end

TXVtPh=[Spike_maze behav.TXVt(nd,2:4) thetaP(round(Spike_maze*Fs))];