function phasehistmap(filename,cellN,PhaseFreqVector_a)
% phasehistmap('FM05_1',1633,5,9);
% cellN=1633;
% PhaseFreqVector_a=5;
% filename='FM05_1';

load([filename '.mat']); 
EEG=readmulti([filename '.lfp'], 128);
PhaseFreqVector_b=PhaseFreqVector_a;
%%

chorder=[1 8 2 7 3 6 4 5];
nsite=8;

reorder=[];
for ii = 1:16
    reorder = [reorder chorder+(ii-1)*nsite];
end
M=reorder(65:128);
W=reshape(M,8,8);
%% Seb method
tic
for kk=1:8

for ll=1:8
ch=W(kk,ll);  
    
phase=[];
eeg=EEG(:,ch);

% [phase,g]=spkphase(spk.i,spk.g,celln,spkinfo.samplerate,eeg,[5 9],spkinfo.samplerate./25);
% n=[];
% n=histc(phase,0:2*pi/18:2*pi);
% n(end)=[];
[MI,n,phase]=cellph(cellN,ch,filename,PhaseFreqVector_a,PhaseFreqVector_b);

binrange=10:20:720;
subplot(8,8,ll+(kk-1)*8)
bar(binrange,[n;n],'b')

xlim([0 720])
set(gca,'XTickLabel',[])
set(gca,'XTickLabel',[0 180 360  540 720])
% set(gca,'Xtick',[0 180 360  540 720])
% xlabel('Phase')
% ylabel('Spik#')

end   
end

set(gcf,'color','w')
dim = [.4 .7 .3 .3];
str = (['Cell#= ' num2str(cellN) ] );
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none');
toc