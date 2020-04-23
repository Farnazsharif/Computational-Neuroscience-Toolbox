smooth=10;
% filename='DJ22_2';
filename='DJ_m8_1';
load([filename '.mat']);
%%
figure


for Celln=1 %142:length(G)%[int_celln]
subplot(1,2,1)    
plotPFSetObj_DJ(filename,Celln,smooth)
subplot(1,2,2) 
waveform_plot(filename,Celln)

% print('-depsc', [' Cell #  '  num2str(Celln) ])

end

%%
smooth=10;
filename='DJ22_2';
Celln=128;
tr_Select=[20 60];
% Fr_max=max(rr);
Fr_max=10;
nPFbin=100;


load([filename '.mat'])
xttsc=xttsc5;

xttsc=xttsc(find(xttsc(:,3)==tr_Select(1) & xttsc(:,1)==1):find(xttsc(:,3)==tr_Select(2) & xttsc(:,1)==100),:);
smoothT=smooth1D(xttsc(:,2),smooth,1);
smoothC=smooth1D(xttsc(:,Celln+4),smooth,1);
ncell=length(smoothC(1,:));
rate=smoothC./repmat(smoothT,1,ncell);
xtsr=[xttsc(:,[1 3 4]) rate]; 
rr=xtsr(:,4)./Fr_max;
matC=reshape(rr,nPFbin,length(rr)/nPFbin)';

figure
%  matC=matnorm(matC,2);
subplot(1,2,1) 
imagesc(matC.^(1))
% subplot(1,2,2) 
% waveform_plot(filename,Celln)

%%
