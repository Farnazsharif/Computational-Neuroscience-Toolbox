function [cmats_sig,strength,pref_phase,mean_phase] = thetamod(LFPtheta,signals,freqs,Fs,band)
%   Detailed explanation goes here

%% initializa
t = (0:(size(LFPtheta,1)-1))'/Fs;
% Find theta epochs  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% seria mejor eliminar esta parte y 
% quitar de todos lados lo de thetaBools ya q el input ya es solo theta
Fs_spec = Fs/10;
t_spec = downsample(t,10);
[thetaSpec,theta_f] = long_wavespec(downsample(LFPtheta, 10), t_spec, 10, [1 20], 1);
%nthetaAvg = 2*round(0.25*Fs_spec/2) - 1;

smooth_sec = 2;
theta_Finds = find_inds(5, theta_f):find_inds(8, theta_f);
delta_Finds = find_inds(1, theta_f):find_inds(4, theta_f);
thetaDeltaRatio = runavg([zeros(round(smooth_sec*Fs_spec/2),1); mean(thetaSpec(theta_Finds,:),1)'./mean(thetaSpec(delta_Finds,:),1)'; ...
    zeros(round(smooth_sec*Fs_spec/2),1)], round(smooth_sec*Fs_spec/2)*2+1);
thetaDeltaInterp = [lin_interp(thetaDeltaRatio, 10); ones(9,1)*thetaDeltaRatio(end)];

TDratio_thresh = 1.2;
thetaBools = thetaDeltaInterp > TDratio_thresh;
%thetaInds = find(thetaBools);
if length(thetaBools) > length(LFPtheta)
thetaBools=thetaBools(1:length(LFPtheta));
end

% Theta parameters using Hilbert transform
theta = filtsig(LFPtheta, 1000/Fs, [2 4 12 14]/1000);
hilb = hilbert(theta);
thetaP = angle(hilb);
%thetaH = abs(hilb);

% Theta parameters using Belluscio's method
%[thetaP, thetaH] = thetaParams_fromExtrema (LFP, Fs, thetaEpochInds);

% Theta phase bins
nbins_theta = 32;
xtheta = 0:(2*pi/nbins_theta):(2*pi);
[~,thetabins] = histc(mod(thetaP, 2*pi), xtheta);

LFP_thetaRUNAvg = zeros(nbins_theta,size(LFPtheta,2));
for i=1:nbins_theta
    LFP_thetaRUNAvg(i,:) = mean(LFPtheta((thetabins == i) & thetaBools,:));
end

%% Spectrograms
fspec = freqs';
cwtScales = Fs./fspec;

signalsCWTs = zeros(length(fspec), size(signals,1), size(signals,2));
for i=1:size(signals,2)
    tmpcwt = abs(cwt(signals(:,i), cwtScales, 'cmor3-1'));
    signalsCWTs(:,:,i) = tmpcwt;
end
clear tmpcwt
ThetaSpectra = zeros(length(fspec), nbins_theta, size(signals,2));

for i=1:nbins_theta
    ThetaSpectra(:,i,:) = mean(signalsCWTs(:,(thetabins == i) & thetaBools,:), 2);
end

%% Plot theta modulation
plotsignals = 1:size(signals,2);
plotfreqs = 1:length(fspec);

% Heat map matrices
cmats_sig = zeros(length(plotfreqs), nbins_theta, length(plotsignals),2);
for j=1:length(plotsignals)
    for k=1:length(plotfreqs)
        cmats_sig(k,:,j,1) = zScore(ThetaSpectra(plotfreqs(k),:,j), mean(signalsCWTs(plotfreqs(k),thetaBools,j)), ...
            std(signalsCWTs(plotfreqs(k),thetaBools,j)));
    end
end
clims = max(maxall(abs(cmats_sig)))*[-1 1];

figure
axs = zeros(length(plotsignals),2);
imgs = zeros(length(plotsignals),2);
lines = zeros(length(plotsignals),2);
axs(1,1) = subplot(1,length(plotsignals),1);
for j=1:length(plotsignals)
    axs(j,1) = subplot(1,length(plotsignals),j);
    imgs(j,1) = imagesc([midpoints(xtheta(1:2)) 4*pi-midpoints(xtheta(1:2))]*180/pi, fspec(plotfreqs([1 end])), cmats_sig(:,[1:end 1:end],j,1), clims);
    set(gca, 'YDir', 'normal', 'XLim', [0 720], 'XTick', 0:180:720);
    hold on;
    lines(j+1,1) = plot([midpoints(xtheta) midpoints(xtheta)+2*pi]*180/pi, (cos([midpoints(xtheta) midpoints(xtheta)]) - 1)*10 + 150, 'k--');
    hold off;
    %ylabel(['Hz']);
    %title([gen_names{plotsignals(j)} ]);
    %cbar = colorbar ('SouthOutside');
end

%% cuantifications
phases=0:11.612:360; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% poner esto en función de thetabins no asi 
fc1v=abs(freqs-band(1));
fc1= find (fc1v==min(fc1v)); 
fc2v=abs(freqs-band(2));
fc2= find (fc2v==min(fc2v)); 
tetmod=cmats_sig(fc1:fc2,:,:,:);
% modulation stregth
strength=zeros(size(signals,2),1);
for i=1:size(signals,2)
    for j=1:length(tetmod(:,1,1,1))
        aux(j,i)=max(tetmod(j,:,i,1));
    end
    strength(i,1)=mean(aux(:,i));
end
% prefered phase
pref_phase=zeros(size(signals,2),1);
aux3=zeros(length(tetmod(:,1,1,1)),size(signals,2));
for i=1:size(signals,2); 
    for j=1:length(tetmod(:,1,1,1)) % n frecs selected
        aux(j,i)=max(tetmod(j,:,i,1));
        aux2(j,i)=find(tetmod(j,:,i,1)==aux(j,i));
        aux3(j,i)=phases(aux2(j,i));
    end
        pref_phase(i,1)= mean (aux3(:,i));
end
% Mean phase 
pow_phase=zeros(size(tetmod(1,:,1,1),size(signals,2))); 
for i=1:size(signals,2)
    for j=1:length(tetmod(1,:,1,1)) % para cada bin de fase
        pow_phase(j,i)=mean(tetmod(:,j,i,1)); % pow_phase is the phase distrib of power
    end
end 
for i=1:size(pow_phase,2)
    z=pow_phase(:,i)'.*exp(1i*phases);
    n=sum(pow_phase(:,i));
    zm=sum(z)./n;
    mean_phase(i)=abs(angle(zm));
end
mean_phase=(mean_phase'.*180)/pi;
end