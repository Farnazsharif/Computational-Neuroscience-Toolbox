% function [ MI , MIh]= MIindex( eeg,EEGsamplerate)
% phase az manfi sort shode

function [MI,MI2]= FarnazMI( eeg,EEGsamplerate)
% phase az manfi sort shode
tic
Ap=[];
bAm=[];
MIh=[];
MI=[];
theta=0;

 for j=2:2:50;  
     
 theta=theta+1;
%  freqband1=[j j+4];
%  [B1,A1]=butter(4,freqband1/EEGsamplerate.*2);
%  Feeg1=filtfilt(B1,A1,eeg);
%  hh1=hilbert(Feeg1);
%  phase1=angle(hilbert(Feeg1));%atan(imag(hh1)./real(hh1))+(Feeg1>0)*pi()-pi();

Feeg1=eegfilt(eeg,EEGsamplerate,j,j+4); % this is just filtering 
 
phase1=angle(hilbert( Feeg1)); % this is getting the phase time series
 

 
gamma=0;
for k=10:5:200;
    
% freqband2=[k k+10];
% [B2,A2]=butter(4,freqband2/EEGsamplerate.*2);
Feeg2=eegfilt(eeg,EEGsamplerate,k,k+10);
hh2=hilbert(Feeg2);
Am2=abs(hh2);
[sphase Ig]=sort(phase1);
sAm=Am2(Ig);
binN=18;
bin=(max(phase1)-min(phase1))/binN;

for i=1:binN
    y=sphase<min(phase1)+bin*i &  sphase>min(phase1)+bin*(i-1);
    n=find(y);
    bAm(i)=mean(sAm(n));
end
gamma=gamma+1;
Ap=bAm./sum(bAm);
MIh(theta,gamma)=(max(Ap)-min(Ap))./max(Ap);
hp= Ap.*log(Ap);
MI(theta,gamma)=1-sum(hp)./(log(max(Ap)).*sum(Ap));
MI2(theta,gamma)=1+sum(hp)./(log(binN));
end
 end
toc
end

% for ll=1:25
% yy(ll,:) = smooth(MI2(ll,:),10,'loess');
% end
% figure(2)
% contourf(yy',30,'lines','none')


