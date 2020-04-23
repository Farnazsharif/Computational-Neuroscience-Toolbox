function [T]=binT(filename,xbinNumber)
% 
% binT(filename, size(eeg,1),EEGsamplerate,nPFbin);
% clear
% clc
% close all
% cellN=1;
% filename='SFA4_S3_TRD1';
% load([filename '.mat']);
% Matrix=behav.TXVt;
% nBins=round(max(behav.TXVt(:,2))*100,1)/2;

load([filename '.mat'])
tic
T=[];
xbin=[];

%% remove negative times
Matrix=behav.TXDTS;
negndx=find(Matrix(:,1)<0);
length(Matrix);
Matrix(negndx,:)=[];
length(Matrix);
%%
ndxzero=find([1;diff(Matrix(:,2))]<-10);
for ii = 2:length(ndxzero)
    T1=[];
    %compute xbin edges
    xstep=Matrix(ndxzero(ii)-1,2)/xbinNumber;
    binedges=0:xstep:Matrix(ndxzero(ii)-1,2);
    xbin=[xbin;(1:xbinNumber)'];
    
    %select cycles from 0 to end
    tmp=Matrix(ndxzero(ii-1):ndxzero(ii)-1,:);
    
    for jj=1:length(binedges)-1
     T1=[];T2=[];
     rndx=find(tmp(:,2)>=binedges(jj)&tmp(:,2)<binedges(jj+1));% Seb Version
%      rndx=find(tmp(:,2)>=binedges(jj)&tmp(:,2)<binedges(jj+1));     
     T_bin(jj,2)=tmp(rndx(end),1);
     if jj==1
     T_bin(jj,1)=tmp(rndx(1),1);
     else
     T_bin(jj,1)=T_bin(jj-1,2); 
     end
    end
    T=cat(1,T,T_bin);
end
toc
%%

% Teeg = (1:data_length)/EEGsamplerate;
% Ndx=[];
% 
% tic
% 
% for j=1: length(T)
%   [~, Ndx(j,1)]=min(abs(T(j,1)-Teeg));
%   [~, Ndx(j,2)]=min(abs(T(j,2)-Teeg));
% end
% Ndx;
% toc

% Ndx=behav.txrts(:,1)*(spkinfo.samplerate./25);


end

