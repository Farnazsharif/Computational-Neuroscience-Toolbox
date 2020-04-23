%function [fmap,TimeSpent]=PF2D(txy,spkCnt,ngrid,smooth)
%
%compute the 2D place field
%txy is the matrix [t,x,y], where t is binned
%spkCnt is the number of spikes per bin of time
%ngrid is the number of bin to divide x and y 
%smooth is the smooth window size in grid bin
%
%fmap is a matrix {ngrid,ngrid} with firing rate values
%TimeSpent is a matrix(ngrid,ngrid) with occupancy time values
%f in Hz if t in s

function [fmap,TimeSpent]=PF2D(txy,spkCnt,ngrid,smooth)
%%
%[fmap,TimeSpent]=PF2D(TXYC(:,1:3),TXYC(:,4),100,10);
% txy=TXYVC(:,1:3);
% spkCnt=TXYVC(:,5);
% ngrid=100;
% smooth=10;
%%
% integrized Pos (in the range 1...nGrid
iPos = 1+floor(ngrid*txy(:,2:3)/(1+eps));

%calculate occupancy t and nspike per xbin
tbin=txy(2,1)-txy(1,1);

TimeSpent = full(sparse(iPos(:,1),iPos(:,2),tbin,ngrid,ngrid));

nspike = full(sparse(iPos(:,1),iPos(:,2),spkCnt,ngrid,ngrid));

%do the smoothing
Fullwin=(-ngrid:ngrid)';
Smoother=exp(-Fullwin.^2/(smooth/2)^2);

TimeSpent = conv2(Smoother, Smoother, TimeSpent, 'same');
nspike = conv2(Smoother, Smoother, nspike, 'same');

%Calculate freq;
fmap=nspike./(TimeSpent+eps);


