
function [Rt]=Ripple_t_bz(rippleT2,sample_win,eegsamplingrate,eeg)

% filename='FM05_2';
% rippleT2=rippleT2(1);
% sample_win=60;
%%
% load([filename '.mat']);
% eeg= readmulti([filename '.lfp'], channelN);
[~,e]=size(eeg);

%%
% Trange=[-0.05 0.05]+rippleT2(j);
% eegndx = round(Trange(1)*eegsamplingrate):round(Trange(2)*eegsamplingrate);
%%
teeg = (1:length(eeg))/eegsamplingrate;
% [~,Trange]=min(abs(rippleT2-teeg))
Trange=round(rippleT2*eegsamplingrate);
eegndx = (Trange-sample_win:Trange+sample_win);
%%
eegS=[];
for i=1:e
    eegS=eeg(eegndx,i);
    [B,A]=butter(4,[120 180]/eegsamplingrate*2);
    eegF(:,i)=filtfilt(B,A,eegS);
end

%%
I=[];
Ienv=[];
V=[];
Venv=[];
Rt=[];

for j=1:e
    
    [maxI,maxV]=LocalMaxima(eegF(:,j));
    [I,V]=LocalMaxima(maxV);
    if isempty (I)==1;
        Ienv(j)=Ienv(j-1);
    elseif length(I)> 1
        
        [~,Is]=max(V);
        Ienv(j)=maxI(I(Is));
        Venv(j)= V(Is);
        
    else
        Ienv(j)=maxI(I);
        Venv(j)=V;
    end
    
end

Rt=eegndx(Ienv);


%%
% for ii=1:64
%     subplot(8,8,ii)
%     plot(eegF(:,ii))
%     hold on
%     plot(Ienv(ii),Venv(ii),'r*')
% end
%
% print('-depsc',['R2_' num2str(E) ])
%
%

end

