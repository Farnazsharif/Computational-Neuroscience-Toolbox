function [Rf_ch_shank,Rf_ch_shank_pow,bigR]=Ref_ch_finder(eeg,Winrange,Rt,EEGsamplerate,ShankN,ChN,chorder)

% filename='FM01_1';
% Winrange=100;
% Rt=Rt4;
% ShankN=6;
% ChN=10;
% chorder=[1 2 3 4 5 6 7 8 9 10];
[~,b]=size(eeg);
%%
[R_n,~]=size(Rt);
incrSize = round(Winrange/50);
bigR=[];
M=[];
for chn=1:b
    eegF=eegfilt(eeg(:,chn)',EEGsamplerate,120,180); 

for i=1:R_n
    eg(i,:)=eegF(Rt(i,chn)-Winrange:Rt(i,chn)+Winrange);
    
    % find bigest power
    egp=eeg(Rt(i,chn)-Winrange:Rt(i,chn)+Winrange,chn);
    [~,~,~,do]=spectrogram(egp,Winrange,Winrange-incrSize,120:180,EEGsamplerate);
     M(i,chn)=mean(mean(do));
end
    eegT(chn,:)=mean(eg);
    
end
   [bigR(:,1),bigR(:,2)]=max(M,[],1);
%% based on the mean amplitude of the ripple
a=max(abs(eegT),[],2);

Rf_ch_shank=[];
for i=1:ShankN
v=i*ChN-(ChN-1):i*ChN;
[s1(i),Rf1(i)]= max(a(v));
cn(i)=v(Rf1(i));
end

Rf_ch_shank(:,1)=cn;
Rf_ch_shank(:,2)=Rf1;
Rf_ch_shank(:,3)=s1;
%% based on the mean power of the ripple
nd1=[];
nd2=[];
for i=1:ShankN
v=i*ChN-(ChN-1):i*ChN;
[nd0(i),nd1(i)]=max(mean(M(:,v(1):v(10))));
for ii=1:R_n
    [~,nd2(ii,i)]=max(M(ii,v(1):v(10)));
end
    [bc1(i),bc2(i)]=max(histcounts(nd2(:,i),[1:ChN+1]));
    cn(i)=v(bc2(i));
end
Rf_ch_shank_pow(:,1)=cn;
Rf_ch_shank_pow(:,2)=bc2;
Rf_ch_shank_pow(:,3)=bc1;
Rf_ch_shank_pow(:,4)=nd1;
Rf_ch_shank_pow(:,5)=nd0;
%% ploting ripple per shanks for acute
% chorder=[1 8 2 7 3 6 4 5];
nsite=ChN;
nshk_bgn=1;
nshk_end=b;
reorder=[];
for ii = 1:ShankN
    reorder = [reorder chorder+(ii-1)*nsite];
end
M=reorder(nshk_bgn:nshk_end);
W=reshape(M,ChN,ShankN);
% eegTF=eegT(reorder,:);
% 
% for i=1:8
%     subplot(1,8,i)
% imagesc(eegTF(i*8-7:i*8,:))
% hold on
% end

%%
figure('Position', [100 100 1000 500] ) 
eegR=[];
for i=1:ShankN
eegR=cat(2,eegR,eegT(W(1,i):W(10,i),:));
end
imagesc(eegR)
colormap(parula)










