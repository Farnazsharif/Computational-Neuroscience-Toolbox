function [ PP,c ] = allpower2(frq1,frq2,ch1,ch2,EEGsamplerate,Ndx,eeg)

% frq1=1;frq2=250;%ch=23;
% [ PP ] = allpower2( 'FM03_2', 40,250,2,64);
% [ PP ] = allpower2( 'FM03_2', 100,250,65,128);

A=(Ndx(end,2)+2000);

freq = frq1:frq2;
scale = frq2scal(freq,EEGsamplerate);


% for jj=ch1:ch2
tic
jj=1; % for one cell
tic
PowerT3=[];
S=[];
eegch=eeg(:,jj);
S = cwt(eegch(1:(A)),scale,'morl');
PowerT3 = (envelop(S.*S))';


Power=PowerT3;
for i=1:length(Ndx)
PP(i,:)=mean(Power(Ndx(i,1):Ndx(i,2),:),1);
end
% cd Power
% save  (['pot1_', num2str(jj)],'PP')
% cd ..
toc
%% b)-2 For cancatted windows
c=[];
 ca=[]; 
 numrows=1000;
for xx=1:length(Ndx)
 
 as=Power(Ndx(xx,1):Ndx(xx,2),:);
 B = imresize(as,[numrows frq2-frq1+1]);
 c=cat(1,B,c);

end


end

% end
%%


