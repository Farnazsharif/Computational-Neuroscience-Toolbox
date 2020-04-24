
function[fase1,amp2]=Amp_Ph_CFC(LFPtheta,signals,frecs1,frecs2,Fs)


Nptos=size(signals,1);
Nsenalesamp=size(signals,2); 
Nsenalesfase=size(LFPtheta,2);
tquitar=0.5; 

wave1=nan(length(frecs1),Nptos,Nsenalesfase);
wave2=nan(length(frecs2),Nptos,Nsenalesamp);
for isen=1:size(LFPtheta,2)
    wave1(:,:,isen) = wavelet_mod(LFPtheta(:,isen),1/Fs,false,frecs1,6);
end
for isen=1:size(signals,2)
    wave2(:,:,isen) = wavelet_mod(signals(:,isen),1/Fs,false,frecs2,6);
end
fase1=angle(wave1(:,round(tquitar*Fs)+1:end-round(tquitar*Fs),:));
amp2=abs(wave2(:,round(tquitar*Fs)+1:end-round(tquitar*Fs),:));

