


L=[];
L=LoadBinary('analogin.dat','nChannels',1,'Frequency',30000);
tic
figure;
plot(L(1:100:end),'.')

%%
L2=L;
L2(L2>0)=[];
figure;
plot(L2(1:100:end),'.')
%%
Laser=bz_getAnalogPulses('data',L,'manualThr',true);
toc
%%
save('Laser','Laser')
