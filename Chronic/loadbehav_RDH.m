%to know the samplerate, type:
%
%read_Intan_RHD2000_file
%frequency_parameters.amplifier_sample_rate

function loadbehav_RDH(filename,samplerate)
%%
% filename='DJ47_1';
% samplerate=30000;
%%

%1) Timestamp, x, lickT, and other digital inputs here
%Timestamp
fileinfo = dir('time.dat');
num_samples = fileinfo.bytes/4; % int32 = 4 bytes
fid = fopen('time.dat', 'r');
timestamp = fread(fid, num_samples, 'int32');
fclose(fid);
timestamp = timestamp / samplerate;
timestamp = timestamp - timestamp(1);

%digital inputs
fileinfo = dir('digitalin.dat');
num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
fid = fopen('digitalin.dat', 'r');
digital_word = fread(fid, num_samples, 'uint16');
fclose(fid);
ndx_change = find(diff(digital_word)~=0)+1;
time_change = timestamp(ndx_change);
digital_word = digital_word(ndx_change);
digital_array = dec2bin(digital_word,16);
%%
if isempty(digital_array)==1
    disp('digitalin is empty')
    Empty=[];
    save('Empty','Empty')
else
 %%
    %lickT
    DIN = 10;
    ndx = 16-DIN;
    lickT = time_change(find(diff(str2num(digital_array(:,ndx)))~=0)+1);
    
    % laser time
    DIN = 4;
    ndx = 16-DIN;
    laserT = time_change(find(diff(str2num(digital_array(:,ndx)))~=0)+1);
    
    
    %positive increment in the belt position
    DIN = 1;
    ndx = 16-DIN;
    pos_incrT = time_change(find(diff(str2num(digital_array(:,ndx)))~=0)+1);
    
    % negative increment in the belt position
    DIN = 2;
    ndx = 16-DIN;
    neg_incrT = time_change(find(diff(str2num(digital_array(:,ndx)))~=0)+1);
    
    % Time when belt position is zero
    DIN = 3;
    ndx = 16-DIN;
    belt0T = time_change(find(diff(str2num(digital_array(:,ndx)))~=0)+1);
    
    
    'step 1 done'
    %%
    TX = Digitalin2tx( pos_incrT, neg_incrT, belt0T );
    
    %2) load treadmill behavior data
    TRDfilename=dir('201*.txt');
    if ~isempty(TRDfilename)
        txlrts=TRDconv_standard(TRDfilename(1).name);
        txlrts(:,1)=txlrts(:,1)/1000; %put time in s
        txlrts(find(diff(txlrts(:,1))==0)+1,1)=txlrts(find(diff(txlrts(:,1))==0)+1,1)+0.0005;
        
        'step 2 done'
        
        %3) sync txrts(:,1) with TDT_Dx_t
        txlrts(:,1)=scaleshift(txlrts(:,1),TX(:,1));
    else
        txlrts = make_txlrts_fromRDH_DJ( TX );  % or make_txlrts_fromRDH( TX )
    end
    'step 3 done'
    %%
    
L=(diff(laserT)>2);

if L(1)==0;
    L=[1;L(:)];
   
else
    L=[0;L(:)];
   
end

k=find(L==1);
laserON=laserT(k(:));



for i=1:length(laserON)
    
   [~,L_ndx(i)]=min(abs(laserON(i)-txlrts(:,1)));
    
end

txlrts(L_ndx,3)=1;

  'step 4 done'  
    
    %%
    TXDTS=EXPANDtxts(txlrts(:,[1 2 5 6]),0.01);
    
    %4) save into .mat file
    behav.txlrts=txlrts;
    behav.TXDTS=TXDTS;
    spkinfo.samplerate=samplerate;
    
    save([filename '.mat'],'-append','behav','spkinfo','laserON')
    
    'all done'
    
end
