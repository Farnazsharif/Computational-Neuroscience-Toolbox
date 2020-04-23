eeg=[];

for i=1:2:3
eeg = readmulti('FM05_1.lfp',128,i);
end
eeg;
