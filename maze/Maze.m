function Maze(filename,TextfilenameTTL,TextfilenameTracking,shankN,chorder,nsamples,peak_position,samplerate)
%
% Maze('DJ52_S3_3O','DJ52_S3_2016_07_01_O_Hardware.txt','DJ52_S3_2016_07_01_O_Track.txt',6,1:10,60,30,30000)
% Maze('DJ52_S1_2Oc','DJ52_S1_2016_06_29_Oc_Hardware.txt','DJ52_S1_2016_06_29_Oc_Track.txt',6,1:10,60,30,30000)

%% 1) making clu

clu_maker(filename,shankN,chorder,nsamples,peak_position,samplerate);

%% 2) making mat file
chN=length(chorder);
loadmulti_maze(filename,TextfilenameTTL,TextfilenameTracking,samplerate,shankN,chN)
