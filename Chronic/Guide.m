%% 1) clu_seperator : first copy .par file from one of the experiments

% Cluseperator(filename,shankN,chN)
Cluseperator('DJ55_S2C',6,10);
'Cluseperator is done'

%% 2) making clu in each folder

% clu_maker(shankN,chorder,nsamples) 
clu_maker(6,1:10,60);
'clu_maker is done'

%% 3) making .mat file in treadmil

% loadmulti_treadmil(shankN,samplingrate,chN,binN,speed_threshold,duration_threshold)
loadmulti_treadmil(6,30000,10,100,5,2)
'loadmulti_treadmil is done'

%% 4) making .mat file in maze

% loadmulti_maze(samplerate,shankN,chN);
loadmulti_maze(30000,6,10);
'loadmulti_maze is done'

%% 5) belt
%first open appendBeltinfo file and fulfil the details
appendBeltinfo('DJ55_S1_3T_new')

%% 6) plot

% plot_PF_rs_W(Celln,smooth,Speed_noise_maze,Speed_threshold_maze)
plot_PF_rs_W(100,10,50,0)

%%  missed cells

G_C(find(ismember(G_C,G)==0))

