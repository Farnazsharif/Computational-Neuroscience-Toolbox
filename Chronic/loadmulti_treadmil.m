function loadmulti_treadmil(samplingrate)

% samplingrate=30000;

samplingrate=30000;
load('TRD1_Belt.mat')
Beltinfo_1=Beltinfo;
load('TRD2_Belt.mat')
Beltinfo_2=Beltinfo;
load('TRD3_Belt.mat')
Beltinfo_3=Beltinfo;
Beltinfo=[];

name='uncombined';
dirinfo1=dir(name);
[ds,~]=size(dirinfo1); %%check the number of folders


fln1=0;
for hh=1:ds
    if isempty(strfind(dirinfo1(hh).name,'SF'))==1 %% detect folder name 'DJ'
        fln1=fln1+1;
    end
end


cd(name)

for j=1:ds-fln1

 
    if isempty(strfind(dirinfo1(j+fln1).name,'T'))==0
        
        cd(dirinfo1(j+fln1).name)
        filename=dirinfo1(j+fln1).name
        
        if isempty(strfind(filename,'TRD1'))==0
        Beltinfo=Beltinfo_1;
        else
        Beltinfo=Beltinfo_2;
        end
        
        save([filename '.mat'],'-append','Beltinfo')
        
        
        %% step 1)
        
        loadbehav_RDH_farn(filename,samplingrate);
        'loadbehav_RDH is done'
        if exist('Empty.mat')==0
        %% step 2 )
        % loadmulti3('DJ57_S1_1T',6,10,5,2)
        loadmulti3_chronic(filename)
        'loadmulti3 is done'
        %% step 3 )
        % filename='DJ56_S2_1T';
        % binN=100;
        % speed_threshold=5;

%         xttsc=makexttsc_DJ(behav.TXDTS,spk,binN,speed_threshold,G_C); for Dajung
        
%       xttsc=makexttsc(TXDTS,spk,xbinNumber,SpeedThreshold,reocrding,filename)% for Farnaz
        load([filename '.mat'])
        xttscT = makexttsc(behav.TXDTS,spk,100,0,'ch',filename);
        
        save([filename 'PlaceField_Total.mat'],'xttscT')
        
        'xttsT is done'
        end
         cd ..
    end
      
end
cd ..

