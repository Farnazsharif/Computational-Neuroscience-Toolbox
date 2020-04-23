clear

% seperating beltinfo

name='uncombined';
dirinfo1=dir(name);
[ds,~]=size(dirinfo1); %%check the number of folders


fln1=0;
for hh=1:ds
    if isempty(strfind(dirinfo1(hh).name,'DJ'))==1 %% detect folder name 'DJ'
        fln1=fln1+1;
    end
end
cd(name)

for j=1:ds-fln1
    
    if isempty(strfind(dirinfo1(j+fln1).name,'T'))==0
        
        cd(dirinfo1(j+fln1).name)
        filename=dirinfo1(j+fln1).name
        load([filename '.mat'])
        save('BeltInfo.mat','BeltInfo')
        cd ..
    end
end
cd ..

%%

