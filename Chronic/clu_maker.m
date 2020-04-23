
function clu_maker(shankN,chorder,nsamples)

% clu_maker(6,1:10,60)
% shankN=6;
% chorder=1:10;
% nsamples=60;

name='uncombined';
dirinfo1=dir(name);
[ds,~]=size(dirinfo1); %%check the number of folders


fln1=0;
for hh=1:ds
    if isempty(strfind(dirinfo1(hh).name,'SF'))==1 %% detect folder name 'DJ'
        fln1=fln1+1;
    end
end

for k=1:(ds-fln1)
    
    cd uncombined
    cd(dirinfo1(k+fln1).name)
    filename=dirinfo1(k+fln1).name;
    
    for shankID=1:shankN
        Cellnumb_pershank=[];
        nchannels = length(chorder);
        fp = fopen([filename '.spk.' num2str(shankID)], 'r');
        spkW = fread(fp, [nchannels, inf], 'short');
        nspike = size(spkW, 2)/nsamples;
        spkW = reshape(spkW, [nchannels, nsamples, nspike]);
        spkW = spkW(chorder,:,:);
        
        load([filename '.mat'])
        load([filename '_clu_' num2str(shankID) '.mat']);
        
        Cellnumb_pershank=length(find(fix(G_C/100)==shankID));
        clu.spk_mean=[];
        clu.spk_std=[];
        clu.g=[];
        CellspkW=[];
        clu.g=g;
        mclu=[];
        
        for i=1: Cellnumb_pershank
            
            CellspkW=spkW(:,:,clu.g==i);
            clu.spk_mean(:,:,i)=mean(CellspkW,3);
            clu.spk_std(:,:,i)=std(CellspkW,0,3);
            
        end
        save([filename '_clu_' num2str(shankID) '.mat'],'-append','clu','mclu')
    end
    
    cd ..
    cd ..
end