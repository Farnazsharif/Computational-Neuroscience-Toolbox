function units= ConvertPhyNPY2Matlab(basepath,basename)
% under construction

% Converts results of manual sorting with KiloSort in Phy into Matlab.
% 
% Inputs:
%   basepath -  directory path to the main recording folder with .dat and .xml
%               as well as shank folders made by makeProbeMapKlusta2.m (default is
%               current directory matlab is pointed to)
%   basename -  shared file name of .dat and .xml (default is last part of
%               current directory path, ie most immediate folder name)

% Outputs:
%   units - a struct (matlab) containing information regarding clustered
%   spike activities
%                 units(j).ids = double array of spike indexes for the unit j;
%                 units(j).ts = double array of spike times for the unit j;
%                 units(j).kwik_id = id of unit j after manual clustring;
%                 units(j).shank = shank number;
%                 units(j).total = total number of firing for unit j;
%                 units(j).amplitudes = double array of amplitudes of all spikes correc\sponding to unit j;
%
%Buzsaki Lab 2018


if ~exist('basepath','var')
    [~,basename] = fileparts(cd);
    basepath = cd;
end


%% code from Peter Peterson
clustering_path=basepath;
disp('Loading Phy clustered data')
    spike_cluster_index = readNPY(fullfile(clustering_path, 'spike_clusters.npy'));
    spike_times = readNPY(fullfile(clustering_path, 'spike_times.npy'));
    spike_amplitudes = readNPY(fullfile(clustering_path, 'amplitudes.npy'));
    spike_channel_shanks = readNPY(fullfile(clustering_path, 'channel_shanks.npy'));
    spike_peak_channel = readNPY(fullfile(clustering_path, 'peak_channel.npy'));
    spike_clusters = unique(spike_cluster_index);
    filename1 = fullfile(clustering_path,'cluster_group.tsv');
    filename2 = fullfile(clustering_path,'cluster_groups.csv');
    if exist(filename1) == 2
        filename = filename1;
    elseif exist(filename2) == 2
        filename = filename2;
    else
        disp('No manual clustering data found')
    end
    delimiter = '\t';
    startRow = 2;
    formatSpec = '%f%s%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    units = [];
    j = 1;
    for i = 1:length(dataArray{1})
            if strcmp(dataArray{2}{i},'good')
                units(j).ids = find(spike_cluster_index == dataArray{1}(i))';
                units(j).ts = double(spike_times(units(j).ids))';
                units(j).kwik_id = dataArray{1}(i);
                units(j).total = length(units(j).ts);
                units(j).amplitudes = double(spike_amplitudes(units(j).ids))';
                units(j).peak_channel=double(spike_peak_channel(find(spike_clusters==units(j).kwik_id)));
                units(j).shank =double(spike_channel_shanks(units(j).peak_channel+1)) ;   %+1 to correct for channel number mis-match
                j = j+1;
            end
    end
    
 