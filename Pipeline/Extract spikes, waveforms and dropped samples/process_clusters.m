function process_clusters
% extracts cluster data from cluster_data.mat, and stores this into two
% variables (SUA for single units, and MUA for multiunits). The variable
% lists the contributing unit id and the timestamps/index of the spike. an
% id conversion variable maps the unit id to the original id in the
% cluster_data.mat file.

% spike times must be in the form of the recording timestamp index (not time)

% input is empty if data is taken from cluster_data.mat, othewise it must
%     be a structure "spikes" with a .SUA and .MUA field containing the data
current_directory=pwd;
if exist('raw/cluster_data.mat','file')==2
    load('raw/cluster_data.mat');
elseif exist('cluster_data.mat','file')==2
    load('cluster_data.mat');
else
    disp('cluster_data.mat cannot be found')
end

if exist('data')==1
cluster_data=data;  %for marta's file format
clear data;
end

if exist('cluster_data')==1
    load extracted_dropped_samples; %loads 3 variables- 'dropped_samples','error_times','Time'
    parameters= list_of_parameters;
    spikes_SUA= cluster_data.spikes.SUA;    spikes_MUA= cluster_data.spikes.MUA;
    spike_id_SUA= spikes_SUA(:,1);          spike_id_MUA= spikes_MUA(:,1);
    spike_index_SUA= spikes_SUA(:,2);       spike_index_MUA= spikes_MUA(:,2);
    
    %determine nunber of unique units and rename cluster id from 1 to n
    clusters_id_SUA= unique(spike_id_SUA);
    clusters_id_MUA= unique(spike_id_MUA);
    spike_id_SUA_ordered= NaN(size(spike_id_SUA));
    
    % SUA-convert cluster id to values 1 to n
    for i=1:length(clusters_id_SUA)
        spike_id_SUA_ordered(find(spike_id_SUA==clusters_id_SUA(i)))=i;
        id_conversion_SUA(i,1:2)=[i clusters_id_SUA(i)];
    end
    
    %remove dropped samples from SUA
    spike_times_SUA= Time(spike_index_SUA);
    dropped_samples_index_SUA=[];
    for i=1:length(dropped_samples)
        dropped_samples_index_SUA=[dropped_samples_index_SUA find(spike_index_SUA==dropped_samples(i))'];
    end
    % for backwards compatibility, keep only SUA in extracted_clusters.mat
    clusters.spike_id=spike_id_SUA_ordered;
    clusters.spike_times=spike_times_SUA;
    clusters.spike_id(dropped_samples_index_SUA)=[];
    clusters.spike_times(dropped_samples_index_SUA)=[];
    clusters.id_conversion= id_conversion_SUA;  %to map new id to old cluster id
    
    %remove dropped samples from MUA
    spike_times_MUA= Time(spike_index_MUA);
    dropped_samples_index_MUA=[];
    for i=1:length(dropped_samples)
        dropped_samples_index_MUA=[dropped_samples_index_MUA find(spike_index_MUA==dropped_samples(i))'];
    end
    spike_times_MUA(dropped_samples_index_MUA)=[];
    
    % Sort SUA by time just in case it hasn't been done before
    [~,index] = sort(clusters.spike_times);
    clusters.spike_id = clusters.spike_id(index);
    clusters.spike_times = clusters.spike_times(index);
    
    % add SUA to MUA .mat
    all_MUA_id= [clusters.spike_id; spike_id_MUA];
    all_MUA_time= [clusters.spike_times; spike_times_MUA];
    [all_MUA_time, idx_sorted]= sort(all_MUA_time);
    MUA.spike_times= all_MUA_time;
    MUA.spike_id= all_MUA_id(idx_sorted);
    
    %smooth MUA activity with a gaussian kernel
    time_step=0.001; %1 ms timestep
    MUA.time_bins_edges= min(Time):(time_step):max(Time);
    w=gausswin(parameters.MUA_filter_length,parameters.MUA_filter_alpha); %41,2
    w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width
    MUA.zscored=zscore(filtfilt(w,1,histcounts(MUA.spike_times,MUA.time_bins_edges)));
    MUA.time_bins=(MUA.time_bins_edges(1:end-1))+time_step/2;
    
    %save extracted SUA and smoothed MUA activity
    save('MUA_clusters.mat','MUA');
    save('extracted_clusters.mat','clusters');
else
    disp('error finding cluster data')
end
end