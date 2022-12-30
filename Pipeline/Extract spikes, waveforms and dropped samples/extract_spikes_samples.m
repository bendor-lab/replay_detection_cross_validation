function extract_spikes_samples
%extracts sample indices of spikes. calls subroutine extract_kwik_file(kwikfile,tet_number)
cluster_data.spikes.SUA=[];
cluster_data.spikes.MUA=[];

if exist('raw')==7
    cd('raw');
end
disp('folder: ')
pwd
tetrode_list= dir('*.kwik');
for i=1:length(tetrode_list)
    st= strfind(tetrode_list(i).name,'_');
    ed= strfind(tetrode_list(i).name,'.');
    tetrode_num(i)= str2num(tetrode_list(i).name(st+1:ed-1));
end

for ii=1:length(tetrode_list)
    kwikfile= ['tetrode_' num2str(tetrode_num(ii)) '.kwik'];
    % Align timestamps from events (as extracted from Neuralynx) and from spikes(as extracted from klustakwik)
    data_tmp= extract_kwik_file(kwikfile,tetrode_num(ii));
    cluster_data.spikes.SUA= [cluster_data.spikes.SUA; data_tmp.spikes.SUA];
    cluster_data.spikes.MUA= [cluster_data.spikes.MUA; data_tmp.spikes.MUA];
    
    if ~isempty(data_tmp.extracted_clusters.SUA_clusters)
        cluster_data.extracted_clusters.(['TT' num2str(tetrode_num(ii))]).SUA_clusters= data_tmp.extracted_clusters.SUA_clusters;
    end
    if ~isempty(data_tmp.extracted_clusters.MUA_clusters) 
        cluster_data.extracted_clusters.(['TT' num2str(tetrode_num(ii))]).MUA_clusters= data_tmp.extracted_clusters.MUA_clusters;
    end
end  

[~,idx_SUA_sorted]= sort(cluster_data.spikes.SUA(:,2)); %rearranging so that samples are in the correct order
cluster_data.spikes.SUA= cluster_data.spikes.SUA(idx_SUA_sorted,:);
[~,idx_MUA_sorted]= sort(cluster_data.spikes.MUA(:,2)); %rearranging so that samples are in the correct order
cluster_data.spikes.MUA= cluster_data.spikes.MUA(idx_MUA_sorted,:);
cluster_data.session= pwd;
cluster_data.clustering= 'klustakwik-klustaviewa'; % 'klustakwik-phy'

if exist('raw')==7
    cd ..
end
save('cluster_data.mat','cluster_data','-v7.3');
end




function data= extract_kwik_file(kwikfile,tet_number)

% load clustered spike timestamps
kwikinfo = h5info(kwikfile);


time_samples = h5read(kwikfile,'/channel_groups/0/spikes/time_samples');
time_samples_clusters = h5read(kwikfile, '/channel_groups/0/spikes/clusters/main');
clusters = unique(time_samples_clusters);

% For each session, get the relevant spikes
 timestamps =  double(time_samples); %
% Where lims is simply an matrix of nx2 numbers that indicate the sample number in the concatenated da file that are at the beginning and the end of the relevant session 


% Build a boolean vector with 1 by a cluster id if it is a good cluster
good_clusters_timestamps=[];
good_clusters_id=[];
MUA_clusters_timestamps=[];
MUA_clusters_id=[];
original_good_cluster_id=[];
original_MUA_cluster_id=[];

for i = 1:length(clusters)
        
    cluster_attr_path = strcat('/channel_groups/0/clusters/main/', num2str(clusters(i)));
    cluster_group = h5readatt(kwikfile, cluster_attr_path, 'cluster_group');
    %cluster_group 2 is Good
    if cluster_group == 2
        indx=find(time_samples_clusters==clusters(i));
        good_clusters_timestamps= [good_clusters_timestamps; timestamps(indx)];
        good_clusters_id= [good_clusters_id; (1000*tet_number+double(clusters(i)))*ones(size(timestamps(indx)))];
        original_good_cluster_id= [original_good_cluster_id; double(clusters(i))];
    elseif cluster_group == 1 %MUA
        indx=find(time_samples_clusters==clusters(i));
        MUA_clusters_timestamps=[MUA_clusters_timestamps; timestamps(indx)];
        MUA_clusters_id= [MUA_clusters_id; (1000*tet_number+double(clusters(i)))*ones(size(timestamps(indx)))];
        original_MUA_cluster_id= [original_MUA_cluster_id; double(clusters(i))];
    end
end
 

if ~isempty(good_clusters_id)
    data.spikes.SUA(:,1)=good_clusters_id;
    data.spikes.SUA(:,2)=good_clusters_timestamps;
    data.extracted_clusters.SUA_clusters= [unique(good_clusters_id) original_good_cluster_id];
else
    data.spikes.SUA=[]; 
    data.extracted_clusters.SUA_clusters= [];
end

if ~isempty(MUA_clusters_id)
    data.spikes.MUA(:,1)=MUA_clusters_id;
    data.spikes.MUA(:,2)=MUA_clusters_timestamps;
    data.extracted_clusters.MUA_clusters= [unique(MUA_clusters_id) original_MUA_cluster_id];
else
    data.spikes.MUA=[];
    data.extracted_clusters.MUA_clusters= [];
end

end






