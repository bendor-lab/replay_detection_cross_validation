
function extract_spikes_phy


rat= input('Rat name: ','s'); % J-GRE, J-BLA
date=input('Date and time: ','s'); % eg. 2017-08-17_15-22-00
session= input('Session (e.g.Day XX): ','s'); %Task, Sleep1, Sleep2 or Day XX
project=input('VC or HIPP: ','s');
if strcmp(project,'HIPP')
    path= (sprintf('%s%s%s','X:\Bendorlab\Drobo\Lab Members\Marta\Behavioural Data\',project,'\',rat,'\CheetahData\',date));
elseif strcmp(project,'VC')
    path= (sprintf('%s%s%s','X:\Bendorlab\Drobo\Lab Members\Marta\Behavioural Data\',project,'\',rat,'\CheetahData\',date));
end
 
analysis_save_path = (sprintf('%s',pwd,'\','Raw\'));

%% Extract spike data
manual_tetrodes=input('Please, enter which tetrodes you would like to process (e.g.[4,5,7]): ');
[classified_clu] = kwik_to_mat(rat,date,session,project,manual_tetrodes); % ONLY PHY: Extracts classification of clusters from excel file

if length(manual_tetrodes)>1
    for ii=1:length(manual_tetrodes)
        kwikfile= [path '\' 'TT' num2str(manual_tetrodes(ii)) '.kwik'];
        tetrode = strcat('TT',mat2str(manual_tetrodes(ii)));
        [data,timestamp_events_start]=extract_spikes_cellID(classified_clu,path,kwikfile,tetrode);
        if ~exist('all_data','var')
            if isfield(data.spikes,'SUA')
                all_data.spikes.SUA=data.spikes.SUA;
            end
            if isfield(data.spikes,'MUA')
                all_data.spikes.MUA=data.spikes.MUA;
            end
            all_data.clustering= data.clustering;
            all_data.positions=data.positions;
            all_data.events=data.events;
        else
            if isfield(data.spikes,'SUA')
                if ~isfield(all_data.spikes,'SUA')
                    all_data.spikes.SUA=data.spikes.SUA;
                else
                    all_data.spikes.SUA=[all_data.spikes.SUA; data.spikes.SUA];
                end
            end
            if isfield(data.spikes,'MUA')
                if ~isfield(all_data.spikes,'SUA')
                    all_data.spikes.MUA=data.spikes.MUA;
                else
                    all_data.spikes.MUA=[all_data.spikes.MUA; data.spikes.MUA];
                end
            end
        end
    end
else
    kwikfile= [path '\' 'TT' num2str(manual_tetrodes) '.kwik'];
    tetrode = strcat('TT',mat2str(manual_tetrodes));
    [data,timestamp_events_start]=extract_spikes_cellID(classified_clu,path,kwikfile,tetrode);
end

if exist('all_data','var')
    data=all_data;
end

% Save in structure
cluster_data = data;
cluster_data.extracted_clusters = classified_clu;
%save(sprintf('%s',save_path,'\',session,'_',date,'\',rat,'_',date,'_','data.mat'),'data','-v7.3');
save(sprintf('%s',analysis_save_path,'cluster_data.mat'),'cluster_data','-v7.3');

end 


function [data,timestamp_events_start]=extract_spikes_cellID(classified_clu,path,kwikfile,tetrode)
%MH_10.10.18

%%%%%%% PARAMETERS TO CHANGE  %%%%%%%%%%

clustering= 'klustakwik-phy';
SR= 30000;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extracts position, events and spikes
%Folder with video and events files

if ~exist('data','var')
    files_cluster= dir(path);
    clear filename
    for i=3:length(files_cluster)
        file_name= files_cluster(i).name;
        if ~isempty(strfind(file_name,'Events'))%finds event file
            evt_file= strcat('\',file_name);
        end
        if ~isempty(strfind(file_name,'.nvt'))  %finds video file
            p_file=strcat('\',file_name);
        end
    end
    position_data= getPositions_MM([path p_file]);
    events_data= getRawTTLs_MM([path evt_file]);
    
    %saves data
    data.events= events_data;
    data.positions= position_data;
else 
    position_data=data.positions;
    events_data=data.events;
end 
% % Align timestamps from events (as extracted from Neuralynx) and from spikes(as extracted from klustakwik)
timestamp_events_start= events_data.timestamps(1)*1e-6;

%% Extracts spikes timestamps and ID from PHY

time_samples = h5read(char(kwikfile),'/channel_groups/0/spikes/time_samples');
time_samples_clusters = h5read(char(kwikfile), '/channel_groups/0/spikes/clusters/main');
clusters = unique(time_samples_clusters);
timestamps = timestamp_events_start + double(time_samples)./SR; 

current_tetrode=char(strsplit(tetrode,'.'));
tt_num=str2num(current_tetrode(3:end));
good_clusters_timestamps=[];
good_clusters_id=[];
good_clusters_samples=[];
MUA_clusters_timestamps=[];
MUA_clusters_id=[];
MUA_samples = [];
for i=1:length(clusters)

    indx=find(time_samples_clusters==clusters(i));
    
    if ~isempty(cell2mat(strfind(fieldnames(classified_clu.(sprintf('%s',tetrode))),'MUA_clusters')))
        if ismember(clusters(i),classified_clu.(sprintf('%s',tetrode)).MUA_clusters)
            MUA_samples = [MUA_samples; time_samples(indx)];
            MUA_clusters_timestamps=[MUA_clusters_timestamps; timestamps(indx)];
            cell_id=(1000*tt_num)+clusters(i); % nomenclature e.g. 7036 -> tetrode 7 cluster 36
            MUA_clusters_id= [MUA_clusters_id; double(cell_id)*ones(size(timestamps(indx)))];
        end
    end
    if ~isempty(cell2mat(strfind(fieldnames(classified_clu.(sprintf('%s',tetrode))),'SUA_clusters')))
        if ismember(clusters(i),classified_clu.(sprintf('%s',tetrode)).SUA_clusters)
            good_clusters_samples = [good_clusters_samples; time_samples(indx)];
            good_clusters_timestamps=[good_clusters_timestamps; timestamps(indx)];
            cell_id=(1000*tt_num)+clusters(i);
            good_clusters_id= [good_clusters_id; double(cell_id)*ones(size(timestamps(indx)))];
        end
end
end

%% Saving data
if ~isfield(data,'spikes')
    if ~isempty(good_clusters_id)
        data.spikes.SUA(:,1)=good_clusters_id;
        data.spikes.SUA(:,2)=good_clusters_samples;
        data.spikes.SUA(:,3)=good_clusters_timestamps;
    end
    if ~isempty(MUA_clusters_id)
        data.spikes.MUA(:,1) = MUA_clusters_id;
        data.spikes.MUA(:,2) = MUA_samples;
        data.spikes.MUA(:,3) = MUA_clusters_timestamps;
    end
    data.clustering= clustering;
end

end