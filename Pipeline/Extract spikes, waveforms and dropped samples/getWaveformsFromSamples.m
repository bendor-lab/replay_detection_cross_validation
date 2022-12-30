%% Waveform extraction from DAT files & Half-width calculation

% Loads samples from each cluster extracted from kwik files and waveform
% informatiion stored in DAT files (which should be saved in 'raw' folder). Calculates mean waveform for each
% cluster,half-with amplitude and peak.

% Output: structure containing cluster information and extracted waveform
% per cluster

function allclusters_waveform = getWaveformsFromSamples
load cluster_data 
load extracted_clusters % need ID conversion
load extracted_events.mat
allclusters_waveform = [];

% Parameters
parameters = list_of_parameters;
nChannels = parameters.nChannels;
nSamplesForSpike = parameters.nSamplesForSpike;
SR = parameters.SR; %30kHz
timestamp_events_start= events_data.timestamps(1)*1e-6;

% Find tetrodes with units
tetrodes = fieldnames(cluster_data.extracted_clusters);

for thisTetrode = 1:length(tetrodes)
    
    TT = cell2mat(tetrodes(thisTetrode));
    TT_num = strsplit(TT,'TT');
    TT_num = cell2mat(TT_num(2));
    if ~isfield(cluster_data.extracted_clusters.(sprintf('%s',TT)), 'SUA_clusters') %if this tetrode doesn't have SUA, skip
        continue
    end
    
    switch cluster_data.clustering
        case 'klustakwik-phy'
            TT_clusters = cluster_data.extracted_clusters.(sprintf('%s',TT)).SUA_clusters+(str2num(TT_num)*1000);
        case 'klustakwik-klustaviewa'
            TT_clusters = cluster_data.extracted_clusters.(sprintf('%s',TT)).SUA_clusters;
        otherwise
            disp('clustering software not found');
            keyboard;
    end
    
    if exist('raw')==7 %find DAT file for this tetrode in 'raw' folder
        OldFolder = cd('raw');
        % check for their name - some directories will have 'tetrode_1.dat', others
        % will have 'TT1.dat'
        these_files= dir;
        these_files= {these_files.name};
        dat_files= find(contains(these_files,'.dat'));
        if contains(these_files(dat_files(1)),'tetrode')
            TT_name= 'tetrode_';
        else
            TT_name= 'TT';
        end
        DAT_file = strcat([TT_name TT_num],'.dat');
        % Create a map of the file(its in int16 format, with 1 sample from each channel ie Ch1:Ch16, Ch1:Ch16 etc)
        thisMemMap = memmapfile(DAT_file,'Format','int16');
        cd(OldFolder)
    end
    TT_clusters_samples = cluster_data.spikes.SUA(find(cluster_data.spikes.SUA(:,1) == TT_clusters(1),1,'first'):find(cluster_data.spikes.SUA(:,1) == TT_clusters(end),1,'last'),2);
    
    timestamps = timestamp_events_start + double(TT_clusters_samples)./SR;
    spikeIndices = nSamplesForSpike(1)*nChannels+1:nSamplesForSpike(2)*nChannels;
    maxSamples = length(thisMemMap.Data)/nChannels;
    
    clear cluster_waveform
    clear widthWaveform
    
    for thisClust = 1:length(TT_clusters)  
        % Get the relevant samples
        clustID = TT_clusters(thisClust);
        clusterID_indices = cluster_data.spikes.SUA(:,1) == clustID; %find cluster ID indices for this cluster
        cluster_samples = cluster_data.spikes.SUA(clusterID_indices,2); %find samples for this cluster
        
        if length(cluster_samples) > 1
            % Extract time samples around spike
            indices1 = repmat(spikeIndices,length(cluster_samples),1);
            indices2 = repmat(cluster_samples*nChannels-nChannels,1,length(spikeIndices));
            addressIndices = indices1+double(indices2);
            addressIndices(addressIndices<1) = 1;
            addressIndices(addressIndices>maxSamples) = maxSamples;
            tdata = thisMemMap.Data(addressIndices); % finds samples in DAT file
            tdata = double(tdata);
            
            % Subtract means from each waveform before averaging
            tdata = tdata-repmat(mean(tdata,2),1,size(tdata,2));
            tmeanWaveform = mean(tdata,1);
            tstdWaveform = std(tdata,[],1);
            meanWaveform{thisClust} = reshape(tmeanWaveform,nChannels,length(tmeanWaveform)/nChannels);
            stdWaveform{thisClust} = reshape(tstdWaveform,nChannels,length(tstdWaveform)/nChannels);
            nWaveform{thisClust} = length(cluster_samples);
            widthWaveform{thisClust}= waveform_width(meanWaveform{thisClust},nChannels);
            
            % Save cluster data in structure
            cluster_waveform(thisClust).clust_ID = clustID; % nomenclature e.g. 7036 -> tetrode 7 cluster 36
            cluster_waveform(thisClust).converted_ID = clusters.id_conversion(find(clusters.id_conversion(:,2)==clustID),1);
            cluster_waveform(thisClust).samples = double(cluster_samples);
            cluster_waveform(thisClust).meanWaveform = meanWaveform{thisClust} ;
            cluster_waveform(thisClust).stdWaveform = stdWaveform{thisClust};
            cluster_waveform(thisClust).nWaveform = nWaveform{thisClust};
            cluster_waveform(thisClust).half_width = widthWaveform{thisClust}.width;
            cluster_waveform(thisClust).max_peak = widthWaveform{thisClust}.max_peak;
            cluster_waveform(thisClust).max_peak_CSC = (widthWaveform{thisClust}.max_peak_CSC)+(4*(str2num(TT_num)-1)); % calculates actual CSC number
            
            
        else
            meanWaveform{thisClust} = zeros(nChannels,length(timestamps));
            stdWaveform{thisClust} = zeros(nChannels,length(timestamps));
            nWaveform{thisClust} = length(cluster_samples);
            widthWaveform{thisClust}= zeros(nChannels,length(timestamps));
        end
    end
    
    allclusters_waveform = [allclusters_waveform cluster_waveform];
        
end

save('extracted_waveforms','allclusters_waveform')
end


function widthWaveform = waveform_width(spike_samples,nChannels)
% Parameters
parameters = list_of_parameters;
SR = parameters.SR; %30kHz
time_btw_samples=1/SR; %in seconds
CSC_peaks = [];
for j = 1:nChannels
    [peak,peak_index] = min(spike_samples(j,:)); %finds the peak of the spike (which is the min since the spike is reversed)
    CSC_peaks = [CSC_peaks; peak,peak_index];
end
[max_peak,CSC_index] = min(CSC_peaks(:,1)); %finds the max peak within the 4 CSC and keeps the CSC index
max_peak_index = CSC_peaks(CSC_index,2);
max_waveform_samples = spike_samples(CSC_index,:);
[~,posttrough_idx] = max(max_waveform_samples(max_peak_index:end)); %finds max value (wave trough) after the peak
posttrough_idx= posttrough_idx+max_peak_index-1; %finds the actual index of the trough
width = (posttrough_idx-max_peak_index)*time_btw_samples; %calculates time (half_width) between peak and trough
widthWaveform.width = width;
widthWaveform.max_peak = max_peak;
widthWaveform.max_peak_CSC = CSC_index;  % CSC with max peak (from the 4 CSC that a tetrode has)
end


