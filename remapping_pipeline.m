% BATCH ANALYSIS PIPELINE
function remapping_pipeline(option,RExp)

% EXTRACT SPIKES, WAVEFORMS AND DROPPED SAMPLES
%use one of these depending on how you clustered
if strcmp(option,'EXTRACT_SPIKES_MARTA')
    extract_spikes_phy; % for clusters extracted from PHY
    extract_events;
    extract_dropped_samples;
    process_clusters;
    getWaveformsFromSamples;
end    
if strcmp(option,'EXTRACT_SPIKES_MARGOT')
    extract_spikes_samples; % for clusters extracted from Klustaviewa
    extract_events;
    extract_dropped_samples;
    process_clusters;
    getWaveformsFromSamples;
end
if strcmp(option,'EXTRACT_VIDEO')
    % EXTRACT VIDEO
    extract_video;   %if too noisy, use "extract_video('targets');
    frame_grab;
end
if strcmp(option,'EXTRACT_POSITION')
    % EXTRACT POSITIONS
    disp('processing position data')
    process_positions_PRE;
    process_positions_POST([]);  %or process_positions_POST([1,3;2,4])
    % plot_cleaning_steps;
end
if strcmp(option,'CSC')
    % EXTRACT CSC
    disp('processing CSC data')
    parallel_extract_PSD('sleep');
    best_channels = determine_best_channel('hpc');
    extract_CSC('hpc');
end
if strcmp(option,'PLACE FIELDS')
    % EXTRACT PLACE FIELDS
    disp('processing place_field data')
    parameters=list_of_parameters;
    calculate_place_fields(parameters.x_bins_width_bayesian);
    calculate_place_fields(parameters.x_bins_width);
    % plot_place_fields;
end
if strcmp(option,'SLEEP')
    % EXTRACT SLEEP
    sleep_stager('manual');
end
if strcmp(option,'extract bayesian and replay')
    disp('Bayesian and extract replay')
    spike_count([],[],[],'Y');
    bayesian_decoding([],[],'Y');
    extract_replay_events; %finds onset and offset of replay events
end
if strcmp(option,'REPLAY')
    % EXTRACT REPLAY EVENTS and BAYESIAN DECODING
    disp('processing replay events')
    replay_decoding; %extract and decodes replay events
    
    % SCORING METHODS: TEST SIGNIFICANCE ON REPLAY EVENTS
    disp('scoring replay events')
    scored_replay = replay_scoring([],[0 1 0 1]); % weighted corr   &  spearman
    save scored_replay scored_replay;
    
    % RUN SHUFFLES
    disp('running shuffles')
    num_shuffles=1000;
    analysis_type=[0 1 0 1];  %just weighted correlation and pacman
    load decoded_replay_events
    p = gcp; % Starting new parallel pool
    shuffle_choice={'PRE spike_train_circular_shift','PRE place_field_circular_shift', 'POST place bin circular shift'};
    tic
    if ~isempty(p)
        for shuffle_id=1:length(shuffle_choice)
            shuffle_type{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events);
        end
    else
        disp('parallel processing not possible');
        for shuffle_id=1:length(shuffle_choice)
            shuffle_type{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events);
        end
    end
    save shuffled_tracks shuffle_type;
    disp('time to run shuffles was...');
    toc
    
    % Evaluate significance
    load('scored_replay.mat');
    load('shuffled_tracks.mat');
    scored_replay= replay_significance(scored_replay, shuffle_type);
    save scored_replay scored_replay
    
    %%%%%%analyze segments%%%%%%%%%%
    % splitting replay events
    tic
    replay_decoding_split_events;
    load decoded_replay_events_segments;
    scored_replay1 = replay_scoring(decoded_replay_events1,[0 1 0 1]);
    scored_replay2 = replay_scoring(decoded_replay_events2,[0 1 0 1]);
    save scored_replay_segments scored_replay1 scored_replay2;
    num_shuffles=1000;
    analysis_type=[0 1 0 1];  %just weighted correlation and spearman
    
    load decoded_replay_events_segments;
    p = gcp; % Starting new parallel pool
    if ~isempty(p)
        for shuffle_id=1:length(shuffle_choice)
            shuffle_type1{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events1);
            shuffle_type2{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events2);
        end
    else
        disp('parallel processing not possible');
        for shuffle_id=1:length(shuffle_choice)
            shuffle_type1{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events1);
            shuffle_type2{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events2);
        end
    end
    
    save shuffled_tracks_segments shuffle_type1 shuffle_type2;
    
    load scored_replay_segments; load shuffled_tracks_segments;
    scored_replay1=replay_significance(scored_replay1, shuffle_type1);
    scored_replay2=replay_significance(scored_replay2, shuffle_type2);
    save scored_replay_segments scored_replay1 scored_replay2
    disp('running segments took...')
    toc
    
    %%%% analyze significant replays from whole events and segmented events
    % only using replay events passing threshold for ripple power
    
    if isempty(RExp)
        % no reexposure
        number_of_significant_replays(0.05,3,'wcorr',[]); % pval, ripple zscore, method, reexposure
        number_of_significant_replays(0.05,3,'spearman',[]); % pval, ripple zscore, method, reexposure
    else
        % reexposure
        number_of_significant_replays(0.05,3,'wcorr',2); % pval, ripple zscore, method, reexposure
        number_of_significant_replays(0.05,3,'spearman',2); % pval, ripple zscore, method, reexposure
    end
end
if strcmp(option,'SHUFFLES') % SHUFFLES & SCORING ONLY
        disp('scoring replay events')
        load('decoded_replay_events.mat');
        scored_replay = replay_scoring(decoded_replay_events,[0 1 0 0]); % weighted corr  
        save scored_replay scored_replay;

        % RUN SHUFFLES
        tic;
        disp('running shuffles')
        num_shuffles=1000;
        analysis_type=[0 1 0 0];  %just weighted correlation 
        p = gcp; % Starting new parallel pool
        shuffle_choice={'PRE spike_train_circular_shift','PRE place_field_circular_shift', 'POST place bin circular shift'};
        tic
        if ~isempty(p)
            for shuffle_id=1:length(shuffle_choice)
                shuffle_type{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events);
            end
        else
            disp('parallel processing not possible');
            for shuffle_id=1:length(shuffle_choice)
                shuffle_type{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events);
            end
        end
        save shuffled_tracks shuffle_type;
        disp('time to run shuffles was...');
        toc
        % Evaluate significance
        load('scored_replay.mat');
        load('shuffled_tracks.mat');
        scored_replay= replay_significance(scored_replay, shuffle_type);
        save scored_replay scored_replay

        %%%%%%analyze segments%%%%%%%%%%
        % splitting replay events
        tic
        replay_decoding_split_events('replay_rate_shuffle');
        load decoded_replay_events_segments;
        scored_replay1 = replay_scoring(decoded_replay_events1,[0 1 0 0]);
        scored_replay2 = replay_scoring(decoded_replay_events2,[0 1 0 0]);
        save scored_replay_segments scored_replay1 scored_replay2;
        num_shuffles=1000;
        analysis_type=[0 1 0 0];  %just weighted correlation

        load decoded_replay_events_segments;
        p = gcp; % Starting new parallel pool
        if ~isempty(p)
            for shuffle_id=1:length(shuffle_choice)
                shuffle_type1{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events1);
                shuffle_type2{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events2);
            end
        else
            disp('parallel processing not possible');
            for shuffle_id=1:length(shuffle_choice)
                shuffle_type1{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events1);
                shuffle_type2{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events2);
            end
        end

        save shuffled_tracks_segments shuffle_type1 shuffle_type2;

        load scored_replay_segments; load shuffled_tracks_segments;
        scored_replay1=replay_significance(scored_replay1, shuffle_type1);
        scored_replay2=replay_significance(scored_replay2, shuffle_type2);
        save scored_replay_segments scored_replay1 scored_replay2
        disp('running segments took...')
        toc
end
if strcmp(option,'ANALYSIS')
    % load scored_replay;
    if isempty(RExp)
        % no reexposure
        sort_replay_events([],'spearman');
        sort_replay_events([],'wcorr');
    else
        % reexposure
        sort_replay_events([1,3;2,4],'spearman');
        sort_replay_events([1,3;2,4],'wcorr');
    end
    
    
    % individual folders
    rate_remapping_TRACK_PAIRS([],'spearman');
    rate_remapping_TRACK_PAIRS([],'wcorr');
    
    % all sessions
%     load('folders_to_process_remapping.mat')
%     rate_remapping_TRACK_PAIRS(folders,'spearman');
%     rate_remapping_TRACK_PAIRS(folders,'wcorr');
%     
    
    plot_rate_remapping('TRACK_PAIRS','spearman');
    plot_rate_remapping('TRACK_PAIRS','wcorr');
  

end


end