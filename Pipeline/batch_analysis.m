function batch_analysis(option)
rexposure_exception=[];  %if you have four tracks but they are not rexposures...
rexposure=[];  %if not empty, assumes that there is rexposure to track


if strcmp(option,'EXTRACT_SPIKES_MARTA')
    extract_spikes_phy; % for clusters extracted from PHY
elseif strcmp(option,'EXTRACT_SPIKES_MARGOT')
    extract_spikes_samples;
elseif strcmp(option,'EXTRACT_VIDEO')
    extract_video;
elseif strcmp(option,'PRE')
    process_positions_PRE; 
elseif strcmp(option,'POST')
    load position_data;
    number_of_tracks=length(position.linear);
    if number_of_tracks<4 | ~isempty(rexposure_exception)
        rexposure=[];
    else
        rexposure=[1,3; 2,4];
    end
    extract_events;
    extract_dropped_samples;
    process_clusters;
    process_positions_POST(rexposure);  %input is empty unless you have track rexposure
    sleep_stager('auto');  %process sleep stages with default thresholds

elseif strcmp(option,'CSC')         % EXTRACT CSC
    disp('processing CSC data')
    p = gcp; % Starting new parallel pool
    if ~isempty(p)
        parallel_extract_PSD('sleep');
    else
        disp('parallel processing not possible');
        extract_PSD([],'sleep');
    end
    % plot_PSDs;
    best_channels = determine_best_channel('hpc');
    extract_CSC('hpc');

elseif strcmp(option,'PLACE FIELDS')       
    % EXTRACT PLACE FIELDS
    disp('processing place_field data')
    parameters=list_of_parameters;
    calculate_place_fields(parameters.x_bins_width_bayesian);
    calculate_place_fields(parameters.x_bins_width);
    spike_count([],[],[],'Y');
    bayesian_decoding([],[],'Y');
    extract_replay_events; %finds onset and offset of replay events    

elseif strcmp(option,'REPLAY')
    % EXTRACT REPLAY EVENTS and BAYESIAN DECODING
    disp('processing replay events')
    replay_decoding; %extract and decodes replay events
    % SCORING METHODS: TEST SIGNIFICANCE ON REPLAY EVENTS
    disp('scoring replay events')
    scored_replay = replay_scoring([],[1 1 1 1]);
    save scored_replay scored_replay;
    
    % RUN SHUFFLES
    disp('running shuffles')
    num_shuffles=1000;
    analysis_type=[1 1 1 0];  %just linear fit, weighted correlation and pacman
    load decoded_replay_events
    p = gcp; % Starting new parallel pool
    shuffle_choice={'PRE spike_train_circular_shift','PRE place_field_circular_shift', 'POST place bin circular shift'};
    
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
    
    % Evaluate significance
    load('scored_replay.mat');
    load('shuffled_tracks.mat');
    scored_replay=replay_significance(scored_replay, shuffle_type);
    save scored_replay scored_replay
    
    %%%%%%analyze segments%%%%%%%%%%
    %splitting replay events
    replay_decoding_split_events;
    load decoded_replay_events_segments;
    scored_replay1 = replay_scoring(decoded_replay_events1,[1 1 1 1]);
    scored_replay2 = replay_scoring(decoded_replay_events2,[1 1 1 1]);
    save scored_replay_segments scored_replay1 scored_replay2;
    num_shuffles=1000;
    analysis_type=[1 1 1 0];  %just weighted correlation and pacman
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

elseif strcmp(option,'SLEEP')
    % EXTRACT SLEEP
    sleep_stager('manual');

elseif strcmp(option,'ANALYSIS')
    current_directory=pwd;
    load scored_replay;
    number_of_tracks=length(scored_replay);
    if number_of_tracks<4 | ~isempty(rexposure_exception)
        rexposure=[];
    else
        rexposure=[1,3; 2,4];
    end
    significant_replay_events=number_of_significant_replays(0.05,3,'wcorr',rexposure);
    [sorted_replay_index,time_range]=sort_replay_events(rexposure);
    cumulative_replay([]);
    cd(current_directory);
elseif strcmp(option,'compare methods')
    current_directory=pwd;
    load scored_replay;
    number_of_tracks=length(scored_replay);
    if number_of_tracks<4 | ~isempty(rexposure_exception)
        rexposure=[];
    else
        rexposure=[1,3; 2,4];
    end
    significant_replay_events=number_of_significant_replays(0.05,3,'linear',rexposure);
    plot_significant_events;
    significant_replay_events=number_of_significant_replays(0.05,3,'wcorr',rexposure);
    plot_significant_events;
    significant_replay_events=number_of_significant_replays(0.05,3,'path',rexposure);
    plot_significant_events;
    significant_replay_events=number_of_significant_replays(0.05,3,'spearman',rexposure);
    plot_significant_events;
    cd(current_directory);
elseif strcmp(option,'plot_cleaning_steps')
    plot_cleaning_steps;
elseif strcmp(option,'plot_place_fields')
    plot_place_fields;
elseif strcmp(option,'plot_significant_events')
    plot_significant_events;
elseif strcmp(option,'plot_cumulative_replay')
    plot_cumulative_replay;
elseif strcmp(option,'plot_rate_remap')
    current_directory=pwd;
    %out=rate_remapping_ONE_TRACK([]);
    %plot_rate_remapping('ONE_TRACK');
    cd(current_directory)
    out=rate_remapping_TRACK_PAIRS([]);
    plot_rate_remapping('TRACK_PAIRS');
    cd(current_directory)
end