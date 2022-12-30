function [] = ground_truth_replay_sequence_analysis(folders,BAYSESIAN_NORMALIZED_ACROSS_TRACKS)
for nfolder = 2:2
    tic
    cd(folders{nfolder})
    load extracted_place_fields_BAYESIAN

    % EXTRACT REPLAY EVENTS and BAYESIAN DECODING
    %     disp('processing replay events')
    replay_decoding_ground_truth(BAYSESIAN_NORMALIZED_ACROSS_TRACKS); %extract and decodes replay events
    % SCORING METHODS: TEST SIGNIFICANCE ON REPLAY EVENTS

    disp('scoring replay events')
    scored_replay = replay_scoring([],[1 1 1 1]);
    save scored_replay scored_replay

    %     % RUN SHUFFLES
    %     disp('running shuffles')
    num_shuffles=1000;
    analysis_type=[1 1 1 0];  %[linear wcorr path spearman]
    load decoded_replay_events
    p = gcp; % Starting new parallel pool
    shuffle_choice={'PRE spike_train_circular_shift','PRE place_field_circular_shift', 'POST place bin circular shift','POST time bin circular shift','POST time bin permutation'};
    
    if ~isempty(p)
        for shuffle_id=1:length(shuffle_choice)
            shuffle_type{shuffle_id}.shuffled_track = ground_truth_parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,...
                decoded_replay_events,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);

%             shuffle_type{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events);

        end
    else
        disp('parallel processing not possible');
        for shuffle_id=1:length(shuffle_choice)
                        shuffle_type{shuffle_id}.shuffled_track = ground_truth_run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,...
                            decoded_replay_events,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);

%             shuffle_type{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events);
        end

    end
    %
    %     tempt = shuffle_type;
    %     load shuffled_tracks
    %     for type = 1:length(shuffle_choice)
    %         for track = 1:length(shuffle_type{type}.shuffled_track)
    %             for event = 1:length(shuffle_type{type}.shuffled_track(track).replay_events)
    %                 shuffle_type{type}.shuffled_track(track).replay_events(event).path_score = ...
    %                     tempt{type}.shuffled_track(track).replay_events(event).path_score;
    %
    %             end
    %         end
    %     end
    %
    save shuffled_tracks shuffle_type;
    %     tempt = [];


    % Evaluate significance
    load('scored_replay.mat');
    load('shuffled_tracks.mat');
    scored_replay=replay_significance(scored_replay, shuffle_type);
    save scored_replay scored_replay

    %%%%%analyze segments%%%%%%%%%%
    % splitting replay events
    replay_decoding_split_events;
    load decoded_replay_events_segments;
    scored_replay1 = replay_scoring(decoded_replay_events1,[1 1 1 1]);
    scored_replay2 = replay_scoring(decoded_replay_events2,[1 1 1 1]);
    %
    %     load scored_replay_segments
    %     for track = 1:length(scored_replay)
    %         for events = 1:length(scored_replay(1).replay_events)
    %             scored_replay1(track).replay_events(events).path_score = tempt1(track).replay_events(events).path_score;
    %             scored_replay2(track).replay_events(events).path_score = tempt2(track).replay_events(events).path_score;
    %         end
    %     end

    save scored_replay_segments scored_replay1 scored_replay2;


    num_shuffles=1000;
    analysis_type=[1 1 1 0];  %just weighted correlation and pacman
    load decoded_replay_events_segments;
    p = gcp; % Starting new parallel pool

    if ~isempty(p)
        for shuffle_id=1:length(shuffle_choice)
            shuffle_type1{shuffle_id}.shuffled_track = ground_truth_parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,...
                decoded_replay_events1,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);
            shuffle_type2{shuffle_id}.shuffled_track = ground_truth_parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,...
                decoded_replay_events2,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);
            %
            %             shuffle_type1{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events1);
            %             shuffle_type2{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events2);

        end
    else
        disp('parallel processing not possible');
        for shuffle_id=1:length(shuffle_choice)
            shuffle_type1{shuffle_id}.shuffled_track = ground_truth_run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,...
                decoded_replay_events1,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);
            shuffle_type2{shuffle_id}.shuffled_track = ground_truth_run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,...
                decoded_replay_events2,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);

            %             shuffle_type1{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events1);
            %             shuffle_type2{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events2);

        end
    end

    % add specific score to scored_replay
    %     tempt1 = shuffle_type1;
    %     tempt2 = shuffle_type2;
    %
    %     load shuffled_tracks_segments
    %     for type = 1:length(shuffle_choice)
    %         for track = 1:length(shuffle_type{type}.shuffled_track)
    %             for event = 1:length(shuffle_type{type}.shuffled_track(track).replay_events)
    %                 shuffle_type1{type}.shuffled_track(track).replay_events(event).path_score = ...
    %                     tempt1{type}.shuffled_track(track).replay_events(event).path_score;
    %                 shuffle_type2{type}.shuffled_track(track).replay_events(event).path_score = ...
    %                     tempt2{type}.shuffled_track(track).replay_events(event).path_score;
    %             end
    %         end
    %     end
    %     tempt1 = [];
    %     tempt2 = [];

    save shuffled_tracks_segments shuffle_type1 shuffle_type2;
    load scored_replay_segments; load shuffled_tracks_segments;
    scored_replay1=replay_significance(scored_replay1, shuffle_type1);
    scored_replay2=replay_significance(scored_replay2, shuffle_type2);
    save scored_replay_segments scored_replay1 scored_replay2

    cd ..
    toc
end
end