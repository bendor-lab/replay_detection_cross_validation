function [] = ground_truth_replay_sequence_analysis_addon_global_remapped(folders,BAYSESIAN_NORMALIZED_ACROSS_TRACKS)

for nfolder = 8:length(folders)
    cd(workingDir)
    for nshuffle = 2:3
        tic
        cd([workingDir,'/',folders{nfolder}])
        copyfile('replayEvents_bayesian_spike_count.mat',...
            [workingDir,'\',folders{nfolder},sprintf('/global_remapped_shuffles/shuffle_%i/replayEvents_bayesian_spike_count.mat',nshuffle)]);
        cd('global_remapped_shuffles')
        cd(sprintf('shuffle_%i',nshuffle))
        
        load extracted_place_fields_BAYESIAN

%         %     % RUN SHUFFLES
%         %     disp('running shuffles')
%         num_shuffles=1000;
%         analysis_type=[1 1 1 0];  %[linear wcorr path spearman]
%         load decoded_replay_events
%         p = gcp; % Starting new parallel pool
%         shuffle_choice={'PRE cell_id_shuffle'};
% 
%         load shuffled_tracks
%         %     shuffle_type{4}.shuffled_track = parallel_shuffles(shuffle_choice,analysis_type,num_shuffles,decoded_replay_events);
%         %     shuffle_type{5}.shuffled_track = parallel_shuffles(shuffle_choice,analysis_type,num_shuffles,decoded_replay_events);
%         shuffle_type{6}.shuffled_track = global_remapped_parallel_shuffles(shuffle_choice,analysis_type,num_shuffles,...
%             decoded_replay_events,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);
%         save shuffled_tracks shuffle_type;
% 
%         % Evaluate significance
%         load('scored_replay.mat');
%         scored_replay=replay_significance(scored_replay, shuffle_type);
%         save scored_replay scored_replay

        %%%%%analyze segments%%%%%%%%%%
        % splitting replay events
        load decoded_replay_events_segments;
        load shuffled_tracks_segments;
        num_shuffles=1000;
        analysis_type=[1 1 1 0];  %just weighted correlation and pacman
        shuffle_choice={'PRE cell_id_shuffle'};
        p = gcp; % Starting new parallel pool
        %     shuffle_type1{4}.shuffled_track = parallel_shuffles(shuffle_choice,analysis_type,num_shuffles,decoded_replay_events1);
        %     shuffle_type2{4}.shuffled_track = parallel_shuffles(shuffle_choice,analysis_type,num_shuffles,decoded_replay_events2);

        %         shuffle_type1{6}.shuffled_track = parallel_shuffles(shuffle_choice,analysis_type,num_shuffles,decoded_replay_events1);

        %         shuffle_type1{6}.shuffled_track = parallel_shuffles(shuffle_choice,analysis_type,num_shuffles,decoded_replay_events1);
        shuffle_type1{6}.shuffled_track = global_remapped_parallel_shuffles(shuffle_choice,analysis_type,num_shuffles,...
            decoded_replay_events1,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);
        %         shuffle_type2{6}.shuffled_track = parallel_shuffles(shuffle_choice,analysis_type,num_shuffles,decoded_replay_events2);
        shuffle_type2{6}.shuffled_track = global_remapped_parallel_shuffles(shuffle_choice,analysis_type,num_shuffles,...
            decoded_replay_events2,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);

        save shuffled_tracks_segments shuffle_type1 shuffle_type2;

        load scored_replay_segments;
        scored_replay1=replay_significance(scored_replay1, shuffle_type1);
        scored_replay2=replay_significance(scored_replay2, shuffle_type2);
        save scored_replay_segments scored_replay1 scored_replay2

        toc
    end
end

end
