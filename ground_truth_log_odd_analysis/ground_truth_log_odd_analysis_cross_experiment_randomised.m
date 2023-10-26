function [] = ground_truth_log_odd_analysis_cross_experiment_randomised(folders,timebin,BAYSESIAN_NORMALIZED_ACROSS_TRACKS)
% Input:
% Folders: Each Folder is data for a session
% Timebin: 0.02 or 1
% 0.02 -> 20ms time bin
% 1 -> 1 event time bin
% posbin: [] or 1
% [] -> 20 position bin per track (will load extracted_place_field_BAYESIAN)
% 1 -> 1 position bin per track

% Output:
% log_odd: save the log odd - log(sum(probability of T1)/sum(probability of
% T2)) of the events that are considered significant according to the original bayesian decoding
% as well as the 1000 T1&T2 turning curve shuffles for each event
% It will also contain the session ID, 'ground truth' replay track ID and
% behavioral epoch (-1: PRE, 0: POST, 1: T1 RUN, 2: T2 RUN) and etc

c = 1;
current_directory=pwd;
load subsets_of_cells;
load place_fields_BAYESIAN_combined;

for f = [6 7]
    cd Tables
    load subsets_of_cells;
    cd ..

    cd(folders{f})
    parameters = list_of_parameters;
    place_fields_BAYESIAN = [];
    load('significant_replay_events_wcorr.mat');
    load('sorted_replay_wcorr');
    load('extracted_place_fields_BAYESIAN');
    load('extracted_position')
    load('extracted_clusters.mat');
    load('extracted_replay_events.mat');
    place_fields_BAYESIAN.track = rmfield(place_fields_BAYESIAN.track,{'global_remapped','global_remapped_raw','rate_fixed_raw','rate_fixed_scaling_factor'})
    
    % Stable common good cells (based on subsets_of_cells)
    %     place_cell_index = subset_of_cells.cell_IDs{8}...
    %         (intersect(find(subset_of_cells.cell_IDs{8} >= 1000*f),find(subset_of_cells.cell_IDs{8} <= 1000*(f+1))))- f*1000;
    place_cell_index = place_fields_BAYESIAN.good_place_cells

    % if place_field_index is not initialised, it affects parallel loop for
    % some reasons
    place_field_index = [];

    %% Decoded Track Replay Structure
    % Extracts spikes for each replay event and decodes it. Saves in a
    % structure called repay_track. Inside each field (each track) the replay events are saved.
    % Loads: extracted_replay_events, extracted_clusters and extracted_place_fields_BAYESIAN

    % REPLAY EVENTS STRUCTURE
    %replay_events is an empty template for replay event analysis.
    %Each track will create its own field

    replay_events = struct('replay_id',{},...%the id of the candidate replay events in chronological order
        'spikes',{});%column 1 is spike id, column 2 is spike time

    % TAKE SPIKES FROM ONLY common good cells on both tracks

    sorted_spikes = zeros(size(clusters.spike_id));
    sorted_spikes(:,1) = clusters.spike_id;
    sorted_spikes(:,2) = clusters.spike_times;
    all_units = unique(clusters.spike_id);

    non_good = setdiff(all_units,place_cell_index);% Find cells that are not good place cells

    for i = 1 : length(non_good)
        non_good_indices = find(sorted_spikes(:,1)== non_good(i));
        sorted_spikes(non_good_indices,:) = [];% remove spikes from unwanted cells
        non_good_indices =[];
    end

    num_spikes = length(sorted_spikes);
    num_units = length(place_cell_index);

    % EXTRACT SPIKES IN REPLAY EVENTS
    num_replay = size(replay.onset, 2);
    current_replay = 1;
    current_replay_spikes = [];

    for i = 1 : num_spikes
        % Collect spike data during replay
        if sorted_spikes(i,2) > replay.offset(current_replay)
            replay_events(current_replay).replay_id = current_replay;
            replay_events(current_replay).spikes = current_replay_spikes;
            current_replay = current_replay + 1;
            if current_replay > num_replay
                break
            end
            current_replay_spikes = [];
        end

        if sorted_spikes(i,2) >= replay.onset(current_replay)
            % If spike happens during replay, records it as replay spike
            current_replay_spikes = [current_replay_spikes; sorted_spikes(i,:)];
        end

        if i == num_spikes && current_replay == num_replay && isempty(current_replay_spikes)
            % for the last replay event, if it is empty, write empty
            replay_events(current_replay).replay_id = current_replay;
            replay_events(current_replay).spikes = current_replay_spikes;
        end

    end

    num_replay_events = length(replay_events);
    msg = [num2str(num_replay_events), ' candidate events.'];
    disp(msg);
    
    load replayEvents_bayesian_spike_count

    %% cross experiment cell id shuffled data
    for shuffle = 1:3

        bayesian_spike_count = 'replayEvents_bayesian_spike_count';
        cross_experiment_shuffled_place_fields_id = [];

        for event = 1:length(replayEvents_bayesian_spike_count.replay_events)
         
            % cross experiment cell id shuffles
            cross_experiment_shuffled_place_fields = [];


            for track_id = 1:2
                cross_experiment_shuffled_place_fields{track_id} = place_fields_BAYESIAN.track(track_id).raw;

                % Find cells that are not from this experiment
                random_cell_id = datasample(find(place_fields_BAYESIAN_combined.track(track_id).session_id~=f), ...
                    length(place_fields_BAYESIAN.track(track_id).sorted_good_cells), 'Replace', false);

                % Original cell id
                original_cell_id = place_fields_BAYESIAN.track(track_id).sorted_good_cells;

                [~,index1] = sort(place_fields_BAYESIAN_combined.track(track_id).centre(random_cell_id));
                %                  [~,index2] = sort(place_fields_BAYESIAN.track(track_id).centre(original_cell_id));
                [~,index2] = sort(index1);

                for j=1:length(original_cell_id) %only swap good cells
                    cross_experiment_shuffled_place_fields{track_id}{original_cell_id(j)}=place_fields_BAYESIAN_combined.track(track_id).raw{random_cell_id(j)};
                end

                place_fields_BAYESIAN.track(track_id).cross_experiment_shuffled{event} = cross_experiment_shuffled_place_fields{track_id};

                % Also save the shuffled place cell id for subsequent spearsman
                % analysis
                cross_experiment_shuffled_place_fields_id{event}{track_id}(1,:) = original_cell_id;
                cross_experiment_shuffled_place_fields_id{event}{track_id}(2,:) = original_cell_id(index2);

            end

        end


        cross_experiment_shuffled_original_probability_ratio = [];
        cross_experiment_shuffled_common_good_probability_ratio = [];
%         save extracted_place_fields_BAYESIAN place_fields_BAYESIAN

        % Sequence decoding
        estimated_sequence_cross_experiment_shuffled = bayesian_decoding_ground_truth(...
            place_fields_BAYESIAN,'replayEvents_bayesian_spike_count','cross_experiment_shuffled',BAYSESIAN_NORMALIZED_ACROSS_TRACKS,'N');
% 
%         estimated_sequence = bayesian_decoding_ground_truth(...
%             place_fields_BAYESIAN,'replayEvents_bayesian_spike_count',[],BAYSESIAN_NORMALIZED_ACROSS_TRACKS,'N');
        %  estimated_sequence_global_remapped = bayesian_decoding(place_fields_BAYESIAN,'replayEvents_bayesian_spike_count','N');

        if ~isfolder({'cross_experiment_shuffles'})
            mkdir('cross_experiment_shuffles')
        end

        cd cross_experiment_shuffles

        shuffle_folder = sprintf('shuffle_%i',shuffle);

        if ~isfolder(shuffle_folder)
            mkdir(shuffle_folder)
        end


        cd(shuffle_folder)
        save cross_experiment_shuffled_place_fields_id cross_experiment_shuffled_place_fields_id
        save extracted_place_fields_BAYESIAN place_fields_BAYESIAN
        cd ..
        cd ..

        for j = 1:length(place_fields_BAYESIAN.track)
            decoded_replay_events(j).replay_events = replay_events;
        end

        % Save in structure
        for j = 1:length(place_fields_BAYESIAN.track)
            for i = 1:length(estimated_sequence_cross_experiment_shuffled(1).replay_events)
                decoded_replay_events(j).replay_events(i).timebins_edges = estimated_sequence_cross_experiment_shuffled(j).replay_events(i).replay_time_edges;
                decoded_replay_events(j).replay_events(i).timebins_centre = estimated_sequence_cross_experiment_shuffled(j).replay_events(i).replay_time_centered;
                decoded_replay_events(j).replay_events(i).timebins_index = 1:length(estimated_sequence_cross_experiment_shuffled(j).replay_events(i).replay_time_centered);
                decoded_replay_events(j).replay_events(i).decoded_position = estimated_sequence_cross_experiment_shuffled(j).replay_events(i).replay; % normalized within one track
            end
        end

        scored_replay = replay_scoring(decoded_replay_events,[1 1 1 1]);

        % RUN SHUFFLES
        disp('running shuffles')
        num_shuffles=1000;
        analysis_type=[1 1 1 0];  %just linear fit, weighted correlation and pacman
        p = gcp; % Starting new parallel pool
        shuffle_choice={'PRE spike_train_circular_shift','PRE place_field_circular_shift', 'POST place bin circular shift',...
            'POST time bin circular shift','POST time bin permutation','PRE cell_id_shuffle'};
        place_field_index = [];

        if ~isempty(p)
            for shuffle_id=1:length(shuffle_choice)
                shuffle_type{shuffle_id}.shuffled_track = randomised_data_parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,'cross_experiment_shuffled',num_shuffles,...
                    decoded_replay_events,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);
            end
        else
            disp('parallel processing not possible');
            for shuffle_id=1:length(shuffle_choice)
                shuffle_type{shuffle_id}.shuffled_track = randomised_data_run_shuffles(shuffle_choice{shuffle_id},analysis_type,'cross_experiment_shuffled',num_shuffles,...
                    decoded_replay_events,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);
            end
        end

        cd cross_experiment_shuffles
        cd(shuffle_folder)
        save shuffled_tracks shuffle_type;
        save decoded_replay_events decoded_replay_events;

        % Evaluate significance
        scored_replay = replay_significance(scored_replay, shuffle_type);
        save scored_replay scored_replay
        cd ..
        cd ..
        clear scored_replay decoded_replay_events shuffle_type estimated_sequence_cross_experiment_shuffled

        destination = ['cross_experiment_shuffles\shuffle_' num2str(shuffle)];
        copyfile('extracted_replay_events.mat',destination)
        copyfile('extracted_position.mat',destination)
        copyfile('extracted_sleep_state.mat',destination)

        %%%%%%analyze segments%%%%%%%%%%
        %splitting replay events
        cd cross_experiment_shuffles
        cd(shuffle_folder)

        replay_decoding_split_events;
        load decoded_replay_events_segments;
        cd ..
        cd ..

        scored_replay1 = replay_scoring(decoded_replay_events1,[1 1 1 1]);
        scored_replay2 = replay_scoring(decoded_replay_events2,[1 1 1 1]);
        cd cross_experiment_shuffles
        cd(shuffle_folder)
        save scored_replay_segments scored_replay1 scored_replay2;

        num_shuffles=1000;
        analysis_type=[1 1 1 0];  %just weighted correlation and pacman
        shuffle_choice={'PRE spike_train_circular_shift','PRE place_field_circular_shift', 'POST place bin circular shift',...
            'POST time bin circular shift','POST time bin permutation','PRE cell_id_shuffle'};
        load decoded_replay_events_segments;
        cd ..
        cd ..


        p = gcp; % Starting new parallel pool
        if ~isempty(p)
            for shuffle_id=1:length(shuffle_choice)
                shuffle_type1{shuffle_id}.shuffled_track = randomised_data_parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,'cross_experiment_shuffled',num_shuffles,...
                    decoded_replay_events1,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);
                shuffle_type2{shuffle_id}.shuffled_track = randomised_data_parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,'cross_experiment_shuffled',num_shuffles,...
                    decoded_replay_events2,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);
                %                 shuffle_type1{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events1);
                %                 shuffle_type2{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events2);
            end
        else
            disp('parallel processing not possible');
            for shuffle_id=1:length(shuffle_choice)
                shuffle_type1{shuffle_id}.shuffled_track = randomised_data_run_shuffles(shuffle_choice{shuffle_id},analysis_type,'cross_experiment_shuffled',num_shuffles,...
                    decoded_replay_events1,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);
                shuffle_type2{shuffle_id}.shuffled_track = randomised_data_run_shuffles(shuffle_choice{shuffle_id},analysis_type,'cross_experiment_shuffled',num_shuffles,...
                    decoded_replay_events2,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);

                %                 shuffle_type1{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events1);
                %                 shuffle_type2{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events2);
            end
        end
        cd cross_experiment_shuffles
        cd(shuffle_folder)
        save shuffled_tracks_segments shuffle_type1 shuffle_type2;
        load scored_replay_segments;
        load shuffled_tracks_segments;
        scored_replay1=replay_significance(scored_replay1, shuffle_type1);
        scored_replay2=replay_significance(scored_replay2, shuffle_type2);
        save scored_replay_segments scored_replay1 scored_replay2

        cd ..
        cd ..


        if isfile('cross_experiment_shuffled_original_probability_ratio.mat')
            load cross_experiment_shuffled_original_probability_ratio;
        end

        % % % log odd analysis % % %

        % % % All Good Cells % % %

        %     shuffle = 1;

        estimated_position_cross_experiment_shuffled_original = [];
        % Calculate probability ratio (100 shuffles)
        place_cell_index = place_fields_BAYESIAN.good_place_cells;
        bayesian_spike_count = 'replayEvents_bayesian_spike_count';
        estimated_position_cross_experiment_shuffled_original = log_odd_bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count,place_cell_index,timebin,[],[],'cross_experiment_shuffled','N');
       
        for event = 1:length(estimated_position_cross_experiment_shuffled_original(1).replay_events)
            cross_experiment_shuffled_original_probability_ratio{shuffle}{1}(1,event,1) = estimated_position_cross_experiment_shuffled_original(1).replay_events(event).probability_ratio;
            cross_experiment_shuffled_original_probability_ratio{shuffle}{1}(2,event,1) = 1/(estimated_position_cross_experiment_shuffled_original(1).replay_events(event).probability_ratio);
            cross_experiment_shuffled_original_probability_ratio{shuffle}{1}(1,event,2) = estimated_position_cross_experiment_shuffled_original(1).replay_events(event).probability_ratio_first;
            cross_experiment_shuffled_original_probability_ratio{shuffle}{1}(2,event,2) = 1/(estimated_position_cross_experiment_shuffled_original(1).replay_events(event).probability_ratio_first);
            cross_experiment_shuffled_original_probability_ratio{shuffle}{1}(1,event,3) = estimated_position_cross_experiment_shuffled_original(1).replay_events(event).probability_ratio_second;
            cross_experiment_shuffled_original_probability_ratio{shuffle}{1}(2,event,3) = 1/(estimated_position_cross_experiment_shuffled_original(1).replay_events(event).probability_ratio_second);

            disp(sprintf('cross experiment shuffled Event : %i (shuffle %i)...',event,shuffle));
        end


        % Calculate the T1-T2 ratemap (turning curve) shuffle 1000 times
        % such that there are 1000 shuffles for each event.
        parfor nshuffle = 1:1000
            estimated_position_ratemap_shuffled = [];
            estimated_position_ratemap_shuffled = log_odd_bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count,place_cell_index,timebin,[],'ratemap shuffle','cross_experiment_shuffled','N');
            for event = 1:length(estimated_position_ratemap_shuffled(1).replay_events)
                ratemap_shuffled_probability_ratio{nshuffle}(1,event,1) = estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio;
                ratemap_shuffled_probability_ratio{nshuffle}(2,event,1) = 1/(estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio);
                ratemap_shuffled_probability_ratio{nshuffle}(1,event,2) = estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio_first;
                ratemap_shuffled_probability_ratio{nshuffle}(2,event,2) = 1/(estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio_first);
                ratemap_shuffled_probability_ratio{nshuffle}(1,event,3) = estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio_second;
                ratemap_shuffled_probability_ratio{nshuffle}(2,event,3) = 1/(estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio_second);
            end
            disp(sprintf('Ratemap shuffling: %i (shuffle %i)...',nshuffle,shuffle));
        end

        cross_experiment_shuffled_original_probability_ratio{shuffle}{2} = ratemap_shuffled_probability_ratio;

        save cross_experiment_shuffled_original_probability_ratio cross_experiment_shuffled_original_probability_ratio
        clear ratemap_shuffled_probability_ratio cross_experiment_shuffled_original_probability_ratio

        % % % Common Good Cells % % %

        if isfile('cross_experiment_shuffled_common_good_probability_ratio.mat')
            load cross_experiment_shuffled_common_good_probability_ratio;
        end

        bayesian_spike_count = 'replayEvents_common_good_cell_bayesian_spike_count';
        place_cell_index = subset_of_cells.cell_IDs{1}...
            (intersect(find(subset_of_cells.cell_IDs{1} >= 1000*f),find(subset_of_cells.cell_IDs{1} <= 1000*(f+1))))- f*1000;
        
        % log odd analysis
        estimated_position_cross_experiment_shuffled_common_good = [];
        % Calculate probability ratio (100 shuffles)
        estimated_position_cross_experiment_shuffled_common_good = log_odd_bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count,place_cell_index,timebin,[],[],'cross_experiment_shuffled','N');
        for event = 1:length(estimated_position_cross_experiment_shuffled_common_good(1).replay_events)
            cross_experiment_shuffled_common_good_probability_ratio{shuffle}{1}(1,event,1) = estimated_position_cross_experiment_shuffled_common_good(1).replay_events(event).probability_ratio;
            cross_experiment_shuffled_common_good_probability_ratio{shuffle}{1}(2,event,1) = 1/(estimated_position_cross_experiment_shuffled_common_good(1).replay_events(event).probability_ratio);
            cross_experiment_shuffled_common_good_probability_ratio{shuffle}{1}(1,event,2) = estimated_position_cross_experiment_shuffled_common_good(1).replay_events(event).probability_ratio_first;
            cross_experiment_shuffled_common_good_probability_ratio{shuffle}{1}(2,event,2) = 1/(estimated_position_cross_experiment_shuffled_common_good(1).replay_events(event).probability_ratio_first);
            cross_experiment_shuffled_common_good_probability_ratio{shuffle}{1}(1,event,3) = estimated_position_cross_experiment_shuffled_common_good(1).replay_events(event).probability_ratio_second;
            cross_experiment_shuffled_common_good_probability_ratio{shuffle}{1}(2,event,3) = 1/(estimated_position_cross_experiment_shuffled_common_good(1).replay_events(event).probability_ratio_second);



            disp(sprintf('Common cell cross experiment shuffled Event : %i (shuffle %i)...',event,shuffle));
        end


        % Calculate the T1-T2 ratemap (turning curve) shuffle 1000 times
        % such that there are 1000 shuffles for each event.
        parfor nshuffle = 1:1000
            estimated_position_ratemap_shuffled = [];
            estimated_position_ratemap_shuffled = log_odd_bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count,place_cell_index,timebin,[],'ratemap shuffle','cross_experiment_shuffled','N');
            for event = 1:length(estimated_position_ratemap_shuffled(1).replay_events)
                ratemap_shuffled_probability_ratio{nshuffle}(1,event,1) = estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio;
                ratemap_shuffled_probability_ratio{nshuffle}(2,event,1) = 1/(estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio);
                ratemap_shuffled_probability_ratio{nshuffle}(1,event,2) = estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio_first;
                ratemap_shuffled_probability_ratio{nshuffle}(2,event,2) = 1/(estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio_first);
                ratemap_shuffled_probability_ratio{nshuffle}(1,event,3) = estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio_second;
                ratemap_shuffled_probability_ratio{nshuffle}(2,event,3) = 1/(estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio_second);  
            end
             disp(sprintf('Common cell ratemap shuffling: %i (shuffle %i)...',nshuffle,shuffle));
        end

        cross_experiment_shuffled_common_good_probability_ratio{shuffle}{2} = ratemap_shuffled_probability_ratio;

        save cross_experiment_shuffled_common_good_probability_ratio cross_experiment_shuffled_common_good_probability_ratio

        clear  cross_experiment_shuffled_original_probability_ratio cross_experiment_shuffled_common_good_probability_ratio ratemap_shuffled_probability_ratio

    end
    
    cd(current_directory)
end

end





