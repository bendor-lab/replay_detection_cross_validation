function [] = ground_truth_log_odd_analysis_spike_train_randomised_addon(folders,timebin,BAYSESIAN_NORMALIZED_ACROSS_TRACKS)
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

for f = 1:length(folders)
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
    %% spike train circular shift data (just redo)
    for shuffle =1:3
        cd spike_train_shuffles
        shuffle_folder = sprintf('shuffle_%i',shuffle);
        cd(shuffle_folder)
        load replayEvents_shifted_bayesian_spike_count % Also save in shuffle folder
        load scored_replay
        load decoded_replay_events
        cd ..
        cd ..

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
                shuffle_type{shuffle_id}.shuffled_track = randomised_data_parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,'spike_train_shifted',num_shuffles,...
                    decoded_replay_events,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);
            end
        else
            disp('parallel processing not possible');
            for shuffle_id=1:length(shuffle_choice)
                shuffle_type{shuffle_id}.shuffled_track = randomised_data_run_shuffles(shuffle_choice{shuffle_id},analysis_type,'spike_train_shifted',num_shuffles,...
                    decoded_replay_events,place_fields_BAYESIAN,BAYSESIAN_NORMALIZED_ACROSS_TRACKS);
            end
        end

        cd(sprintf('spike_train_shuffles/%s',shuffle_folder))
        save shuffled_tracks shuffle_type;

        % Evaluate significance
        scored_replay = replay_significance(scored_replay, shuffle_type);
        save scored_replay scored_replay
        clear scored_replay decoded_replay_events shuffle_type estimated_sequence_spike_train_shifted

    cd ..
    cd ..
    end
    cd(current_directory)
end

