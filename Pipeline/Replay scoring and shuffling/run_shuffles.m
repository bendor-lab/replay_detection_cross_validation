
function shuffled_track = run_shuffles(shuffle_choice,analysis_type,num_shuffles,decoded_replay_events)
% computes shuffles on replay events
% shuffle_choice- determines what type of shuffles will be performed, entered as a string
%
% analysis type- vector of 0's and 1's for which replaying
% scoring method is used [line fitting, weighted
%correlation, "pac-man" path finding, and spearman correlation coefficient]
% (e.g. [1 0 0 1] would be line fitting and spearman correlation)
%
% num_shuffles= number of shuffles
%decoded_replay_events is data obtained from "load decoded_replay_events"
analysis_type(4)=0;  %don't do shuffles on spearman correlation

if isempty(decoded_replay_events)
    load decoded_replay_events;
end
if strcmp(shuffle_choice, 'PRE spike_train_circular_shift') | strcmp(shuffle_choice, 'PRE place_field_circular_shift')
    load('extracted_place_fields_BAYESIAN.mat');
end
if strcmp(shuffle_choice, 'PRE spike_train_circular_shift')
    load('replayEvents_bayesian_spike_count.mat');
    replay_indices = replayEvents_bayesian_spike_count.replay_events_indices;
    spike_count_structure = replayEvents_bayesian_spike_count;
end

num_tracks = length(decoded_replay_events);
num_replay_events = length(decoded_replay_events(1).replay_events);

for s=1:num_shuffles
    shuffled_struct = [];
    %%%%%%%%
    %%%POST place bin circular shift
    %%%%%%%
    if strcmp(shuffle_choice, 'POST place bin circular shift') %
        % For each shuffle, creates a new structure where the position rows within each time bin in each decoded event have been shuffled
        for track = 1 : num_tracks
            for event = 1: num_replay_events
                pixel_size = size(decoded_replay_events(track).replay_events(event).decoded_position);
                shuffled_struct(track).replay_events(event)=decoded_replay_events(track).replay_events(event);
                index=decoded_replay_events(track).replay_events(event).timebins_index;  %index is needed when analyzing segments of replay event
                if length(index)>=5
                    for k = 1 : pixel_size(2)
                        shuffled_struct(track).replay_events(event).decoded_position(:,k) = circshift(shuffled_struct(track).replay_events(event).decoded_position(:,k),ceil(rand*pixel_size(1)));
                    end
                else
                    shuffled_struct(track).replay_events(event).decoded_position=NaN;
                end
            end
        end
    %%%%%%%%
    %%%POST time bin circular shift
    %%%%%%%        
    elseif strcmp(shuffle_choice, 'POST time bin circular shift') %
        % For each shuffle, creates a new structure where the time bins within each position bin in each decoded event have been shuffled
        for track = 1 : num_tracks
            for event = 1: num_replay_events
                pixel_size = size(decoded_replay_events(track).replay_events(event).decoded_position);
                shuffled_struct(track).replay_events(event)=decoded_replay_events(track).replay_events(event);
                index=decoded_replay_events(1).replay_events(event).timebins_index;  %index is needed when analyzing segments of replay event, same value across all tracks
                shuffled_index = circshift(index,ceil(rand*pixel_size(2)));

                if length(index)>=5
                    shuffled_struct(track).replay_events(event).decoded_position = decoded_replay_events(track).replay_events(event).decoded_position(:,shuffled_index);
                else
                    shuffled_struct(track).replay_events(event).decoded_position=NaN;
                end
            end
        end        
        

        %%%%%%%%
        %%%POST time bin permutation
        %%%%%%%
    elseif strcmp(shuffle_choice, 'POST time bin permutation') %
        % For each shuffle, creates a new structure where the time bins within each position bin in each decoded event have been shuffled
        for track = 1 : num_tracks
            for event = 1: num_replay_events
                pixel_size = size(decoded_replay_events(track).replay_events(event).decoded_position);
                shuffled_struct(track).replay_events(event)=decoded_replay_events(track).replay_events(event);
                index=decoded_replay_events(1).replay_events(event).timebins_index;  %index is needed when analyzing segments of replay event, same value across all tracks
                index = index-index(1)+1; % Start with 1
                shuffled_index = index(randperm(length(index)));
                

                if length(index)>=5
                    shuffled_struct(track).replay_events(event).decoded_position = decoded_replay_events(track).replay_events(event).decoded_position(:,shuffled_index);
                else
                    shuffled_struct(track).replay_events(event).decoded_position=NaN;
                end
            end
        end


        %%%%%%%%
        %%%PRE spike_train_circular_shift
        %%%%%%%
    elseif strcmp(shuffle_choice, 'PRE spike_train_circular_shift')
        replay_events_spike_count = spike_count_structure.n.replay;
        % For each replay event, takes spike count matrix and does circular shuffle on the spikes of each cell
        for event = 1 : length(replayEvents_bayesian_spike_count.replay_events)
            thisReplay_indxs = find(replay_indices == event);
            index=decoded_replay_events(1).replay_events(event).timebins_index;  %index is needed when analyzing segments of replay event, same value across all tracks
            for i=1:size(replay_events_spike_count,1)
                %thisReplay_indxs(index) is used below to handle whole replay events and segments, for whole events- thisReplay_indxs(index) == thisReplay_indxs
                replay_events_spike_count(i,thisReplay_indxs(index)) = circshift(replay_events_spike_count(i,thisReplay_indxs(index)),ceil(rand*length(index)),2);
            end
        end
        
        % Save in structure to input in bayesian decoding code
        shuffled_spike_count.n.replay = replay_events_spike_count;
        shuffled_spike_count.replay_events_indices = replay_indices;
        shuffled_spike_count.replay_events = replayEvents_bayesian_spike_count.replay_events;
        
        % Each replay event is then decoded using the new spike count structure
        estimated_position = bayesian_decoding([],shuffled_spike_count,'N');
        
        for track = 1:num_tracks
            for event = 1 : num_replay_events
                shuffled_struct(track).replay_events(event)=decoded_replay_events(track).replay_events(event);
                index=decoded_replay_events(track).replay_events(event).timebins_index;  %index is needed when analyzing segments of replay event
                if length(index)>=5
                    shuffled_struct(track).replay_events(event).decoded_position = estimated_position(track).replay_events(event).replay(:,index); % normalized by all tracks
                else
                    shuffled_struct(track).replay_events(event).decoded_position = NaN;
                end
            end
        end
        
        %%%%%%%%
        %%%PRE place_field_circular_shift
        %%%%%%%
        
    elseif strcmp(shuffle_choice, 'PRE place_field_circular_shift')
        shuffled_place_fields=place_fields_BAYESIAN;
        for track = 1 : num_tracks
            num_cells = size(shuffled_place_fields.track(track).raw,2);
            %circular shift place field bins
            for i = 1 : num_cells
                field =  cell2mat(shuffled_place_fields.track(track).raw(i));
                shuffled_place_fields.track(track).raw(i) = {circshift(field, ceil(rand*length(field)))};
            end
        end
        estimated_position = bayesian_decoding(shuffled_place_fields,'replayEvents_bayesian_spike_count','N');
        %calculate decoded position with shuffled place fields
        for track = 1:num_tracks
            for event = 1 : num_replay_events
                shuffled_struct(track).replay_events(event)=decoded_replay_events(track).replay_events(event);
                index=decoded_replay_events(track).replay_events(event).timebins_index;  %index is needed when analyzing segments of replay event
                if length(index)>=5
                    shuffled_struct(track).replay_events(event).decoded_position = estimated_position(track).replay_events(event).replay(:,index); % normalize probability across all tracks
                else
                    shuffled_struct(track).replay_events(event).decoded_position = NaN;
                end
            end
        end
    end
    %score replay for decoded events with shuffled place fields
    shuffle_output=replay_scoring(shuffled_struct,analysis_type);  %don't do shuffle for spearman
    for track = 1 : num_tracks
        for event = 1: num_replay_events
            shuffled_track(track).replay_events(event).linear_score(s) = shuffle_output(track).replay_events(event).linear_score;
            shuffled_track(track).replay_events(event).weighted_corr_score(s) = shuffle_output(track).replay_events(event).weighted_corr_score;
            shuffled_track(track).replay_events(event).path_score(s) = shuffle_output(track).replay_events(event).path_score;
        end
        
    end
end
%save shuffle choice
for track = 1 : num_tracks
    shuffled_track(track).shuffle_choice=shuffle_choice;
end
end