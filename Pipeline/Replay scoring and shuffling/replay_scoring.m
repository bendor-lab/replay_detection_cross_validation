
function scored_replay = replay_scoring(data,analysis_type)
%scores replay events using several measures- line fitting, weighted
%correlation, "pac-man" path finding, and spearman correlation coefficient.
%modify_decoding.m subroutine is run to zero probability if no spikes from
%good place fields occur during time bin (which creates unwanted noise in
%these scoring algorithms
%only analyzes replay if there are no NaNs and at least 5 time bins in duration
%
%
% data- can be empty of you want to use "decoded_replay_events.mat"
%, or pass in the decoded replay of shuffle events
%
% analysis type- vector of 0's and 1's for which replaying
% scoring method is used [line fitting, weighted
%correlation, "pac-man" path finding, and spearman correlation coefficient]
% (e.g. [1 0 0 1] would be line fitting and spearman correlation)

if isempty(data)
    load('decoded_replay_events.mat');
    data=decoded_replay_events;
    clear decoded_replay_events;
end
num_tracks = length(data);
num_replay_events = length(data(1).replay_events);

for track = 1 : num_tracks
    for event = 1 : num_replay_events
        scored_replay(track).replay_events(event).linear_score=NaN;
        scored_replay(track).replay_events(event).weighted_corr_score=NaN;
        scored_replay(track).replay_events(event).spearman_score=NaN;
        scored_replay(track).replay_events(event).spearman_p=NaN;
        scored_replay(track).replay_events(event).path_score=NaN;
    end
end


% Line fitting
if analysis_type(1) == 1
    for track = 1 : num_tracks
        for event = 1 : num_replay_events
            decoded_event = data(track).replay_events(event);
            if length(find(isnan(decoded_event.decoded_position)))>0 | length(decoded_event.timebins_centre)<5 | size(decoded_event.spikes,1)==0 %NaN or event too short
                scored_replay(track).replay_events(event).linear_score=NaN;
            else
                decoded_position=modify_decoding(decoded_event);  % make column zero if there aren't any spikes
                [scored_replay(track).replay_events(event).linear_score,~] = line_fitting(decoded_position);
            end
        end
    end
end


% Weighted correlation
if analysis_type(2) == 1
    for track = 1 : num_tracks
        for event = 1 : num_replay_events
            decoded_event = data(track).replay_events(event);
            if length(find(isnan(decoded_event.decoded_position)))>0 | length(decoded_event.timebins_centre)<5 | size(decoded_event.spikes,1)==0 %NaN or event too short
                scored_replay(track).replay_events(event).weighted_corr_score=NaN;
            else
                decoded_position=modify_decoding(decoded_event);  % make column zero if there aren't any spikes
                scored_replay(track).replay_events(event).weighted_corr_score = weighted_correlation(decoded_position);
            end
        end
    end
end
% "pacman" path finding
if analysis_type(3) == 1
    for track = 1 : num_tracks
        for event = 1 : num_replay_events
            decoded_event = data(track).replay_events(event);
            if length(find(isnan(decoded_event.decoded_position)))>0 | length(decoded_event.timebins_centre)<5 | size(decoded_event.spikes,1)==0
                scored_replay(track).replay_events(event).path_score=NaN;
            else
                decoded_position=modify_decoding(decoded_event);  % make column zero if there aren't any spikes
%                 [scored_replay(track).replay_events(event).path_score,~] = pacman(decoded_position);
                [scored_replay(track).replay_events(event).path_score] = fast_pacman(decoded_position);
            end
        end
    end
    
end
% spearman correlation
if analysis_type(4) == 1
    load extracted_place_fields_BAYESIAN
    for track = 1 : num_tracks
        sorted_place_fields=place_fields_BAYESIAN.track(track).sorted_good_cells;
        decoded_event = data(track).replay_events(event);
        for event = 1 : num_replay_events
            if length(find(isnan(decoded_event.decoded_position)))>0 | length(decoded_event.timebins_centre)<5 | size(decoded_event.spikes,1)==0
                scored_replay(track).replay_events(event).spearman_score=NaN;
                scored_replay(track).replay_events(event).spearman_p=1;
            else
                spike_id=data(track).replay_events(event).spikes(:,1);
                spike_times=data(track).replay_events(event).spikes(:,2);
                [scored_replay(track).replay_events(event).spearman_score,scored_replay(track).replay_events(event).spearman_p] = spearman_median(spike_id, spike_times, sorted_place_fields);
            end
        end
    end
end
end


function modified_decoded_event = modify_decoding(events)
modified_decoded_event = events.decoded_position;
for i = 1 : length(events.timebins_edges)-1
    spikes = find(events.spikes(:,2)>=events.timebins_edges(i) & events.spikes(:,2)<events.timebins_edges(i+1),1);  %>= and < used, to match histcount function
    if isempty(spikes)
        modified_decoded_event(:,i) = zeros(size(modified_decoded_event(:,i)));
    end
end
end