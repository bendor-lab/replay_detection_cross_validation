%calculate whether the events are significant, based on shuffles

function scored_replay=replay_significance(scored_replay, shuffle_type)  
%extracts the replay event's p value by comparing the replay score with the distrbution of shuffled scores
%
%scored_replay is a variable of replay scores, 
%shuffle_type is a variable containing the replay scores for each shuffle type
%
%

num_shuffle_types = length(shuffle_type);
num_tracks = length(scored_replay);
num_replay_events = length(scored_replay(1).replay_events);
%extract replay scores for each replay event, across all tracks, across all methods
    for track = 1: num_tracks
        for event = 1 : num_replay_events
            for s=1:num_shuffle_types
                scored_replay(track).replay_events(event).p_value_wcorr(s)=...
                    get_p_value(scored_replay(track).replay_events(event).weighted_corr_score,...
                    shuffle_type{s}.shuffled_track(track).replay_events(event).weighted_corr_score);
                scored_replay(track).replay_events(event).p_value_linear(s)=...
                    get_p_value(scored_replay(track).replay_events(event).linear_score,...
                    shuffle_type{s}.shuffled_track(track).replay_events(event).linear_score);
                scored_replay(track).replay_events(event).p_value_path(s)=...
                    get_p_value(scored_replay(track).replay_events(event).path_score,...
                    shuffle_type{s}.shuffled_track(track).replay_events(event).path_score);
            end
        end
    end
end 


function out = get_p_value(event_score,shuffles_scores)
% Finds the proportion of scores in the shuffled distribution greater than the score of the candidate trajectory
out = length(find(shuffles_scores>=event_score))/length(shuffles_scores); 
if ~isempty(find(isnan(shuffles_scores),1)) %if shuffle scores are NaNs
    out = NaN;
end
end