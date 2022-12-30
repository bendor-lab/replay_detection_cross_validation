
function [p_values, replay_scores] = extract_score_and_pvalue_ground_truth(scored_replay, scored_replay1, scored_replay2, method, replay_above_rippleThresh_index)
% INPUTS:
%     Scored_replay, scored_replay1, scored_replay2: structures loaded. Contain score and pvalue of each replay event for each method analysed.
%     Method: method of analysis. 'path' for path finding/pacman, 'wcorr' for weighted correlation; 'linear' for linear fit and 'spearman' for spearman correlation 

number_of_tracks = length(scored_replay);
number_of_events = length(scored_replay(1).replay_events);

if strcmp(method,'path')   %path finding, pac-man method
    for track = 1 :  number_of_tracks % for each track
        for event = 1 : number_of_events %for each event
            p_values.WHOLE(track,event,:) = scored_replay(track).replay_events(event).p_value_path;
            replay_scores.WHOLE(track,event) = scored_replay(track).replay_events(event).path_score;
            p_values.FIRST_HALF(track,event,:) = scored_replay1(track).replay_events(event).p_value_path;
            replay_scores.FIRST_HALF(track,event) = scored_replay1(track).replay_events(event).path_score;
            p_values.SECOND_HALF(track,event,:) = scored_replay2(track).replay_events(event).p_value_path;
            replay_scores.SECOND_HALF(track,event) = scored_replay2(track).replay_events(event).path_score;
        end
    end
    
elseif strcmp(method,'wcorr') %weighted correlation method
    for track = 1 :  number_of_tracks % for each track
        for event = 1 : number_of_events %for each event
            p_values.WHOLE(track,event,:) = scored_replay(track).replay_events(event).p_value_wcorr; % pvalue whole event for the three shuffles
            replay_scores.WHOLE(track,event) = scored_replay(track).replay_events(event).weighted_corr_score; % wcorr score whole event
            p_values.FIRST_HALF(track,event,:) = scored_replay1(track).replay_events(event).p_value_wcorr; % pvalue first half of event
            replay_scores.FIRST_HALF(track,event) = scored_replay1(track).replay_events(event).weighted_corr_score; % wcorr score first half of event for the three shuffles
            p_values.SECOND_HALF(track,event,:) = scored_replay2(track).replay_events(event).p_value_wcorr; % pvalue second half of event
            replay_scores.SECOND_HALF(track,event) = scored_replay2(track).replay_events(event).weighted_corr_score; % wcorr score second half of event for the three shuffles
        end
    end
    
elseif strcmp(method,'linear') %linear fit  method
    for track = 1 :  number_of_tracks % for each track
        for event = 1 : number_of_events %for each event
            p_values.WHOLE(track,event,:) = scored_replay(track).replay_events(event).p_value_linear;
            replay_scores.WHOLE(track,event) = scored_replay(track).replay_events(event).linear_score;
            p_values.FIRST_HALF(track,event,:) = scored_replay1(track).replay_events(event).p_value_linear;
            replay_scores.FIRST_HALF(track,event) = scored_replay1(track).replay_events(event).linear_score;
            p_values.SECOND_HALF(track,event,:) = scored_replay2(track).replay_events(event).p_value_linear;
            replay_scores.SECOND_HALF(track,event) = scored_replay2(track).replay_events(event).linear_score;
        end
    end
    
elseif strcmp(method,'spearman') %spearman correlation coefficient
    for track = 1 :  number_of_tracks % for each track
        for event = 1 : number_of_events %for each event
            p_values.WHOLE(track,event,:) = scored_replay(track).replay_events(event).spearman_p*[1 1 1];
            replay_scores.WHOLE(track,event) = scored_replay(track).replay_events(event).spearman_score;
            p_values.FIRST_HALF(track,event,:) = scored_replay1(track).replay_events(event).spearman_p*[1 1 1];
            replay_scores.FIRST_HALF(track,event) = scored_replay1(track).replay_events(event).spearman_score;
            p_values.SECOND_HALF(track,event,:) = scored_replay2(track).replay_events(event).spearman_p*[1 1 1];
            replay_scores.SECOND_HALF(track,event) = scored_replay2(track).replay_events(event).spearman_score;
        end
    end

elseif strcmp(method,'spearman_all_spikes') %spearman correlation coefficient (All spikes)
    for track = 1 :  number_of_tracks % for each track
        for event = 1 : number_of_events %for each event
            p_values.WHOLE(track,event,:) = scored_replay(track).replay_events(event).spearman_all_p*[1 1 1];
            replay_scores.WHOLE(track,event) = scored_replay(track).replay_events(event).spearman_all_score;
            p_values.FIRST_HALF(track,event,:) = scored_replay1(track).replay_events(event).spearman_all_p*[1 1 1];
            replay_scores.FIRST_HALF(track,event) = scored_replay1(track).replay_events(event).spearman_all_score;
            p_values.SECOND_HALF(track,event,:) = scored_replay2(track).replay_events(event).spearman_all_p*[1 1 1];
            replay_scores.SECOND_HALF(track,event) = scored_replay2(track).replay_events(event).spearman_all_score;
        end
    end       
    
else disp('ERROR: scoring method not recognized')
    
end

%remove NANs for first and second half events - convert pvalue to 1 (to make it non significant) and score to 0 (minimum score)
p_values.FIRST_HALF(find(isnan(p_values.FIRST_HALF))) = 1;
replay_scores.FIRST_HALF(find(isnan(replay_scores.FIRST_HALF))) = 0;
p_values.SECOND_HALF(find(isnan(p_values.SECOND_HALF))) = 1;
replay_scores.SECOND_HALF(find(isnan(replay_scores.SECOND_HALF))) = 0;

% If we have applied the threshold for ripple power, then select the correct indices
if ~isempty(replay_above_rippleThresh_index)
    p_values.WHOLE = p_values.WHOLE(:,replay_above_rippleThresh_index,:);
    p_values.FIRST_HALF = p_values.FIRST_HALF(:,replay_above_rippleThresh_index,:);
    p_values.SECOND_HALF = p_values.SECOND_HALF(:,replay_above_rippleThresh_index,:);
    replay_scores.WHOLE = replay_scores.WHOLE(:,replay_above_rippleThresh_index);
    replay_scores.FIRST_HALF = replay_scores.FIRST_HALF(:,replay_above_rippleThresh_index);
    replay_scores.SECOND_HALF = replay_scores.SECOND_HALF(:,replay_above_rippleThresh_index);
end

end