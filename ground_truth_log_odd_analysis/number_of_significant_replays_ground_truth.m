function [significant_replay_events sig_event_info] = number_of_significant_replays_ground_truth(p_value_threshold,ripple_zscore_threshold, method, shuffles,rexposure_option)
% Establish significance for the replay events that are above a set ripple power threshold. 
% The output will therefore have less replay events than the structures that have been loaded. The output will be two variables:
    % sig_event_info - contains information about details of replay event
    % significant_replay_event - info about replay events that are significant on one track
% INPUTS
%     p_value_threshold: INT. Default is 0.05
%     ripple_zscore_threshold: INT. Default is 3
%     method: method of analysis. 'path' for path finding/pacman, 'wcorr' for weighted correlation; 'linear' for linear fit and 'spearman' for spearman correlation 
%     rexposure_option: 


% Load and set variables
load scored_replay_segments
load scored_replay
load extracted_replay_events
load extracted_position
load extracted_sleep_state
load decoded_replay_events
load decoded_replay_events_segments

parameters = list_of_parameters;
significant_replay_events.bin_size = 1;
significant_replay_events.smoothing_bin_size = 61*5;
significant_replay_events.bin_size_BIG = 61;

if isempty(ripple_zscore_threshold)
    ripple_zscore_threshold = 3; % 3 SD
end
if isempty(p_value_threshold)
    p_value_threshold = 0.05;
end

if ~isempty(rexposure_option)
    significant_replay_events.track1_and_1R_significant = [];
    significant_replay_events.track2_and_2R_significant = [];
end

% Find indices of replay events with ripple power above threshold
replay_above_rippleThresh_index = find(replay.ripple_peak >= ripple_zscore_threshold);

%extract replay score (e.g. wcorr value) and p value for each shuffle and event segment (whole, first, second) and for each track
[p_values, replay_scores] = extract_score_and_pvalue(scored_replay, scored_replay1, scored_replay2, method, replay_above_rippleThresh_index);
number_of_tracks = length(scored_replay);
number_of_events = length(replay_above_rippleThresh_index);

%create variable significant_replay_event that contains all information related to replay events judged significant (for each track) with this function
significant_replay_events.p_value_threshold = p_value_threshold;
significant_replay_events.ripple_zscore_threshold = ripple_zscore_threshold;
significant_replay_events.method = method;
significant_replay_events.pre_ripple_threshold_index = replay_above_rippleThresh_index;  %index of events with signficant ripple power (to use if referencing other .mat files)

% Find midpoint time for each event
significant_replay_events.all_event_times = (replay.onset(replay_above_rippleThresh_index)+replay.offset(replay_above_rippleThresh_index))/2;

% Calculate amount of replay events for each track based on scoring method
significant_replay_events.time_bin_edges = min(position.t):significant_replay_events.bin_size:max(position.t);  %five second bins
significant_replay_events.time_bin_centres = (min(position.t)+(significant_replay_events.bin_size/2)):significant_replay_events.bin_size:(max(position.t)-(significant_replay_events.bin_size/2));
significant_replay_events.HIST = histcounts(significant_replay_events.all_event_times,significant_replay_events.time_bin_edges);

significant_replay_events.time_bin_edges_BIG = min(position.t):significant_replay_events.bin_size_BIG:max(position.t);  %one minute bins
significant_replay_events.time_bin_centres_BIG = (min(position.t)+(significant_replay_events.bin_size_BIG/2)):significant_replay_events.bin_size_BIG:(max(position.t)-(significant_replay_events.bin_size_BIG/2));


% Check if the event against the 3 shuffles is significant
for track = 1 :  number_of_tracks % for each track
    
    for event = 1 : number_of_events %for each event
        
        if shuffles == 3
            % For each event and each scoring method, check if the highest p value within the 3 shuffles is significant (check for the entire whole and for both segments)
            sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,1:3))<p_value_threshold | max(p_values.FIRST_HALF(track,event,1:3))<p_value_threshold/2 | max(p_values.SECOND_HALF(track,event,1:3))<p_value_threshold/2; % indices - 0 or 1 for significant event
            [sig_event_info.p_value(track,event) sig_event_info.segment_id(track,event)] = min([max(p_values.WHOLE(track,event,1:3)) max(p_values.FIRST_HALF(track,event,1:3)) max(p_values.SECOND_HALF(track,event,1:3))]); %min pvalue between the 3 max pvalues (not the same than the min pvalue)

            %             sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,1:3) < p_value_threshold) ; % indices - 0 or 1 for significant event
            %             sig_event_info.p_value(track,event) = max(p_values.WHOLE(track,event,1:3)) ; %min pvalue between the 3 max pvalues (not the same than the min pvalue)
        
        elseif shuffles == 1 % 1 shufffle (POST place bin shuffle)
            % If only POST place bin circular shift
            sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,3))<p_value_threshold | max(p_values.FIRST_HALF(track,event,3))<p_value_threshold/2 | max(p_values.SECOND_HALF(track,event,3))<p_value_threshold/2; % indices - 0 or 1 for significant event
            [sig_event_info.p_value(track,event) sig_event_info.segment_id(track,event)] = min([max(p_values.WHOLE(track,event,3)) max(p_values.FIRST_HALF(track,event,3)) max(p_values.SECOND_HALF(track,event,3))]); %min pvalue between the 3 max pvalues (not the same than the min pvalue)

            %             sig_event_info.significance(track,event) = p_values.WHOLE(track,event,3)<p_value_threshold; % indices - 0 or 1 for significant event
            %             sig_event_info.p_value(track,event) = max(p_values.WHOLE(track,event,3)); %min pvalue between the 3 max pvalues (not the same than the min pvalue)

        elseif strcmp(shuffles,'PRE place and POST place')
            % if Only PRE place field circular shift and POST place bin
            % shift
            sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,2:3))<p_value_threshold | max(p_values.FIRST_HALF(track,event,2:3))<p_value_threshold/2 | max(p_values.SECOND_HALF(track,event,2:3))<p_value_threshold/2; % indices - 0 or 1 for significant event
            [sig_event_info.p_value(track,event) sig_event_info.segment_id(track,event)] = min([max(p_values.WHOLE(track,event,2:3)) max(p_values.FIRST_HALF(track,event,2:3)) max(p_values.SECOND_HALF(track,event,2:3))]); %min pvalue between the 3 max pvalues (not the same than the min pvalue)
            %             sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,2:3) < p_value_threshold); % indices - 0 or 1 for significant event
            %             sig_event_info.p_value(track,event) = max(p_values.WHOLE(track,event,2:3)) ; %min pvalue between the 3 max pvalues (not the same than the min pvalue)

        elseif strcmp(shuffles,'PRE place and POST time')
            % if Only PRE place field circular shift and POST place bin
            % shift
            sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,[2 5]))<p_value_threshold | max(p_values.FIRST_HALF(track,event,[2 5]))<p_value_threshold/2 | max(p_values.SECOND_HALF(track,event,[2 5]))<p_value_threshold/2; % indices - 0 or 1 for significant event
            [sig_event_info.p_value(track,event) sig_event_info.segment_id(track,event)] = min([max(p_values.WHOLE(track,event,[2 5])) max(p_values.FIRST_HALF(track,event,[2 5])) max(p_values.SECOND_HALF(track,event,[2 5]))]); %min pvalue between the 3 max pvalues (not the same than the min pvalue)
            %             sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,2:3) < p_value_threshold); % indices - 0 or 1 for significant event
            %             sig_event_info.p_value(track,event) = max(p_values.WHOLE(track,event,2:3)) ; %min pvalue between the 3 max pvalues (not the same than the min pvalue)

        elseif strcmp(shuffles,'PRE spike and POST place')
            % If only PRE spike train circular shift and PRE place bin
            % circular shift
            sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,[1 3]))<p_value_threshold | max(p_values.FIRST_HALF(track,event,[1 3]))<p_value_threshold/2 | max(p_values.SECOND_HALF(track,event,[1 3]))<p_value_threshold/2; % indices - 0 or 1 for significant event
            [sig_event_info.p_value(track,event) sig_event_info.segment_id(track,event)] = min([max(p_values.WHOLE(track,event,[1 3])) max(p_values.FIRST_HALF(track,event,[1 3])) max(p_values.SECOND_HALF(track,event,[1 3]))]); %min pvalue between the 3 max pvalues (not the same than the min pvalue)
            %             sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,3:4))<p_value_threshold; % indices - 0 or 1 for significant event
            %             sig_event_info.p_value(track,event) = max(p_values.WHOLE(track,event,3:4)); %min pvalue between the 3 max pvalues (not the same than the min pvalue)

        elseif strcmp(shuffles,'two POST')
            % If only POST place bin circular shift and POST time bin
            % circular shift
            sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,[3 5]))<p_value_threshold | max(p_values.FIRST_HALF(track,event,[3 5]))<p_value_threshold/2 | max(p_values.SECOND_HALF(track,event,[3 5]))<p_value_threshold/2; % indices - 0 or 1 for significant event
            [sig_event_info.p_value(track,event) sig_event_info.segment_id(track,event)] = min([max(p_values.WHOLE(track,event,[3 5])) max(p_values.FIRST_HALF(track,event,[3 5])) max(p_values.SECOND_HALF(track,event,[3 5]))]); %min pvalue between the 3 max pvalues (not the same than the min pvalue)
%             sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,3:4))<p_value_threshold; % indices - 0 or 1 for significant event
%             sig_event_info.p_value(track,event) = max(p_values.WHOLE(track,event,3:4)); %min pvalue between the 3 max pvalues (not the same than the min pvalue)
        
        elseif strcmp(shuffles,'two PRE')
            % If only PRE place field circular shift and PRE spike train
            % circular shift
            sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,1:2))<p_value_threshold | max(p_values.FIRST_HALF(track,event,1:2))<p_value_threshold/2 | max(p_values.SECOND_HALF(track,event,1:2))<p_value_threshold/2; % indices - 0 or 1 for significant event
            [sig_event_info.p_value(track,event) sig_event_info.segment_id(track,event)] = min([max(p_values.WHOLE(track,event,1:2)) max(p_values.FIRST_HALF(track,event,1:2)) max(p_values.SECOND_HALF(track,event,1:2))]); %min pvalue between the 3 max pvalues (not the same than the min pvalue)            
%   
%             sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,1:2))<p_value_threshold; % indices - 0 or 1 for significant event
%             sig_event_info.p_value(track,event) = max(p_values.WHOLE(track,event,1:2)); %min pvalue between the 3 max pvalues (not the same than the min pvalue)
                    
        elseif strcmp(shuffles,'PRE spike_train_circular_shift')
           % If only PRE spike train circular shift
            sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,1))<p_value_threshold | max(p_values.FIRST_HALF(track,event,1))<p_value_threshold/2 | max(p_values.SECOND_HALF(track,event,1))<p_value_threshold/2; % indices - 0 or 1 for significant event
            [sig_event_info.p_value(track,event) sig_event_info.segment_id(track,event)] = min([max(p_values.WHOLE(track,event,1)) max(p_values.FIRST_HALF(track,event,1)) max(p_values.SECOND_HALF(track,event,1))]); %min pvalue between the 3 max pvalues (not the same than the min pvalue)
%             sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,1))<p_value_threshold; % indices - 0 or 1 for significant event
%             sig_event_info.p_value(track,event) = max(p_values.WHOLE(track,event,1)); %min pvalue between the 3 max pvalues (not the same than the min pvalue)
                    
        elseif strcmp(shuffles,'PRE place_field_circular_shift') 
           % If only PRE place field circular shift
            sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,2))<p_value_threshold | max(p_values.FIRST_HALF(track,event,2))<p_value_threshold/2 | max(p_values.SECOND_HALF(track,event,2))<p_value_threshold/2; % indices - 0 or 1 for significant event
            [sig_event_info.p_value(track,event) sig_event_info.segment_id(track,event)] = min([max(p_values.WHOLE(track,event,2)) max(p_values.FIRST_HALF(track,event,2)) max(p_values.SECOND_HALF(track,event,2))]); %min pvalue between the 3 max pvalues (not the same than the min pvalue)
           
%             sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,2))<p_value_threshold; % indices - 0 or 1 for significant event
%             sig_event_info.p_value(track,event) = max(p_values.WHOLE(track,event,2)); %min pvalue between the 3 max pvalues (not the same than the min pvalue)
                    
        elseif strcmp(shuffles,'POST place bin circular shift')
            % If only POST place bin circular shift
            sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,3))<p_value_threshold | max(p_values.FIRST_HALF(track,event,3))<p_value_threshold/2 | max(p_values.SECOND_HALF(track,event,3))<p_value_threshold/2; % indices - 0 or 1 for significant event
            [sig_event_info.p_value(track,event) sig_event_info.segment_id(track,event)] = min([max(p_values.WHOLE(track,event,3)) max(p_values.FIRST_HALF(track,event,3)) max(p_values.SECOND_HALF(track,event,3))]); %min pvalue between the 3 max pvalues (not the same than the min pvalue)

%             sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,3))<p_value_threshold; % indices - 0 or 1 for significant event
%             sig_event_info.p_value(track,event) = max(p_values.WHOLE(track,event,3)); %min pvalue between the 3 max pvalues (not the same than the min pvalue)
%                                 
        elseif strcmp(shuffles,'POST time bin circular shift')
            sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,4))<p_value_threshold | max(p_values.FIRST_HALF(track,event,4))<p_value_threshold/2 | max(p_values.SECOND_HALF(track,event,4))<p_value_threshold/2; % indices - 0 or 1 for significant event
            [sig_event_info.p_value(track,event) sig_event_info.segment_id(track,event)] = min([max(p_values.WHOLE(track,event,4)) max(p_values.FIRST_HALF(track,event,4)) max(p_values.SECOND_HALF(track,event,4))]); %min pvalue between the 3 max pvalues (not the same than the min pvalue)
%             sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,4))<p_value_threshold; % indices - 0 or 1 for significant event
%             sig_event_info.p_value(track,event) = max(p_values.WHOLE(track,event,4)); %min pvalue between the 3 max pvalues (not the same than the min pvalue)
%          

        elseif strcmp(shuffles,'POST time bin permutation')
            sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,5))<p_value_threshold | max(p_values.FIRST_HALF(track,event,5))<p_value_threshold/2 | max(p_values.SECOND_HALF(track,event,5))<p_value_threshold/2; % indices - 0 or 1 for significant event
            [sig_event_info.p_value(track,event) sig_event_info.segment_id(track,event)] = min([max(p_values.WHOLE(track,event,5)) max(p_values.FIRST_HALF(track,event,5)) max(p_values.SECOND_HALF(track,event,5))]); %min pvalue between the 3 max pvalues (not the same than the min pvalue)
%             sig_event_info.significance(track,event) = max(p_values.WHOLE(track,event,4))<p_value_threshold; % indices - 0 or 1 for significant event
%             sig_event_info.p_value(track,event) = max(p_values.WHOLE(track,event,4)); %min pvalue between the 3 max pvalues (not the same than the min pvalue)
%               
        end
        
        
        for i = 1 : size(p_values.WHOLE,3) %checks if event is significant with a specific type of shuffle (for either WHOLE and/or event segments)
            sig_event_info.method_idx(track,event,i) = [(p_values.WHOLE(track,event,i))<p_value_threshold | (p_values.FIRST_HALF(track,event,i))<p_value_threshold/2 | (p_values.SECOND_HALF(track,event,i))<p_value_threshold/2];
        end
        
        % Is the highest p value within the 3 shuffles less than the p value threshold (for the whole event and for segments)
        sig_event_info.segment_idx(track,event,:) = [max(p_values.WHOLE(track,event,:))<p_value_threshold  max(p_values.FIRST_HALF(track,event,:))<p_value_threshold/2  max(p_values.SECOND_HALF(track,event,:))<p_value_threshold/2];
        
        % Save the score within segments and whole event for each scoring method
        sig_event_info.replay_scores(track,event,:) = [replay_scores.WHOLE(track,event) replay_scores.FIRST_HALF(track,event) replay_scores.SECOND_HALF(track,event)];
        
        %sig_event_info.replay_segment_index   --    0 if not signficant, 1 if whole event is significant, 2 if first half, 3 if second half
        %sig_event_info.best_replay_score  --   will be highest replay score on a significant segment (WHOLE or first/second half of event)
        if sig_event_info.significance(track,event) == 0
            sig_event_info.replay_segment_index(track,event) = 0;
            sig_event_info.best_replay_score(track,event) = 0;
        else
            %find maximum score from significant whole/segmented events
            [sig_event_info.best_replay_score(track,event),sig_event_info.replay_segment_index(track,event)] = max(sig_event_info.segment_idx(track,event,:).*sig_event_info.replay_scores(track,event,:));
        end
        
    end
end


% Find events that are significant for more than one track (multi_tracks doesn't consider events significant for first and second exposure, as this will be assigned later on)
if isempty(rexposure_option)
    multi_tracks_index = find(sum(sig_event_info.significance ,1)>1); %is more than one track significant
elseif rexposure_option == 1
    multi_tracks_index = [find(sum(sig_event_info.significance(1:2,:) ,1)>1) find(sum(sig_event_info.significance(3:4,:) ,1)>1) find(sum(sig_event_info.significance([1 4],:) ,1)>1) ...
        find(sum(sig_event_info.significance(2:3,:) ,1)>1)];    %check if overlap of tracks in initial exposure or rexposure
    multi_tracks_index = unique(multi_tracks_index);
elseif rexposure_option == 2  %you will only compare first exposure or second, but not cross comparisions.  so you can allow track A1 and track B2 to both be significant without calling it a multievent
    multi_tracks_index = [find(sum(sig_event_info.significance(1:2,:) ,1)>1) find(sum(sig_event_info.significance(3:4,:) ,1)>1)];    %check if overlap of tracks in initial exposure or rexposure
    multi_tracks_index = unique(multi_tracks_index);
end

% If only two significant tracks with one being the first segment and the second being the second segment, and no track has higher significant score for the whole replay event,
% this is an exception and both are simultaneously signficant (no need to remove) - WHOLE/SEGMENT based on highest significant score
% exceptions = [];
% for i = 1:length(multi_tracks_index)
%     if (length(find(sig_event_info.replay_segment_index(:,multi_tracks_index(i))==1))==0 & ...
%             length(find(sig_event_info.replay_segment_index(:,multi_tracks_index(i))==2))==1 & ...
%             length(find(sig_event_info.replay_segment_index(:,multi_tracks_index(i))==3))==1)
%         exceptions = [exceptions i];
%     end
% end

% significant_replay_events.multi_track_BUT_diff_segments_exception_index = multi_tracks_index(exceptions);
% multi_tracks_index(exceptions) = [];
significant_replay_events.multi_tracks_index = multi_tracks_index;

% Take list of significant replay events, and deal with events where more than one track is significant
sig_event_info.significance_NO_MULTI = sig_event_info.significance;
sig_event_info.significance_NO_MULTI(:,multi_tracks_index) = 0; % events significant just for one track from the start
sig_event_info.significance_MULTI = zeros(size(sig_event_info.significance));
sig_event_info.significance_MULTI(:,multi_tracks_index) = sig_event_info.significance(:,multi_tracks_index);


%HISTOGRAMS OF REPLAY ACTIVITY 
for track=1:number_of_tracks
    sig_event_info.track(track).HIST = histcounts(significant_replay_events.all_event_times(sig_event_info.significance(track,:)==1),significant_replay_events.time_bin_edges); %significant events for each track at each time bin (includes multiple events)
%     sig_event_info.track(track).HIST_BEST_BAYESIAN = histcounts(significant_replay_events.all_event_times(sig_event_info.significance_BEST_BAYESIAN(track,:)==1),significant_replay_events.time_bin_edges);%significant events for each track at each time bin (events ONLY sig for one track)
    sig_event_info.track(track).HIST_MULTI = histcounts(significant_replay_events.all_event_times(sig_event_info.significance_MULTI(track,:)==1),significant_replay_events.time_bin_edges); % events sig for multiple tracks at each time bin
    sig_event_info.track(track).HIST_NO_MULTI = histcounts(significant_replay_events.all_event_times(sig_event_info.significance_NO_MULTI(track,:)==1),significant_replay_events.time_bin_edges); % events sig for one track only at each time bin
end


%REPLAY INFORMATION - output of function
for track = 1 :  number_of_tracks
    significant_replay_events.track(track).HIST = sig_event_info.track(track).HIST_NO_MULTI;% Instead of best bayesian, we use no multi
    significant_replay_events.track(track).index = find(sig_event_info.significance(track,:)==1); %index of sig events for this track
    significant_replay_events.track(track).ref_index = replay_above_rippleThresh_index(sig_event_info.significance(track,:)==1);
    significant_replay_events.track(track).ref_index_NO_MULTI = replay_above_rippleThresh_index(sig_event_info.significance_NO_MULTI(track,:)==1);
    significant_replay_events.track(track).ref_index_MULTI = replay_above_rippleThresh_index(sig_event_info.significance_MULTI(track,:)==1);
    
    significant_replay_events.track(track).event_times = significant_replay_events.all_event_times(significant_replay_events.track(track).index);
    significant_replay_events.track(track).replay_score = sig_event_info.best_replay_score(track,significant_replay_events.track(track).index); %score of sig replay events
    significant_replay_events.track(track).p_value =  sig_event_info.p_value(track,significant_replay_events.track(track).index); %pvalue
    significant_replay_events.track(track).event_segment_best_score = sig_event_info.replay_segment_index(track,significant_replay_events.track(track).index); % info about if events is sig for whole or segment
        
    for i = 1:length(significant_replay_events.track(track).index)
        if (significant_replay_events.track(track).event_segment_best_score(i)==1)  %if best score in WHOLE event
            significant_replay_events.track(track).spikes{i} = decoded_replay_events(track).replay_events(replay_above_rippleThresh_index(significant_replay_events.track(track).index(i))).spikes;
            significant_replay_events.track(track).decoded_position{i} = decoded_replay_events(track).replay_events(replay_above_rippleThresh_index(significant_replay_events.track(track).index(i))).decoded_position;
            significant_replay_events.track(track).event_duration(i) = range(decoded_replay_events(track).replay_events(replay_above_rippleThresh_index(significant_replay_events.track(track).index(i))).timebins_edges);
       
        elseif (significant_replay_events.track(track).event_segment_best_score(i)==2) %if best score in first half of event
            significant_replay_events.track(track).spikes{i} = decoded_replay_events1(track).replay_events(replay_above_rippleThresh_index(significant_replay_events.track(track).index(i))).spikes;
            significant_replay_events.track(track).decoded_position{i} = decoded_replay_events1(track).replay_events(replay_above_rippleThresh_index(significant_replay_events.track(track).index(i))).decoded_position;
            significant_replay_events.track(track).event_duration(i) = range(decoded_replay_events1(track).replay_events(replay_above_rippleThresh_index(significant_replay_events.track(track).index(i))).timebins_edges);
        
        elseif (significant_replay_events.track(track).event_segment_best_score(i)==3) %if best score in second half of event
            significant_replay_events.track(track).spikes{i} = decoded_replay_events2(track).replay_events(replay_above_rippleThresh_index(significant_replay_events.track(track).index(i))).spikes;
            significant_replay_events.track(track).decoded_position{i} = decoded_replay_events2(track).replay_events(replay_above_rippleThresh_index(significant_replay_events.track(track).index(i))).decoded_position;
            significant_replay_events.track(track).event_duration(i) = range(decoded_replay_events2(track).replay_events(replay_above_rippleThresh_index(significant_replay_events.track(track).index(i))).timebins_edges);
        else
            significant_replay_events.track(track).spikes{i} = NaN;
            significant_replay_events.track(track).decoded_position{i} = NaN;
            significant_replay_events.track(track).event_duration(i) = NaN;
        end
    end
end


% For plotting
for track = 1 : number_of_tracks
    for j = 1 : length(significant_replay_events.time_bin_edges_BIG)-11 %for each time bin
        index = find(significant_replay_events.track(track).event_times >= significant_replay_events.time_bin_edges_BIG(j) & significant_replay_events.track(track).event_times < significant_replay_events.time_bin_edges_BIG(j+11));
        if isempty(index)
            significant_replay_events.track(track).mean_best_score(j) = NaN;
            significant_replay_events.track(track).mean_bayesian_bias(j) = NaN;
            significant_replay_events.track(track).median_log_pvalue(j) = NaN;
            significant_replay_events.track(track).min_log_pvalue(j) = NaN;
        else
            significant_replay_events.track(track).mean_best_score(j) = mean(significant_replay_events.track(track).replay_score(index));
%             significant_replay_events.track(track).mean_bayesian_bias(j) = mean(significant_replay_events.track(track).bayesian_bias(index));
            significant_replay_events.track(track).median_log_pvalue(j) = median(log10(1e-5+significant_replay_events.track(track).p_value(index)));  %limit p-value to 1e-4
            significant_replay_events.track(track).min_log_pvalue(j) = min(log10(1e-5+significant_replay_events.track(track).p_value(index)));  %limit p-value to 1e-4
            
        end
    end
end

% Excluded events because they are significant for more than one track and can't be tell a part
% number_of_excluded_events_from_BAYESIAN_BIAS = length(significant_replay_events.BAYESIAN_BIAS_excluded_events_index)
% disp(number_of_excluded_events_from_BAYESIAN_BIAS)

% Save
switch method
    case 'wcorr'
        if ~isempty(rexposure_option) && rexposure_option == 2
            save significant_replay_events_wcorr_individual_exposures significant_replay_events
        else
            save significant_replay_events_wcorr significant_replay_events
        end
    case 'spearman'
        if ~isempty(rexposure_option) && rexposure_option == 2
            save significant_replay_events_spearman_individual_exposures significant_replay_events
        else
            save significant_replay_events_spearman significant_replay_events
        end
    case 'spearman_all_spikes'
        if ~isempty(rexposure_option) && rexposure_option == 2
            save significant_replay_events_spearman_all_spikes_individual_exposures significant_replay_events
        else
            save significant_replay_events_spearman_all_spikes significant_replay_events
        end        
    case 'path'
        if ~isempty(rexposure_option) && rexposure_option == 2
            save significant_replay_events_path_individual_exposures significant_replay_events
        else
            save significant_replay_events_path significant_replay_events
        end        
    case 'linear'
        if ~isempty(rexposure_option) && rexposure_option == 2
            save significant_replay_events_linear_individual_exposures significant_replay_events
        else
            save significant_replay_events_linear significant_replay_events
        end        
end

end


function [p_values, replay_scores] = extract_score_and_pvalue(scored_replay, scored_replay1, scored_replay2, method, replay_above_rippleThresh_index)
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
    
elseif strcmp(method,'spearman_all_spikes') %spearman correlation coefficient
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