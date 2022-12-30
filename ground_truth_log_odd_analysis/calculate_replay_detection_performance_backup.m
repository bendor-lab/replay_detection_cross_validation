function out = calculate_replay_detection_performance_backup(data,epoch_index,segment_id,log_pval,p_val_threshold,session_index,event_index,total_number,option)
%                 multi_event_index = find(log_odd_compare{nshuffle}.multi_track == 1);
out = [];
log_odd_difference = [];
percent_multi_events = [];
percent_sig_events = [];
if strcmp(option,'bootstrap')
    parfor threshold = 1:1
        % find the event index with p value lower than the threshold
        %             current_index = find(log_pval{nmethod}(index{nmethod}{epoch})< p_val_threshold(threshold));
        difference = [];
        percent = [];
        track_2_index = [];
        track_1_index = [];

        for nboot = 1:1000 % Bootstrapping 1000 times
            resampled_event = datasample([epoch_index{1} epoch_index{2}],...
                length([epoch_index{1} epoch_index{2}]));

            this_segment_event = resampled_event(ismember(resampled_event,find(segment_id == 1)));
            thresholded_event = this_segment_event(find(log_pval(this_segment_event) <= p_val_threshold(threshold)));
            this_segment_event = resampled_event(ismember(resampled_event,find(segment_id > 1)));
            thresholded_event = [thresholded_event this_segment_event(find(log_pval(this_segment_event) <= p_val_threshold(threshold)/2))];

            %                         thresholded_event = p_value_thresholding(resampled_event,segment_id{nmethod}{nshuffle},log_pval{nmethod}{nshuffle},p_val_threshold(threshold));

            track_1_index = thresholded_event(ismember(thresholded_event,epoch_index{1}));
            track_2_index = thresholded_event(ismember(thresholded_event,epoch_index{2}));


            log_odd_difference(threshold,nboot) = mean(data(track_1_index)) ...
                - mean(data(track_2_index));

            %                         [multi_event_percent(nboot) multi_event_number(nboot)] = calculate_multitrack_event_percentage(session_index{nmethod}{nshuffle},event_index{nmethod}{nshuffle},track_1_index,track_2_index,epoch)
            [multi_event_percent multi_event_number] = calculate_multitrack_event_percentage(session_index,event_index,track_1_index,track_2_index);
            percent_multi_events(threshold,nboot) = multi_event_percent;
            % Significant event proportion (minus multitrack event number to avoid double counting)
            percent_sig_events(threshold,nboot) = (length(track_1_index)+length(track_2_index)-multi_event_number)/total_number;
        end
    end

elseif strcmp(option,'original')
    
    for threshold = 1:length(p_val_threshold)
        all_events = [epoch_index{1} epoch_index{2}];
        this_segment_event = all_events(ismember(all_events,find(segment_id == 1)));
        thresholded_event = this_segment_event(find(log_pval(this_segment_event) <= p_val_threshold(threshold)));
        this_segment_event = all_events(ismember(all_events,find(segment_id > 1)));
        thresholded_event = [thresholded_event this_segment_event(find(log_pval(this_segment_event) <= p_val_threshold(threshold)/2))];

        %                         thresholded_event = p_value_thresholding(resampled_event,segment_id{nmethod}{nshuffle},log_pval{nmethod}{nshuffle},p_val_threshold(threshold));

        track_1_index = thresholded_event(ismember(thresholded_event,epoch_index{1}));
        track_2_index = thresholded_event(ismember(thresholded_event,epoch_index{2}));


        log_odd_difference(threshold) = mean(data(track_1_index)) ...
            - mean(data(track_2_index));
        %                         [multi_event_percent(nboot) multi_event_number(nboot)] = calculate_multitrack_event_percentage(session_index{nmethod}{nshuffle},event_index{nmethod}{nshuffle},track_1_index,track_2_index,epoch)
        [multi_event_percent multi_event_number] = calculate_multitrack_event_percentage(session_index,event_index,track_1_index,track_2_index);
        percent_multi_events(threshold) = multi_event_percent;
        % Significant event proportion (minus multitrack event number to avoid double counting)
        percent_sig_events(threshold) = (length(track_1_index)+length(track_2_index)-multi_event_number)/total_number;
    end
end


out.log_odd_difference = log_odd_difference;
% out.log_odd_difference_CI= log_odd_difference_CI;
out.percent_sig_events = percent_sig_events;
% out.percent_sig_events_CI = percent_sig_events_CI;
out.percent_multi_events = percent_multi_events;
% out.percent_multi_events_CI = percent_multi_events_CI;
end
% 

for nboot = 1:1 % Bootstrapping 1000 times
%     resampled_event = datasample([epoch_index{nmethod}{nshuffle}{epoch}{1} epoch_index{nmethod}{nshuffle}{epoch}{2}],...
%         length([epoch_index{nmethod}{nshuffle}{epoch}{1} epoch_index{nmethod}{nshuffle}{epoch}{2}]));
%     resampled_event = [datasample(epoch_index{nmethod}{nshuffle}{epoch}{1},...
%         length(epoch_index{nmethod}{nshuffle}{epoch}{1})) datasample(epoch_index{nmethod}{nshuffle}{epoch}{2},...
%         length(epoch_index{nmethod}{nshuffle}{epoch}{2}))];
     resampled_event = [epoch_index{nmethod}{nshuffle}{epoch}{1} epoch_index{nmethod}{nshuffle}{epoch}{2}];

    this_segment_event = resampled_event(ismember(resampled_event,find(segment_id{nmethod}{nshuffle} == 1)));
    thresholded_event = this_segment_event(find(log_pval{nmethod}{nshuffle}(this_segment_event) <= p_val_threshold(threshold)));
    this_segment_event = resampled_event(ismember(resampled_event,find(segment_id{nmethod}{nshuffle} > 1)));
    thresholded_event = [thresholded_event this_segment_event(find(log_pval{nmethod}{nshuffle}(this_segment_event) <= p_val_threshold(threshold)/2))];

    %                         thresholded_event = p_value_thresholding(resampled_event,segment_id{nmethod}{nshuffle},log_pval{nmethod}{nshuffle},p_val_threshold(threshold));

    track_1_index = thresholded_event(ismember(thresholded_event,epoch_index{nmethod}{nshuffle}{epoch}{1}));
    track_2_index = thresholded_event(ismember(thresholded_event,epoch_index{nmethod}{nshuffle}{epoch}{2}));


%     log_odd_difference(threshold,nboot) = mean(data(track_1_index)) ...
%         - mean(data(track_2_index));

    %                         [multi_event_percent(nboot) multi_event_number(nboot)] = calculate_multitrack_event_percentage(session_index{nmethod}{nshuffle},event_index{nmethod}{nshuffle},track_1_index,track_2_index,epoch)
    [multi_event_percent multi_event_number] = calculate_multitrack_event_percentage(session_index{nmethod}{nshuffle},event_index{nmethod}{nshuffle},track_1_index,track_2_index);
    percent_multi_events(threshold,nboot) = multi_event_percent;
    % Significant event proportion (minus multitrack event number to avoid double counting)
    percent_sig_events(threshold,nboot) = (length(track_1_index)+length(track_2_index)-multi_event_number)/total_number(epoch);
end