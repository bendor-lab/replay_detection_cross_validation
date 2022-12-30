function out = calculate_replay_detection_performance(data,epoch_index,segment_id,log_pval,p_val_threshold,total_number,option)
%                 multi_event_index = find(log_odd_compare{nshuffle}.multi_track == 1);
out = [];
log_odd_difference = [];
percent_multi_events = [];
percent_sig_events = [];
if strcmp(option,'bootstrap')
    parfor threshold = 1:length(p_val_threshold)
        % find the event index with p value lower than the threshold
        %             current_index = find(log_pval{nmethod}(index{nmethod}{epoch})< p_val_threshold(threshold));
        difference = [];
        percent = [];
        track_2_index = [];
        track_1_index = [];

        for nboot = 1:1000 % Bootstrapping 1000 times
            resampled_event = datasample(epoch_index,length(epoch_index));

            %     for track = 1:2
            this_segment_event = find(segment_id(1,resampled_event) == 1);
            track_1_index = this_segment_event(find(log_pval(this_segment_event) <= p_val_threshold(threshold)));
            this_segment_event = find(segment_id(1,resampled_event) > 1);
            track_1_index = [track_1_index this_segment_event(find(log_pval(this_segment_event) <= p_val_threshold(threshold)/2))];

            this_segment_event = find(segment_id(2,resampled_event) == 1);
            track_2_index = this_segment_event(find(log_pval(this_segment_event) <= p_val_threshold(threshold)));
            this_segment_event = find(segment_id(2,resampled_event) > 1);
            track_2_index = [track_2_index this_segment_event(find(log_pval(this_segment_event) <= p_val_threshold(threshold)/2))];
            %     end
            multi_event_number = sum(ismember(track_1_index,track_2_index));
            percent_multi_events(threshold,nboot) = multi_event_number/(length(track_1_index)+length(track_2_index)-multi_event_number);
            
            log_odd_difference(threshold,nboot) = mean(data(1,track_1_index)) ...
                - mean(data(2,track_2_index));

            % Significant event proportion (minus multitrack event number to avoid double counting)
            percent_sig_events(threshold,nboot) = (length(track_1_index)+length(track_2_index)-multi_event_number)/total_number;
        end
    end

elseif strcmp(option,'original')
    
    for threshold = 1:length(p_val_threshold)
        all_events = epoch_index;

        %     for track = 1:2
        this_segment_event = find(segment_id(1,all_events) == 1);
        track_1_index = this_segment_event(find(log_pval(this_segment_event) <= p_val_threshold(threshold)));
        this_segment_event = find(segment_id(1,all_events) > 1);
        track_1_index = [track_1_index this_segment_event(find(log_pval(this_segment_event) <= p_val_threshold(threshold)/2))];

        this_segment_event = find(segment_id(2,all_events) == 1);
        track_2_index = this_segment_event(find(log_pval(this_segment_event) <= p_val_threshold(threshold)));
        this_segment_event = find(segment_id(2,all_events) > 1);
        track_2_index = [track_2_index this_segment_event(find(log_pval(this_segment_event) <= p_val_threshold(threshold)/2))];
        %     end
        multi_event_number = sum(ismember(track_1_index,track_2_index));
        percent_multi_events(threshold) = multi_event_number/(length(track_1_index)+length(track_2_index)-multi_event_number);

        log_odd_difference(threshold) = mean(data(1,track_1_index)) ...
            - mean(data(2,track_2_index));

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

% for nboot = 1:1 % Bootstrapping 1000 times
%     %     resampled_event = datasample([epoch_index{nmethod}{nshuffle}{epoch}{1} epoch_index{nmethod}{nshuffle}{epoch}{2}],...
%     %         length([epoch_index{nmethod}{nshuffle}{epoch}{1} epoch_index{nmethod}{nshuffle}{epoch}{2}]));
%     %     resampled_event = [datasample(epoch_index{nmethod}{nshuffle}{epoch}{1},...
%     %         length(epoch_index{nmethod}{nshuffle}{epoch}{1})) datasample(epoch_index{nmethod}{nshuffle}{epoch}{2},...
%     %         length(epoch_index{nmethod}{nshuffle}{epoch}{2}))];
%     resampled_event = epoch_index{nmethod}{nshuffle}{epoch};
% 
%     %     for track = 1:2
%     this_segment_event = find(segment_id{nmethod}{nshuffle}(1,resampled_event) == 1);
%     track_1_index = this_segment_event(find(log_pval{nmethod}{nshuffle}(this_segment_event) <= p_val_threshold(threshold)));
%     this_segment_event = find(segment_id{nmethod}{nshuffle}(1,resampled_event) > 1);
%     track_1_index = [track_1_index this_segment_event(find(log_pval{nmethod}{nshuffle}(this_segment_event) <= p_val_threshold(threshold)/2))];
% 
%     this_segment_event = find(segment_id{nmethod}{nshuffle}(2,resampled_event) == 1);
%     track_2_index = this_segment_event(find(log_pval{nmethod}{nshuffle}(this_segment_event) <= p_val_threshold(threshold)));
%     this_segment_event = find(segment_id{nmethod}{nshuffle}(2,resampled_event) > 1);
%     track_2_index = [track_2_index this_segment_event(find(log_pval{nmethod}{nshuffle}(this_segment_event) <= p_val_threshold(threshold)/2))];
%     %     end
%     multi_event_number = sum(ismember(track_1_index,track_2_index));
%     multi_event_percent = multi_event_number/(length(track_1_index)+length(track_2_index)-multi_event_number);
% 
%     log_odd_difference(threshold) = mean(data{nmethod}{nshuffle}(1,track_1_index)) ...
%         - mean(data{nmethod}{nshuffle}(2,track_2_index));
% 
%     percent_multi_events(threshold,nboot) = multi_event_percent;
%     % Significant event proportion (minus multitrack event number to avoid double counting)
%     percent_sig_events(threshold,nboot) = (length(track_1_index)+length(track_2_index)-multi_event_number)/total_number(epoch);
% end
