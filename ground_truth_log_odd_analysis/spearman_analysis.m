function [] = spearman_analysis(folders)


for nfolder = 1:length(folders)
    
    cd(folders{nfolder})
    
    % spearman correlation
    load extracted_place_fields_BAYESIAN
    load scored_replay
    load('decoded_replay_events.mat');
    data=decoded_replay_events;

%     % Median Spike Spearman Correlation
%     for track = 1 : length(place_fields_BAYESIAN.track)
%         sorted_place_fields=place_fields_BAYESIAN.track(track).sorted_good_cells;
% 
%         for event = 1 : length(data(track).replay_events)
%             decoded_event = data(track).replay_events(event);
%             if length(find(isnan(decoded_event.decoded_position)))>0 | length(decoded_event.timebins_centre)<5 | size(decoded_event.spikes,1)==0
%                 scored_replay(track).replay_events(event).spearman_score=NaN;
%                 scored_replay(track).replay_events(event).spearman_p=1;
%             else
%                 spike_id=data(track).replay_events(event).spikes(:,1);
%                 spike_times=data(track).replay_events(event).spikes(:,2);
%                 [scored_replay(track).replay_events(event).spearman_score,scored_replay(track).replay_events(event).spearman_p] = spearman_median(spike_id, spike_times, sorted_place_fields);
%             end
%         end
%     end
    
    % Median and All Spikes Spearman Correlation (whole)
    for track = 1 : length(place_fields_BAYESIAN.track)
        sorted_place_fields=place_fields_BAYESIAN.track(track).sorted_good_cells;
        
        for event = 1 : length(data(track).replay_events)
            decoded_event = data(track).replay_events(event);
            if length(find(isnan(decoded_event.decoded_position)))>0 | length(decoded_event.timebins_centre)<5 | size(decoded_event.spikes,1)==0
                scored_replay(track).replay_events(event).spearman_all_score=NaN;
                scored_replay(track).replay_events(event).spearman_all_p=1;

                scored_replay(track).replay_events(event).spearman_score=NaN;
                scored_replay(track).replay_events(event).spearman_p=1;
            else
                spike_id=data(track).replay_events(event).spikes(:,1);
                spike_times=data(track).replay_events(event).spikes(:,2);
                [scored_replay(track).replay_events(event).spearman_all_score,scored_replay(track).replay_events(event).spearman_all_p] = spearman_all_spikes(spike_id, spike_times, sorted_place_fields);
                [scored_replay(track).replay_events(event).spearman_score,scored_replay(track).replay_events(event).spearman_p] = spearman_median(spike_id, spike_times, sorted_place_fields);

            end


        end
    end
    
    save scored_replay scored_replay
    clear scored_replay
    clear decoded_replay_events
    
    load  scored_replay_segments
    load decoded_replay_events_segments
    data1=decoded_replay_events1;
    data2=decoded_replay_events2;
    
    % Median and All Spikes Spearman Correlation (first half and second half of the 'whole' event)
    for track = 1 : length(place_fields_BAYESIAN.track)
        sorted_place_fields=place_fields_BAYESIAN.track(track).sorted_good_cells;
        
        for event = 1 : length(data(track).replay_events)
            decoded_event1 = data1(track).replay_events(event);
            decoded_event2 = data2(track).replay_events(event);
                        
            if length(find(isnan(decoded_event1.decoded_position)))>0 | length(decoded_event1.timebins_centre)<5 | size(decoded_event1.spikes,1)==0
                scored_replay1(track).replay_events(event).spearman_all_score=NaN;
                scored_replay1(track).replay_events(event).spearman_all_p=1;

                scored_replay1(track).replay_events(event).spearman_score=NaN;
                scored_replay1(track).replay_events(event).spearman_p=1;                
            else
                spike_id=data1(track).replay_events(event).spikes(:,1);
                spike_times=data1(track).replay_events(event).spikes(:,2);
                [scored_replay1(track).replay_events(event).spearman_all_score,scored_replay1(track).replay_events(event).spearman_all_p] = spearman_all_spikes(spike_id, spike_times, sorted_place_fields);
                [scored_replay1(track).replay_events(event).spearman_score,scored_replay1(track).replay_events(event).spearman_p] = spearman_median(spike_id, spike_times, sorted_place_fields);
            
            end
            
            if length(find(isnan(decoded_event2.decoded_position)))>0 | length(decoded_event2.timebins_centre)<5 | size(decoded_event2.spikes,1)==0
                scored_replay2(track).replay_events(event).spearman_all_score=NaN;
                scored_replay2(track).replay_events(event).spearman_all_p=1;

                scored_replay2(track).replay_events(event).spearman_score=NaN;
                scored_replay2(track).replay_events(event).spearman_p=1;                  
            else
                spike_id=data2(track).replay_events(event).spikes(:,1);
                spike_times=data2(track).replay_events(event).spikes(:,2);
                [scored_replay2(track).replay_events(event).spearman_all_score,scored_replay2(track).replay_events(event).spearman_all_p] = spearman_all_spikes(spike_id, spike_times, sorted_place_fields);
                [scored_replay2(track).replay_events(event).spearman_score,scored_replay2(track).replay_events(event).spearman_p] = spearman_median(spike_id, spike_times, sorted_place_fields);
                        
            end
        end
    end
    save decoded_replay_events_segments decoded_replay_events1 decoded_replay_events2
    save scored_replay_segments scored_replay1 scored_replay2
    clear scored_replay1 scored_replay2
    clear decoded_replay_events1 decoded_replay_events2
    cd ..

end
    
    




end

