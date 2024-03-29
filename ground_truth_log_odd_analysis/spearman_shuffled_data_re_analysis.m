function [] = spearman_shuffled_data_re_analysis(folders,option)

if strcmp(option,'global_remapped')
    data_folders = 'global_remapped';
elseif strcmp(option,'place_field_shifted')
    data_folders = 'place_field';
elseif strcmp(option,'spike_train_shifted')
    data_folders = 'spike_train';
elseif strcmp(option,'cross_experiment_shuffled')
    data_folders = 'cross_experiment';
end

for nfolder = 1:length(folders)

    cd(folders{nfolder})
    disp(folders{nfolder})
    % spearman correlation
    load extracted_place_fields_BAYESIAN

    cd([data_folders,'_shuffles'])
    DIR = dir('shuffle_*');
    DataPath = natsortfiles({DIR.name})'


    for nshuffle = 1:length(DataPath)
        
        shuffle_folder = sprintf('shuffle_%i',nshuffle);
        disp(shuffle_folder)
        cd(shuffle_folder)
        if strcmp(option,'global_remapped')
            load global_remapped_place_fields_id
            shuffled_data_place_fields_id = global_remapped_place_fields_id;
        else
            load cross_experiment_shuffled_place_fields_id
            shuffled_data_place_fields_id = cross_experiment_shuffled_place_fields_id;
        end
        load decoded_replay_events
        load scored_replay
        data = decoded_replay_events;
        scored_replay(1).replay_events(length(data(1).replay_events)).spearman_score = [];
        scored_replay(2).replay_events(length(data(1).replay_events)).spearman_score = [];

        for event = 1:length(data(1).replay_events)
            % Median Spike Spearman Correlation
            for track = 1 : length(place_fields_BAYESIAN.track)
                % Find sorted place cell id (but global remapped)
                place_cell_id = shuffled_data_place_fields_id{event}{track};
                sorted_place_fields = place_cell_id(2,:); % apply this index to the global remapped place cells (cell id shuffled).

                decoded_event = data(track).replay_events(event);
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
%         save scored_replay scored_replay
%         clear scored_replay
%         clear decoded_replay_events


        % Median Spikes Spearman Correlation (first half and second half of the 'whole' event)
        load  scored_replay_segments
        load decoded_replay_events_segments
        data1=decoded_replay_events1;
        data2=decoded_replay_events2;

        for track = 1 : length(place_fields_BAYESIAN.track)

            place_cell_id = shuffled_data_place_fields_id{event}{track};
            sorted_place_fields = place_cell_id(2,:); % apply this index to the global remapped place cells (cell id shuffled).

            for event = 1 : length(data(track).replay_events)
                decoded_event1 = data1(track).replay_events(event);
                decoded_event2 = data2(track).replay_events(event);

                if length(find(isnan(decoded_event1.decoded_position)))>0 | length(decoded_event1.timebins_centre)<5 | size(decoded_event1.spikes,1)==0
                    scored_replay1(track).replay_events(event).spearman_score=NaN;
                    scored_replay1(track).replay_events(event).spearman_p=1;
                else
                    spike_id=data1(track).replay_events(event).spikes(:,1);
                    spike_times=data1(track).replay_events(event).spikes(:,2);
                    [scored_replay1(track).replay_events(event).spearman_score,scored_replay1(track).replay_events(event).spearman_p] = spearman_median(spike_id, spike_times, sorted_place_fields);
                end

                if length(find(isnan(decoded_event2.decoded_position)))>0 | length(decoded_event2.timebins_centre)<5 | size(decoded_event2.spikes,1)==0
                    scored_replay2(track).replay_events(event).spearman_score=NaN;
                    scored_replay2(track).replay_events(event).spearman_p=1;
                else
                    spike_id=data2(track).replay_events(event).spikes(:,1);
                    spike_times=data2(track).replay_events(event).spikes(:,2);
                    [scored_replay2(track).replay_events(event).spearman_score,scored_replay2(track).replay_events(event).spearman_p] = spearman_median(spike_id, spike_times, sorted_place_fields);
                end
            end
        end

%         save scored_replay_segments scored_replay1 scored_replay2
%         clear scored_replay1 scored_replay2
%         clear decoded_replay_events1 decoded_replay_events2
        

        % All Spikes Spearman Correlation (whole)

%         load scored_replay
%         data = decoded_replay_events;
        scored_replay(1).replay_events(length(data(1).replay_events)).all_spearman_score = [];
        scored_replay(2).replay_events(length(data(1).replay_events)).all_spearman_score = [];
        for track = 1 : length(place_fields_BAYESIAN.track)

            place_cell_id = shuffled_data_place_fields_id{event}{track};
            sorted_place_fields = place_cell_id(2,:); % apply this index to the global remapped place cells (cell id shuffled).


            for event = 1 : length(data(track).replay_events)
                decoded_event = data(track).replay_events(event);
                if length(find(isnan(decoded_event.decoded_position)))>0 | length(decoded_event.timebins_centre)<5 | size(decoded_event.spikes,1)==0
                    scored_replay(track).replay_events(event).spearman_all_score=NaN;
                    scored_replay(track).replay_events(event).spearman_all_p=1;
                else
                    spike_id=data(track).replay_events(event).spikes(:,1);
                    spike_times=data(track).replay_events(event).spikes(:,2);
                    [scored_replay(track).replay_events(event).spearman_all_score,scored_replay(track).replay_events(event).spearman_all_p] = spearman_all_spikes(spike_id, spike_times, sorted_place_fields);
                end
            end
        end

        save scored_replay scored_replay
        clear scored_replay
%         clear decoded_replay_events

%         load scored_replay_segments
%         load decoded_replay_events_segments
%         data1=decoded_replay_events1;
%         data2=decoded_replay_events2;

        % All Spikes Spearman Correlation (Part 1 and Part 2)
        for track = 1 : length(place_fields_BAYESIAN.track)
            
            place_cell_id = shuffled_data_place_fields_id{event}{track};
            sorted_place_fields = place_cell_id(2,:); % apply this index to the global remapped place cells (cell id shuffled).

            for event = 1 : length(data(track).replay_events)
                decoded_event1 = data1(track).replay_events(event);
                decoded_event2 = data2(track).replay_events(event);

                if length(find(isnan(decoded_event1.decoded_position)))>0 | length(decoded_event1.timebins_centre)<5 | size(decoded_event1.spikes,1)==0
                    scored_replay1(track).replay_events(event).spearman_all_score=NaN;
                    scored_replay1(track).replay_events(event).spearman_all_p=1;
                else
                    spike_id=data1(track).replay_events(event).spikes(:,1);
                    spike_times=data1(track).replay_events(event).spikes(:,2);
                    [scored_replay1(track).replay_events(event).spearman_all_score,scored_replay1(track).replay_events(event).spearman_all_p] = spearman_all_spikes(spike_id, spike_times, sorted_place_fields);
                end

                if length(find(isnan(decoded_event2.decoded_position)))>0 | length(decoded_event2.timebins_centre)<5 | size(decoded_event2.spikes,1)==0
                    scored_replay2(track).replay_events(event).spearman_all_score=NaN;
                    scored_replay2(track).replay_events(event).spearman_all_p=1;
                else
                    spike_id=data2(track).replay_events(event).spikes(:,1);
                    spike_times=data2(track).replay_events(event).spikes(:,2);
                    [scored_replay2(track).replay_events(event).spearman_all_score,scored_replay2(track).replay_events(event).spearman_all_p] = spearman_all_spikes(spike_id, spike_times, sorted_place_fields);
                end
            end
        end

        save scored_replay_segments scored_replay1 scored_replay2
        clear scored_replay1 scored_replay2
        clear decoded_replay_events1 decoded_replay_events2
        cd ..

    end
    cd ..
    cd ..
end
cd ..
end


