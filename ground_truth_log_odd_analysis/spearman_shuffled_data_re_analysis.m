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
    load replayEvents_bayesian_spike_count

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
        elseif strcmp(option,'cross_experiment_shuffled')
            load cross_experiment_shuffled_place_fields_id
            shuffled_data_place_fields_id = cross_experiment_shuffled_place_fields_id;
        elseif strcmp(option,'spike_train_shifted')
            load replayEvents_shifted_bayesian_spike_count
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
                if strcmp(option,'global_remapped') | strcmp(option,'cross_experiment_shuffled')
                    place_cell_id = shuffled_data_place_fields_id{event}{track};
                    sorted_place_fields = place_cell_id(2,:); % apply this index to the global remapped place cells (cell id shuffled).
                elseif strcmp(option,'spike_train_shifted')
                    sorted_place_fields = place_fields_BAYESIAN.track(track).sorted_good_cells;
                    unsorted_place_cell_id = place_fields_BAYESIAN.good_place_cells(find(ismember(place_fields_BAYESIAN.good_place_cells,sorted_place_fields)));

                    this_shuffled_event_spike_count = replayEvents_shifted_bayesian_spike_count.n.replay(ismember(place_fields_BAYESIAN.good_place_cells,sorted_place_fields),...
                        replayEvents_shifted_bayesian_spike_count.replay_events_indices==event);
                    this_event_spike_count = replayEvents_bayesian_spike_count.n.replay(ismember(place_fields_BAYESIAN.good_place_cells,sorted_place_fields),...
                        replayEvents_shifted_bayesian_spike_count.replay_events_indices==event);
                    spikes = data(track).replay_events(event).spikes;

                    for ncell = 1:length(unsorted_place_cell_id)
                        if sum(this_event_spike_count(ncell,:))~=0
                            %                             for nbin = 1:size(this_event_spike_count,2)
                            %                                 conv(circshift(this_event_spike_count(ncell,:),nbin),this_shuffled_event_spike_count(ncell,:));
                            %                             cir_xcorr(ncell,:) = cconv(this_event_spike_count(ncell,:),this_event_spike_count(ncell,:));
                            cir_xcorr(ncell,:) = xcorr(this_event_spike_count(ncell,:),this_shuffled_event_spike_count(ncell,:));
                            %                             [~,index(ncell)]=max(cir_xcorr(ncell,:));
                            [~,index(ncell)]=max(cir_xcorr(ncell,:));
                            index(ncell)= -(index(ncell)-size(this_event_spike_count,2));
                            this_cell_spike_times = spikes(spikes(:,1)==unsorted_place_cell_id(ncell),2);
                            for nspike = 1:length(this_cell_spike_times)

                                if this_cell_spike_times(nspike)+index(ncell)*0.02 > data(track).replay_events(event).timebins_edges(2)
                                    this_cell_spike_times(nspike) = this_cell_spike_times(nspike)+index(ncell)*0.02-data(track).replay_events(event).timebins_edges(2) + ...
                                        data(track).replay_events(event).timebins_edges(1);
                                else
                                    this_cell_spike_times(nspike) = this_cell_spike_times(nspike)+index(ncell)*0.02;
                                end
                            end

                            spikes(spikes(:,1)==unsorted_place_cell_id(ncell),2) = this_cell_spike_times;
                            %                         index(ncell)=

                            %                             end
                        else
                            index(ncell)=nan;
                        end
                    end
                    [~,index]= sort(spikes(:,2));
                    spikes(:,2)=spikes(index,2);
                    spikes(:,1)=spikes(index,1);

                    data(track).replay_events(event).spikes = spikes;

                elseif strcmp(option,'place_field_shifted')
                    sorted_good_cells = place_fields_BAYESIAN.track(track).sorted_good_cells;
                    for ncell = 1:length(sorted_good_cells)
                        [~, shifted_place_fields_peak(ncell)] = max(place_fields_BAYESIAN.track(1).place_field_shifted{event}{sorted_good_cells(ncell)});
                    end
                    [shifted_place_fields_peak,index] = sort(shifted_place_fields_peak);
                    sorted_place_fields = sorted_good_cells(index);
                end

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

            % Find sorted place cell id (but global remapped)
            if strcmp(option,'global_remapped') | strcmp(option,'cross_experiment_shuffled')
                place_cell_id = shuffled_data_place_fields_id{event}{track};
                sorted_place_fields = place_cell_id(2,:); % apply this index to the global remapped place cells (cell id shuffled).
            elseif strcmp(option,'spike_train_shifted')
                sorted_place_fields = place_fields_BAYESIAN.track(track).sorted_good_cells;
                unsorted_place_cell_id = place_fields_BAYESIAN.good_place_cells(find(ismember(place_fields_BAYESIAN.good_place_cells,sorted_place_fields)));

                this_shuffled_event_spike_count = replayEvents_shifted_bayesian_spike_count.n.replay(ismember(place_fields_BAYESIAN.good_place_cells,sorted_place_fields),...
                    replayEvents_shifted_bayesian_spike_count.replay_events_indices==event);
                this_event_spike_count = replayEvents_bayesian_spike_count.n.replay(ismember(place_fields_BAYESIAN.good_place_cells,sorted_place_fields),...
                    replayEvents_shifted_bayesian_spike_count.replay_events_indices==event);
                spikes = data(track).replay_events(event).spikes;

                for ncell = 1:length(unsorted_place_cell_id)
                    if sum(this_event_spike_count(ncell,:))~=0
                        %                             for nbin = 1:size(this_event_spike_count,2)
                        %                                 conv(circshift(this_event_spike_count(ncell,:),nbin),this_shuffled_event_spike_count(ncell,:));
                        %                             cir_xcorr(ncell,:) = cconv(this_event_spike_count(ncell,:),this_event_spike_count(ncell,:));
                        cir_xcorr(ncell,:) = xcorr(this_event_spike_count(ncell,:),this_shuffled_event_spike_count(ncell,:));
                        %                             [~,index(ncell)]=max(cir_xcorr(ncell,:));
                        [~,index(ncell)]=max(cir_xcorr(ncell,:));
                        index(ncell)= -(index(ncell)-size(this_event_spike_count,2));
                        this_cell_spike_times = spikes(spikes(:,1)==unsorted_place_cell_id(ncell),2);
                        for nspike = 1:length(this_cell_spike_times)

                            if this_cell_spike_times(nspike)+index(ncell)*0.02 > data(track).replay_events(event).timebins_edges(2)
                                this_cell_spike_times(nspike) = this_cell_spike_times(nspike)+index(ncell)*0.02-data(track).replay_events(event).timebins_edges(2) + ...
                                    data(track).replay_events(event).timebins_edges(1);
                            else
                                this_cell_spike_times(nspike) = this_cell_spike_times(nspike)+index(ncell)*0.02;
                            end
                        end

                        spikes(spikes(:,1)==unsorted_place_cell_id(ncell),2) = this_cell_spike_times;
                        %                         index(ncell)=

                        %                             end
                    else
                        index(ncell)=nan;
                    end
                end
                [~,index]= sort(spikes(:,2));
                spikes(:,2)=spikes(index,2);
                spikes(:,1)=spikes(index,1);

                data(track).replay_events(event).spikes = spikes;
                %                 data1(track).replay_events(event).

                data1(track).replay_events(event).spikes = spikes(spikes(:,2)<data1(track).replay_events(event).midpoint,:);
                data2(track).replay_events(event).spikes = spikes(spikes(:,2)>data1(track).replay_events(event).midpoint,:);

            elseif strcmp(option,'place_field_shifted')
                sorted_good_cells = place_fields_BAYESIAN.track(track).sorted_good_cells;
                for ncell = 1:length(sorted_good_cells)
                    [~, shifted_place_fields_peak(ncell)] = max(place_fields_BAYESIAN.track(1).place_field_shifted{event}{sorted_good_cells(ncell)});
                end
                [shifted_place_fields_peak,index] = sort(shifted_place_fields_peak);
                sorted_place_fields = sorted_good_cells(index);
            end

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

            if strcmp(option,'global_remapped')
                load global_remapped_place_fields_id
                shuffled_data_place_fields_id = global_remapped_place_fields_id;
            elseif strcmp(option,'cross_experiment_shuffled')
                load cross_experiment_shuffled_place_fields_id
                shuffled_data_place_fields_id = cross_experiment_shuffled_place_fields_id;
            elseif strcmp(option,'spike_train_shifted')
                load replayEvents_shifted_bayesian_spike_count
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
                    if strcmp(option,'global_remapped') | strcmp(option,'cross_experiment_shuffled')
                        place_cell_id = shuffled_data_place_fields_id{event}{track};
                        sorted_place_fields = place_cell_id(2,:); % apply this index to the global remapped place cells (cell id shuffled).
                    elseif strcmp(option,'spike_train_shifted')
                        sorted_place_fields = place_fields_BAYESIAN.track(track).sorted_good_cells;
                        unsorted_place_cell_id = place_fields_BAYESIAN.good_place_cells(find(ismember(place_fields_BAYESIAN.good_place_cells,sorted_place_fields)));

                        this_shuffled_event_spike_count = replayEvents_shifted_bayesian_spike_count.n.replay(ismember(place_fields_BAYESIAN.good_place_cells,sorted_place_fields),...
                            replayEvents_shifted_bayesian_spike_count.replay_events_indices==event);
                        this_event_spike_count = replayEvents_bayesian_spike_count.n.replay(ismember(place_fields_BAYESIAN.good_place_cells,sorted_place_fields),...
                            replayEvents_shifted_bayesian_spike_count.replay_events_indices==event);
                        spikes = data(track).replay_events(event).spikes;

                        for ncell = 1:length(unsorted_place_cell_id)
                            if sum(this_event_spike_count(ncell,:))~=0
                                %                             for nbin = 1:size(this_event_spike_count,2)
                                %                                 conv(circshift(this_event_spike_count(ncell,:),nbin),this_shuffled_event_spike_count(ncell,:));
                                %                             cir_xcorr(ncell,:) = cconv(this_event_spike_count(ncell,:),this_event_spike_count(ncell,:));
                                cir_xcorr(ncell,:) = xcorr(this_event_spike_count(ncell,:),this_shuffled_event_spike_count(ncell,:));
                                %                             [~,index(ncell)]=max(cir_xcorr(ncell,:));
                                [~,index(ncell)]=max(cir_xcorr(ncell,:));
                                index(ncell)= -(index(ncell)-size(this_event_spike_count,2));
                                this_cell_spike_times = spikes(spikes(:,1)==unsorted_place_cell_id(ncell),2);
                                for nspike = 1:length(this_cell_spike_times)

                                    if this_cell_spike_times(nspike)+index(ncell)*0.02 > data(track).replay_events(event).timebins_edges(2)
                                        this_cell_spike_times(nspike) = this_cell_spike_times(nspike)+index(ncell)*0.02-data(track).replay_events(event).timebins_edges(2) + ...
                                            data(track).replay_events(event).timebins_edges(1);
                                    else
                                        this_cell_spike_times(nspike) = this_cell_spike_times(nspike)+index(ncell)*0.02;
                                    end
                                end

                                spikes(spikes(:,1)==unsorted_place_cell_id(ncell),2) = this_cell_spike_times;
                                %                         index(ncell)=

                                %                             end
                            else
                                index(ncell)=nan;
                            end
                        end
                        [~,index]= sort(spikes(:,2));
                        spikes(:,2)=spikes(index,2);
                        spikes(:,1)=spikes(index,1);

                        data(track).replay_events(event).spikes = spikes;

                    elseif strcmp(option,'place_field_shifted')
                        sorted_good_cells = place_fields_BAYESIAN.track(track).sorted_good_cells;
                        for ncell = 1:length(sorted_good_cells)
                            [~, shifted_place_fields_peak(ncell)] = max(place_fields_BAYESIAN.track(1).place_field_shifted{event}{sorted_good_cells(ncell)});
                        end
                        [shifted_place_fields_peak,index] = sort(shifted_place_fields_peak);
                        sorted_place_fields = sorted_good_cells(index);
                    end

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

                    if strcmp(option,'global_remapped')
                        load global_remapped_place_fields_id
                        shuffled_data_place_fields_id = global_remapped_place_fields_id;
                    elseif strcmp(option,'cross_experiment_shuffled')
                        load cross_experiment_shuffled_place_fields_id
                        shuffled_data_place_fields_id = cross_experiment_shuffled_place_fields_id;
                    elseif strcmp(option,'spike_train_shifted')
                        load replayEvents_shifted_bayesian_spike_count
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
                            if strcmp(option,'global_remapped') | strcmp(option,'cross_experiment_shuffled')
                                place_cell_id = shuffled_data_place_fields_id{event}{track};
                                sorted_place_fields = place_cell_id(2,:); % apply this index to the global remapped place cells (cell id shuffled).
                            elseif strcmp(option,'spike_train_shifted')
                                sorted_place_fields = place_fields_BAYESIAN.track(track).sorted_good_cells;
                                unsorted_place_cell_id = place_fields_BAYESIAN.good_place_cells(find(ismember(place_fields_BAYESIAN.good_place_cells,sorted_place_fields)));

                                this_shuffled_event_spike_count = replayEvents_shifted_bayesian_spike_count.n.replay(ismember(place_fields_BAYESIAN.good_place_cells,sorted_place_fields),...
                                    replayEvents_shifted_bayesian_spike_count.replay_events_indices==event);
                                this_event_spike_count = replayEvents_bayesian_spike_count.n.replay(ismember(place_fields_BAYESIAN.good_place_cells,sorted_place_fields),...
                                    replayEvents_shifted_bayesian_spike_count.replay_events_indices==event);
                                spikes = data(track).replay_events(event).spikes;

                                for ncell = 1:length(unsorted_place_cell_id)
                                    if sum(this_event_spike_count(ncell,:))~=0
                                        %                             for nbin = 1:size(this_event_spike_count,2)
                                        %                                 conv(circshift(this_event_spike_count(ncell,:),nbin),this_shuffled_event_spike_count(ncell,:));
                                        %                             cir_xcorr(ncell,:) = cconv(this_event_spike_count(ncell,:),this_event_spike_count(ncell,:));
                                        cir_xcorr(ncell,:) = xcorr(this_event_spike_count(ncell,:),this_shuffled_event_spike_count(ncell,:));
                                        %                             [~,index(ncell)]=max(cir_xcorr(ncell,:));
                                        [~,index(ncell)]=max(cir_xcorr(ncell,:));
                                        index(ncell)= -(index(ncell)-size(this_event_spike_count,2));
                                        this_cell_spike_times = spikes(spikes(:,1)==unsorted_place_cell_id(ncell),2);
                                        for nspike = 1:length(this_cell_spike_times)

                                            if this_cell_spike_times(nspike)+index(ncell)*0.02 > data(track).replay_events(event).timebins_edges(2)
                                                this_cell_spike_times(nspike) = this_cell_spike_times(nspike)+index(ncell)*0.02-data(track).replay_events(event).timebins_edges(2) + ...
                                                    data(track).replay_events(event).timebins_edges(1);
                                            else
                                                this_cell_spike_times(nspike) = this_cell_spike_times(nspike)+index(ncell)*0.02;
                                            end
                                        end

                                        spikes(spikes(:,1)==unsorted_place_cell_id(ncell),2) = this_cell_spike_times;
                                        %                         index(ncell)=

                                        %                             end
                                    else
                                        index(ncell)=nan;
                                    end
                                end
                                [~,index]= sort(spikes(:,2));
                                spikes(:,2)=spikes(index,2);
                                spikes(:,1)=spikes(index,1);

                                data(track).replay_events(event).spikes = spikes;

                            elseif strcmp(option,'place_field_shifted')
                                sorted_good_cells = place_fields_BAYESIAN.track(track).sorted_good_cells;
                                for ncell = 1:length(sorted_good_cells)
                                    [~, shifted_place_fields_peak(ncell)] = max(place_fields_BAYESIAN.track(1).place_field_shifted{event}{sorted_good_cells(ncell)});
                                end
                                [shifted_place_fields_peak,index] = sort(shifted_place_fields_peak);
                                sorted_place_fields = sorted_good_cells(index);
                            end
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


