function estimated_position = log_odd_bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count,place_field_index,timebin,posbin,option,modify,save_option)
% INPUTS:
% place fields: matrix or []. Place fields to use as a template in the decoding. If empty, it will load place fields for the whole session (extracted_place_fields_BAYESIAN.mat)
% bayesian_spike_count: Maybe original replay spike count, common good replay spike count or rate-fixed common good cell replay spike count
% place_field_index: Original (all good cells) or Stable common good cells (based on subsets_of_cells)
% timebin: If 0.02, then it means 20ms. If 1, then it means 1 timebin for 1 entire event
% posbin: If 20, then it means one track is divided into 20 position bins. If 1, then 1 position bin covers the entirety of the track.
% Save_option: enter 'Y' for saving. Else, won't save.
% OUTPUTS:
% estimated_position structure.


parameters=list_of_parameters;
% load extracted_position

if strcmp(bayesian_spike_count, 'replayEvents_bayesian_spike_count')
    load replayEvents_bayesian_spike_count  %if decoding replay events (which will have different time bins)
    bayesian_spike_count = replayEvents_bayesian_spike_count;
    
elseif strcmp(bayesian_spike_count, 'replayEvents_common_good_cell_bayesian_spike_count')
    load replayEvents_common_good_cell_bayesian_spike_count
    bayesian_spike_count = replayEvents_common_good_cell_bayesian_spike_count;
elseif strcmp(bayesian_spike_count, 'replayEvents_shifted_bayesian_spike_count')
    load replayEvents_shifted_bayesian_spike_count
    bayesian_spike_count = replayEvents_shifted_bayesian_spike_count;
elseif strcmp(bayesian_spike_count, 'replayEvents_common_cell_shifted_bayesian_spike_count')
    load replayEvents_common_cell_shifted_bayesian_spike_count
    bayesian_spike_count = replayEvents_common_cell_shifted_bayesian_spike_count;

elseif strcmp(bayesian_spike_count, 'circular_shift_spike_count') % For noise in replay model
    load circular_shift_spike_count
    bayesian_spike_count = circular_shift_spike_count;
elseif strcmp(bayesian_spike_count, 'jitter_noise_spike_count') % For noise in replay model
    load jitter_noise_spike_count
    bayesian_spike_count = jitter_noise_spike_count;
elseif strcmp(bayesian_spike_count,'replayEvents_FIXED_spike_count_20ms')
    cd 'rate_fixed'
    load replayEvents_FIXED_spike_count_20ms
    bayesian_spike_count = replayEvents_FIXED_spike_count;
    cd ..
elseif strcmp(bayesian_spike_count,'replayEvents_FIXED_spike_count_one_bin')
    cd 'rate_fixed'
    load replayEvents_FIXED_spike_count_one_bin
    bayesian_spike_count = replayEvents_FIXED_spike_count;
    cd ..
elseif strcmp(bayesian_spike_count,'shuffled_spike_count')
    load shuffled_spike_count
    bayesian_spike_count = shuffled_spike_count;
elseif isempty(bayesian_spike_count) %if decoding the whole session
    load bayesian_spike_count
    place_field_index = [];
end

if isempty(timebin)
    bin_width = parameters.replay_bin_width;
elseif isnumeric(timebin) && timebin >= 1
    % If time bin is n time bins per event
    bin_width = [];
    for i = 1:length(bayesian_spike_count.replay_events)
        % Get the bin widths for the time bin in each event bin
        % (if timebin = 1, then only one bin width per event)
        for m = 1:length(bayesian_spike_count.replay_events(i).replay_time_edges)-1
            bin_width{i}(m) = bayesian_spike_count.replay_events(i).replay_time_edges(m) - bayesian_spike_count.replay_events(i).replay_time_edges(m+1);
        end
    end
    m = [];
    
else
    % if timebin is less than 1, then it is specifing ms time lenghth
    bin_width = timebin;
end

%%%%% BAYESIAN DECODING  %%%%%%

track_no = 2; % Track 1 and Track 2 (and/or Shuffled 1 and Shuffled 2)

if isnumeric(posbin) && ~isempty(posbin)
    estimated_position(1).position_bin_centres = [];
    position_bin_centres = [];
    posbin_edges = linspace(0,200*track_no,posbin*track_no+1);
    
    for bin = 1:length(posbin_edges)-1
        estimated_position(1).position_bin_centres(bin) = (posbin_edges(bin)+posbin_edges(bin+1))/2;
    end
    
%     estimated_position(1).discrete_position = NaN(size(position.linear(1).linear));
%     position_bin_edges = linspace(0,200,posbin+1);
%     
%     for track_id = 1:track_no
%         % Bin position for decoding error
%         discrete_position = discretize(position.linear(track_id).linear,position_bin_edges); %group position points in bins delimited by edges
%         for bin = 1:length(position_bin_edges)-1
%             position_bin_centres(bin) = (position_bin_edges(bin)+position_bin_edges(bin+1))/2;
%         end
%         
%         if track_id == 2
%             position_bin_centres = 200 + position_bin_centres;
%         end
%         
%         index = find(~isnan(discrete_position));
%         estimated_position(1).discrete_position(index) = position_bin_centres(discrete_position(index)); %creates new positions based on centre of bins
%     end
    
    %     estimated_position(2) = estimated_position(1);
    
elseif isempty(posbin)
    
    posbin = 20;
    estimated_position(1).position_bin_centres = [];
    position_bin_centres = [];
    posbin_edges = linspace(0,200*track_no,posbin*track_no+1);
    
    for bin = 1:length(posbin_edges)-1
        estimated_position(1).position_bin_centres(bin) = (posbin_edges(bin)+posbin_edges(bin+1))/2;
    end
    
%     estimated_position(1).discrete_position = NaN(size(position.linear(1).linear));
%     position_bin_edges = linspace(0,200,posbin+1);
%     
%     for track_id = 1:track_no
%         % Bin position for decoding error
%         discrete_position = discretize(position.linear(track_id).linear,position_bin_edges); %group position points in bins delimited by edges
%         for bin = 1:length(position_bin_edges)-1
%             position_bin_centres(bin) = (position_bin_edges(bin)+position_bin_edges(bin+1))/2;
%         end
%         
%         if track_id == 2
%             position_bin_centres = 200 + position_bin_centres;
%         end
%         
%         index = find(~isnan(discrete_position));
%         estimated_position(1).discrete_position(index) = position_bin_centres(discrete_position(index)); %creates new positions based on centre of bins
%     end
    
    posbin = [];
end


%     if strcmp(option,'track vs shuffle') | strcmp(option,'track and shuffle')
%         estimated_position(track_id+2).position_bin_centres = estimated_position(track_id).position_bin_centres
%         estimated_position(track_id+2).discrete_position = estimated_position(track_id).discrete_position
%     end


if isfield(bayesian_spike_count,'replay_events')   % When running replay events separately
    estimated_position(1).replay_events = bayesian_spike_count.replay_events;
    
    for i = 1 : length(bayesian_spike_count.replay_events)
        estimated_position(1).replay_events(i).replay = zeros(length(estimated_position(1).position_bin_centres),length(estimated_position(1).replay_events(i).replay_time_centered));
    end
    
    n.replay = bayesian_spike_count.n.replay;
    estimated_position(1).replay = zeros(length(estimated_position(1).position_bin_centres),length(n.replay));
    
    %         if strcmp(option,'track vs shuffle') | strcmp(option,'track and shuffle')
    %             estimated_position(track_id+2).replay_events = estimated_position(track_id).replay_events;
    %             estimated_position(track_id+2).replay = estimated_position(track_id).replay;
    %         end
    %     estimated_position(2) = estimated_position(1);
end


% Find ratemaps
if isempty(place_field_index)
    if isfield(place_fields_BAYESIAN,'good_place_cells')
        place_field_index = place_fields_BAYESIAN.good_place_cells;  %use all place cells
    elseif ~isfield(place_fields_BAYESIAN,'good_place_cells') && isfield(place_fields_BAYESIAN.track,'good_cells')
        place_field_index = place_fields_BAYESIAN.track.good_cells;  % when using output from get_place_fields_laps
    else
        disp('ERROR- field good_place_cells missing');
        place_field_index = place_fields_BAYESIAN.pyramidal_cells;
    end
end


all_place_fields = [];

if isnumeric(posbin) && ~isempty(posbin) % If the size of the position bin specified
    load Track_fields
    if isempty(modify) % Original data
%         Track 1 and Track 2 place fields concatenated into one
        all_place_fields{1} = [Track_fields(1).mean_rate(:,place_field_index)' Track_fields(2).mean_rate(:,place_field_index)'];
        
        for k = 1:length(place_field_index) % Turning curve shuffled for each cell (across two tracks)
            track_id = randperm(2);
            all_place_fields{2}(k,:) = [Track_fields(track_id(1)).mean_rate(:,place_field_index(k))' Track_fields(track_id(2)).mean_rate(:,place_field_index(k))'];
        end
        
    elseif strcmp(modify,'rate fixed') % Rate fixed place fields
%         Track 1 and Track 2 place fields concatenated into one
        all_place_fields{1} = [Track_fields(1).rate_fixed_mean_rate(:,place_field_index)' Track_fields(2).rate_fixed_mean_rate(:,place_field_index)'];
        
        for k = 1:length(place_field_index) % Turning curve shuffled for each cell (across two tracks)
            track_id = randperm(2);
            all_place_fields{2}(k,:) = [Track_fields(track_id(1)).rate_fixed_mean_rate(:,place_field_index(k))' Track_fields(track_id(2)).rate_fixed_mean_rate(:,place_field_index(k))'];
        end
        
        
        %             elseif strcmp(modify,'global remapping') % global remapped
        %                 % Track 1 and Track 2 place fields concatenated into one
        %                 all_place_fields{1} = [Track_fields(1).global_remapped_mean_rate(:,place_field_index)'...
        %                     Track_fields(2).global_remapped_mean_rate(:,place_field_index)'];
        %
        %                 for k = 1:length(place_field_index) % Turning curve shuffled for each cell (across two tracks)
        %                     track_id = randperm(2);
        %                     all_place_fields{2}(k,:) = [Track_fields(track_id(1)).global_remapped_mean_rate(:,place_field_index(k))'...
        %                         Track_fields(track_id(2)).global_remapped_mean_rate(:,place_field_index(k))'];
        %                 end
        
        
    elseif strcmp(modify,'rate remapping') % rate remapped
        all_place_fields = [];
        for event = 1:length(bayesian_spike_count.replay_events)
            % Track 1 and Track 2 place fields concatenated into one
            all_place_fields{1}{event} = [Track_fields(1).rate_remapped_mean_rate{event}(:,place_field_index)'...
                Track_fields(2).rate_remapped_mean_rate{event}(:,place_field_index)'];
            
            for k = 1:length(place_field_index) % Turning curve shuffled for each cell (across two tracks)
                track_id = randperm(2);
                all_place_fields{2}{event}(k,:) = [Track_fields(track_id(1)).rate_remapped_mean_rate{event}(:,place_field_index(k))'...
                    Track_fields(track_id(2)).rate_remapped_mean_rate{event}(:,place_field_index(k))'];
            end
        end
    end
    
else % if position bin not specified, then use original 20 position bins
    
    if isempty(modify) % Original data
        
        for track_id = 1:2
            for k=1:length(place_field_index)
                single_place_field = place_fields_BAYESIAN.track(track_id).raw{place_field_index(k)}; %get raw place field
                single_place_field(find(isnan(single_place_field))) = 0; % remove NaNs in place field and replace by 0
                % single_place_field = smooth(single_place_field,parameters.smoothing_number_of_bins); %smooth place field
                
                if min(single_place_field)<0
                    disp('error- spike rate of place field less than zero')
                end
                tempt{track_id}(k,:) = single_place_field;
            end
        end
        
        all_place_fields{1} = [tempt{1} tempt{2}];
        
        
        
    elseif strcmp(modify,'rate fixed') % Rate fixed place fields
        
        for track_id = 1:2
            for k=1:length(place_field_index)
                single_place_field = place_fields_BAYESIAN.track(track_id).rate_fixed_raw{place_field_index(k)}; %get raw place field
                single_place_field(find(isnan(single_place_field))) = 0; % remove NaNs in place field and replace by 0
                % single_place_field = smooth(single_place_field,parameters.smoothing_number_of_bins); %smooth place field
                
                if min(single_place_field)<0
                    disp('error- spike rate of place field less than zero')
                end
                tempt{track_id}(k,:) = single_place_field;
            end
        end
        
        all_place_fields{1} = [tempt{1} tempt{2}];
        %
        %     elseif strcmp(modify,'PRE merged_place_field_circular_shift')
        %
        %         num_cells = size(place_fields_BAYESIAN.track(1).rate_fixed_raw,2);
        %         %circular shift place field bins (two track place fields merged)
        %         for k = 1 : num_cells
        %             field =  [cell2mat(place_fields_BAYESIAN.track(1).rate_fixed_raw(k)) cell2mat(place_fields_BAYESIAN.track(2).rate_fixed_raw(k))];
        %             shuffled_place_fields(k,:) = {circshift(field, ceil(rand*length(field)))};
        %         end
        %
        %         all_place_fields{1} = shuffled_place_fields(place_field_index,:);
    
 
elseif strcmp(modify,'global remapping')
    
    for track_id = 1:2
        for k=1:length(place_field_index)
            single_place_field = place_fields_BAYESIAN.track(track_id).global_remapped_raw{place_field_index(k)}; %get raw place field
            single_place_field(find(isnan(single_place_field))) = 0; % remove NaNs in place field and replace by 0
            % single_place_field = smooth(single_place_field,parameters.smoothing_number_of_bins); %smooth place field
            
            if min(single_place_field)<0
                disp('error- spike rate of place field less than zero')
            end
            tempt{track_id}(k,:) = single_place_field;
        end
    end
    
    all_place_fields{1} = [tempt{1} tempt{2}];
    
    elseif strcmp(modify,'place_field_shifted')
        % for circular shift, get place fields for each event
        for event = 1:length(bayesian_spike_count.replay_events)
            for track_id = 1:2
                for k=1:length(place_field_index)
                    single_place_field = place_fields_BAYESIAN.track(track_id).place_field_shifted{event}{place_field_index(k)}; %get raw place field
                    single_place_field(find(isnan(single_place_field))) = 0; % remove NaNs in place field and replace by 0
                    % single_place_field = smooth(single_place_field,parameters.smoothing_number_of_bins); %smooth place field

                    if min(single_place_field)<0
                        disp('error- spike rate of place field less than zero')
                    end
                    tempt{track_id}(k,:) = single_place_field;
                end
            end

            all_place_fields{1}{event} = [tempt{1} tempt{2}];
        end


    elseif strcmp(modify,'cross_experiment_shuffled')
        % for circular shift, get place fields for each event
        for event = 1:length(bayesian_spike_count.replay_events)
            for track_id = 1:2
                for k=1:length(place_field_index)
                    single_place_field = place_fields_BAYESIAN.track(track_id).cross_experiment_shuffled{event}{place_field_index(k)}; %get raw place field
                    single_place_field(find(isnan(single_place_field))) = 0; % remove NaNs in place field and replace by 0
                    % single_place_field = smooth(single_place_field,parameters.smoothing_number_of_bins); %smooth place field

                    if min(single_place_field)<0
                        disp('error- spike rate of place field less than zero')
                    end
                    tempt{track_id}(k,:) = single_place_field;
                end
            end

            all_place_fields{1}{event} = [tempt{1} tempt{2}];
        end

    elseif strcmp(modify,'place field shift')
        % for circular shift, get place fields for each event
        for event = 1:length(bayesian_spike_count.replay_events)
            for track_id = 1:2
                for k=1:length(place_field_index)
                    single_place_field = place_fields_BAYESIAN.track(track_id).raw_shifted{event}{place_field_index(k)}; %get raw place field
                    single_place_field(find(isnan(single_place_field))) = 0; % remove NaNs in place field and replace by 0
                    % single_place_field = smooth(single_place_field,parameters.smoothing_number_of_bins); %smooth place field
                    
                    if min(single_place_field)<0
                        disp('error- spike rate of place field less than zero')
                    end
                    tempt{track_id}(k,:) = single_place_field;
                end
            end
            
            all_place_fields{1}{event} = [tempt{1} tempt{2}];
        end

    elseif strcmp(modify,'track label swap')
        % for circular shift, get place fields for each event
        for event = 1:length(bayesian_spike_count.replay_events)
            for track_id = 1:2
                for k=1:length(place_field_index)
                    single_place_field = place_fields_BAYESIAN.track(track_id).raw_swapped{event}{place_field_index(k)}; %get raw place field
                    single_place_field(find(isnan(single_place_field))) = 0; % remove NaNs in place field and replace by 0
                    % single_place_field = smooth(single_place_field,parameters.smoothing_number_of_bins); %smooth place field

                    if min(single_place_field)<0
                        disp('error- spike rate of place field less than zero')
                    end
                    tempt{track_id}(k,:) = single_place_field;
                end
            end

            all_place_fields{1}{event} = [tempt{1} tempt{2}];
        end

    elseif strcmp(modify,'global remapping per event')
        %         cd 'global_remapped'
        %         load extracted_place_fields_BAYESIAN
        %         cd ..
        
        % for global remapping, get place fields for each event
        for event = 1:length(bayesian_spike_count.replay_events)
            for track_id = 1:2
                for k=1:length(place_field_index)
                    single_place_field = place_fields_BAYESIAN.track(track_id).global_remapped{event}{place_field_index(k)}; %get raw place field
                    single_place_field(find(isnan(single_place_field))) = 0; % remove NaNs in place field and replace by 0
                    % single_place_field = smooth(single_place_field,parameters.smoothing_number_of_bins); %smooth place field
                    
                    if min(single_place_field)<0
                        disp('error- spike rate of place field less than zero')
                    end
                    tempt{track_id}(k,:) = single_place_field;
                end
            end
            
            all_place_fields{1}{event} = [tempt{1} tempt{2}];
        end
        
    elseif strcmp(modify,'rate remapping')
        cd 'rate_remapped'
        load extracted_place_fields_BAYESIAN
        cd ..
        
        for track_id = 1:2
            for k=1:length(place_field_index)
                single_place_field = place_fields_BAYESIAN.track(track_id).raw{place_field_index(k)}; %get raw place field
                single_place_field(find(isnan(single_place_field))) = 0; % remove NaNs in place field and replace by 0
                % single_place_field = smooth(single_place_field,parameters.smoothing_number_of_bins); %smooth place field
                
                if min(single_place_field)<0
                    disp('error- spike rate of place field less than zero')
                end
                tempt{track_id}(k,:) = single_place_field;
            end
        end
        
        all_place_fields{1} = [tempt{1} tempt{2}];
    end
    
    % Turning curve shuffle (T1 and T2 label)
    if strcmp(modify,'global remapping per event') | strcmp(modify,'PRE place_field_circular_shift') | strcmp(modify,'place_field_shifted') | strcmp(modify,'cross_experiment_shuffled') 
        % If global remapping or PRE place field circular shift, then do it
        % for each event.
        for event = 1:length(bayesian_spike_count.replay_events)
            for k = 1:length(place_field_index) % Turning curve shuffled for each cell (across two tracks)
                track_id = randperm(2);
                all_place_fields{2}{event}(k,:) = [tempt{track_id(1)}(k,:) tempt{track_id(2)}(k,:)];
            end
        end
    else
        for k = 1:length(place_field_index) % Turning curve shuffled for each cell (across two tracks)
            track_id = randperm(2);
            all_place_fields{2}(k,:) = [tempt{track_id(1)}(k,:) tempt{track_id(2)}(k,:)];
        end
    end
end

load decoded_replay_events
load decoded_replay_events_segments
% Apply formula of bayesian decoding
replay_id = bayesian_spike_count.replay_events_indices;

if strcmp(modify,'place field shift') | strcmp(modify,'global remapping per event')...
        | strcmp(modify,'rate remapping') | strcmp(modify,'track label swap') | strcmp(modify,'place_field_shifted') | strcmp(modify,'cross_experiment_shuffled') 
    if isempty(option)
        estimated_position(1).replay = track_reconstruct(decoded_replay_events,n.replay,all_place_fields{1},place_field_index,replay_id,bin_width,'per event');
    elseif strcmp(option,'ratemap shuffle')
        estimated_position(1).replay = track_reconstruct(decoded_replay_events,n.replay,all_place_fields{2},place_field_index,replay_id,bin_width,'per event');
        
    end
else
    if isempty(option)
        estimated_position(1).replay = track_reconstruct(decoded_replay_events,n.replay,all_place_fields{1},place_field_index,replay_id,bin_width,'');
    elseif strcmp(option,'ratemap shuffle')
        estimated_position(1).replay = track_reconstruct(decoded_replay_events,n.replay,all_place_fields{2},place_field_index,replay_id,bin_width,'');
    end
end

decoded_replay_events = [];
% save estimated_position estimated_position


%%%%%% NORMALIZING  %%%%%%%
%       columns need to sum to 1 (total probability across positions).
for track_id = 1
    
    summed_probability_replay = zeros(1,size(estimated_position(1).replay,2));
    total_posbins = length(estimated_position(1).position_bin_centres);
    
    % normalize probability across all position bins (both tracks )at each time bin would sum to 1
    
    summed_probability_replay = sum(estimated_position(track_id).replay,1);
    estimated_position(track_id).replay_Normalized_per_time_bin = estimated_position(track_id).replay...
        ./summed_probability_replay;
    estimated_position(track_id).replay_Normalized = estimated_position(track_id).replay...
        ./summed_probability_replay;
end

% If running replay events separately, also extract individual events from the estimated position matrix
if isfield(bayesian_spike_count,'replay_events')
    for track_id = 1

        for event = 1 : length(bayesian_spike_count(1).replay_events)
            thisReplay_indxs = find(bayesian_spike_count.replay_events_indices == event);
            % Normalize probability by the total sum of the probability
            % within each event such that probability should sum to 1
            estimated_position(track_id).replay_Normalized(:,thisReplay_indxs) = estimated_position(track_id).replay_Normalized(:,thisReplay_indxs)./...
                sum(sum(estimated_position(track_id).replay_Normalized(:,thisReplay_indxs)));
            estimated_position(track_id).replay_events(event).replay = estimated_position(track_id).replay_Normalized(:,thisReplay_indxs);
            % Calculate whole event log odd ((sum of probability of T1)/(sum of probability of T2))
            estimated_position(track_id).replay_events(event).probability_ratio = sum(sum(estimated_position(track_id).replay_Normalized(1:total_posbins/2,thisReplay_indxs)))...
                /sum(sum(estimated_position(track_id).replay_Normalized(total_posbins/2+1:end,thisReplay_indxs)));

            % Calculate half event log odd
            %             estimated_position(track_id).replay_events(event).probability_ratio_first = sum(sum(estimated_position(track_id).replay_Normalized(1:total_posbins/2,thisReplay_indxs(1):thisReplay_indxs(round(length(thisReplay_indxs)/2)))))...
            %                 /sum(sum(estimated_position(track_id).replay_Normalized(total_posbins/2+1:end,thisReplay_indxs(1):thisReplay_indxs(round(length(thisReplay_indxs)/2)))));
            %
            %             estimated_position(track_id).replay_events(event).probability_ratio_second = sum(sum(estimated_position(track_id).replay_Normalized(1:total_posbins/2,thisReplay_indxs(end)-round(length(thisReplay_indxs)/2):thisReplay_indxs(end))))...
            %                 /sum(sum(estimated_position(track_id).replay_Normalized(total_posbins/2+1:end,thisReplay_indxs(end)-round(length(thisReplay_indxs)/2):thisReplay_indxs(end))));
            timebin_index1 = thisReplay_indxs(decoded_replay_events1(1).replay_events(event).timebins_index);
            
            estimated_position(track_id).replay_events(event).probability_ratio_first = sum(sum(estimated_position(track_id).replay_Normalized(1:total_posbins/2,timebin_index1)))...
                /sum(sum(estimated_position(track_id).replay_Normalized(total_posbins/2+1:end,timebin_index1)));

            timebin_index2 = thisReplay_indxs(decoded_replay_events2(1).replay_events(event).timebins_index);

            estimated_position(track_id).replay_events(event).probability_ratio_second = sum(sum(estimated_position(track_id).replay_Normalized(1:total_posbins/2,timebin_index2)))...
                /sum(sum(estimated_position(track_id).replay_Normalized(total_posbins/2+1:end,timebin_index2)));

        end
    end
end


if strcmp(save_option, 'Y')
    save estimated_position estimated_position
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function estimated_position = track_reconstruct(decoded_replay_events,n,all_place_fields,place_field_index,replay_id,bin_width,option)

if strcmp(option,'per event')

    for event = 1 : length(decoded_replay_events(1).replay_events)
        all_place_fields_this_event = all_place_fields{event};
        thisReplay_indxs = find(replay_id == event);

        % Creates matrix where rows are cells and columns are position bins
        bin_length = size(all_place_fields_this_event,2); %columns (position bins)
        number_of_cells = size(all_place_fields_this_event,1); %rows (cells)
        parameters.bayesian_threshold=10.^(log2(number_of_cells)-log2(400)); % small value multiplied to all values to get rid of zeros
        all_place_fields_this_event(find(all_place_fields_this_event<parameters.bayesian_threshold)) = parameters.bayesian_threshold;
        sum_of_place_fields = sum(all_place_fields_this_event,1);  % adds up spikes per posiiton bin (used later for exponential)

        for j = 1: length(thisReplay_indxs)
            n_spikes = n(:,thisReplay_indxs(j))*ones(1,bin_length); %number of spikes in time bin
            %             n_spikes = n_spikes(active_cell_index,:);
            pre_product = all_place_fields_this_event.^n_spikes; % pl field values raised to num of spikes
            pre_product(find(pre_product<parameters.bayesian_threshold)) = parameters.bayesian_threshold;
            product_of_place_fields = prod(pre_product,1); %product of pl fields
            if length(bin_width) == 1 % if specifying only one universla time bin (i.e 20ms) for all replay events
                estimated_position(:,thisReplay_indxs(j)) = product_of_place_fields.*(exp(-bin_width*sum_of_place_fields)); % bayesian formula
                %NOTE- columns do not sum to 1.  this is done at a later stage to allow normalization within a track or across tracks
            else
                for m = 1: length(bin_width{event}) % if different time bin for each replay event
                    estimated_position(:,thisReplay_indxs(j)) = product_of_place_fields.*(exp(-bin_width{event}(m)*sum_of_place_fields)); % bayesian formula
                end
            end
        end

        thisReplay_indxs = [];
    end

else
    % Creates matrix where rows are cells and columns are position bins
    bin_length = size(all_place_fields,2); %columns (position bins)
    number_of_cells = size(all_place_fields,1); %rows (cells)
    parameters.bayesian_threshold=10.^(log2(number_of_cells)-log2(400)); % small value multiplied to all values to get rid of zeros
    all_place_fields(find(all_place_fields<parameters.bayesian_threshold)) = parameters.bayesian_threshold;
    sum_of_place_fields = sum(all_place_fields,1);  % adds up spikes per posiiton bin (used later for exponential)

    for event = 1 : length(decoded_replay_events(1).replay_events)
        thisReplay_indxs = find(replay_id == event);


        for j = 1: length(thisReplay_indxs)
            n_spikes = n(:,thisReplay_indxs(j))*ones(1,bin_length); %number of spikes in time bin
            %             n_spikes = n_spikes(active_cell_index,:);
            pre_product = all_place_fields.^n_spikes; % pl field values raised to num of spikes
            pre_product(find(pre_product<parameters.bayesian_threshold)) = parameters.bayesian_threshold;
            product_of_place_fields = prod(pre_product,1); %product of pl fields
            if length(bin_width) == 1 % if specifying only one universla time bin (i.e 20ms) for all replay events
                estimated_position(:,thisReplay_indxs(j)) = product_of_place_fields.*(exp(-bin_width*sum_of_place_fields)); % bayesian formula
                %NOTE- columns do not sum to 1.  this is done at a later stage to allow normalization within a track or across tracks
            else
                for m = 1: length(bin_width{event}) % if different time bin for each replay event
                    estimated_position(:,thisReplay_indxs(j)) = product_of_place_fields.*(exp(-bin_width{event}(m)*sum_of_place_fields)); % bayesian formula
                end
            end
        end

        thisReplay_indxs = [];
    end
end

end


function estimated_position = reconstruct(n,all_place_fields,bin_width)
% Creates matrix where rows are cells and columns are position bins
bin_length = size(all_place_fields,2); %columns
number_of_cells = size(all_place_fields,1); %rows
parameters.bayesian_threshold=10.^(log2(number_of_cells)-log2(400)); % small value multiplied to all values to get rid of zeros
all_place_fields(find(all_place_fields<parameters.bayesian_threshold)) = parameters.bayesian_threshold;
sum_of_place_fields = sum(all_place_fields,1);  % adds up spikes per bin (used later for exponential)
for j = 1: size(n,2)
    n_spikes = n(:,j)*ones(1,bin_length); %number of spikes in time bin
    pre_product = all_place_fields.^n_spikes; % pl field values raised to num of spikes
    pre_product(find(pre_product<parameters.bayesian_threshold)) = parameters.bayesian_threshold;
    product_of_place_fields = prod(pre_product,1); %product of pl fields
    estimated_position(:,j) = product_of_place_fields.*(exp(-bin_width*sum_of_place_fields)); % bayesian formula
    %NOTE- columns do not sum to 1.  this is done at a later stage to allow normalization within a track or across tracks
end
end


