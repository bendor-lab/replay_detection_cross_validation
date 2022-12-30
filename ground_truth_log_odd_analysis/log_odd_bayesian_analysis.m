function [log_odd] = log_odd_bayesian_analysis(folders,timebin,posbin)
% Input:
% Folders: Each Folder is data for a session
% Timebin: 0.02 or 1
% 0.02 -> 20ms time bin
% 1 -> 1 event time bin
% posbin: [] or 1
% [] -> 20 position bin per track (will load extracted_place_field_BAYESIAN)
% 1 -> 1 position bin per track

% Output:
% log_odd: save the log odd - log(sum(probability of T1)/sum(probability of
% T2)) of the events that are considered significant according to the original bayesian decoding
% as well as the 1000 T1&T2 turning curve shuffles for each event
% It will also contain the session ID, 'ground truth' replay track ID and
% behavioral epoch (-1: PRE, 0: POST, 1: T1 RUN, 2: T2 RUN) and etc

c = 1;
current_directory=pwd;
log_odd=[];

for f = 1:length(folders)
    cd Tables
    load subsets_of_cells;
    cd ..
    
    cd(folders{f})
    parameters = list_of_parameters;
    load('significant_replay_events_wcorr.mat');
    load('sorted_replay_wcorr');
    load('extracted_place_fields_BAYESIAN');
    load('extracted_position')
    load('extracted_clusters.mat');
    load('extracted_replay_events.mat');
    
    % Stable common good cells (based on subsets_of_cells)
    place_cell_index = subset_of_cells.cell_IDs{1}  ...
        (intersect(find(subset_of_cells.cell_IDs{1} >= 1000*f),find(subset_of_cells.cell_IDs{1} <= 1000*(f+1))))- f*1000;
    
    
    %% Get Original Replay Spike Count for the place cells of interest
    load decoded_replay_events
    load extracted_replay_events
    replay_starts = replay.onset;
    replay_ends = replay.offset;
    
    % Get time vectors for bayesian decoding and matrix with spike count
    disp('Spike count...');
    a = track_spike_count([],replay_starts,replay_ends,timebin,place_cell_index,'N');
    
    
%% Track Field (if one bin)
%     
%     if isnumeric(posbin) && ~isempty(posbin) % if time bin is not 20ms, make 'Track' Field
%         Track_fields = [];
%         x_bin_centres = []
%         time_bin_width=position.t(2)-position.t(1);
%         
%         for track_id = 1:length(position.linear)
%             
%             position_index = isnan(position.linear(track_id).linear);
%             position_speed = abs(position.v_cm);
%             position_speed(position_index) = NaN;  %make sure speed is NaN if position is NaN
%             
%             position_during_spike = interp1(position.t,position.linear(track_id).linear,clusters.spike_times,'nearest'); %interpolates position into spike time
%             speed_during_spike = interp1(position.t,position_speed,clusters.spike_times,'nearest');
%             
%             x_bin_centres = [];
%             x_bin_edges = linspace(0,200,posbin+1);
%             
%             for bin = 1:length(x_bin_edges)-1
%                 x_bin_centres(bin) = (x_bin_edges(bin)+x_bin_edges(bin+1))/2;
%             end
%             
%             % Time spent at each x_bin (speed filtered)
%             x_hist = time_bin_width.*histcounts(position.linear(track_id).linear(find(position.t>place_fields_BAYESIAN.track(track_id).time_window(1) &...
%                 position.t<place_fields_BAYESIAN.track(track_id).time_window(2) & position_speed>parameters.speed_threshold_laps...
%                 & position_speed<parameters.speed_threshold_max)),x_bin_edges); % Changed bin_centre to bin_edges
%             
%             Track_fields(track_id).dwell_map = x_hist;
%             
%             for j = 1 : max(clusters.id_conversion(:,1))
%                 % Number of spikes per bin within time window (speed filtered)
%                 Track_fields(track_id).spike_hist(:,j) = histcounts(position_during_spike(find(clusters.spike_id==j & ...
%                     clusters.spike_times>place_fields_BAYESIAN.track(track_id).time_window(1) & ...
%                     clusters.spike_times<place_fields_BAYESIAN.track(track_id).time_window(2) & ...
%                     speed_during_spike>parameters.speed_threshold_laps &...
%                     speed_during_spike<parameters.speed_threshold_max)),x_bin_edges); % Changed bin_centre to bin_edges
%                 
%                 Track_fields(track_id).mean_rate(:,j) = Track_fields(track_id).spike_hist(:,j)./x_hist'; % place field calculation
%                 if Track_fields(track_id).mean_rate(:,j) == nan
%                     Track_fields(track_id).mean_rate(:,j) = 0;
%                 end
%                 
%                 % Number of spikes per bin within time window (speed filtered)
%                 [FR, ID] = max(Track_fields(track_id).mean_rate(:,j)); % raw pl field peak
%                 Track_fields(track_id).peak(j,1) = ID;
%                 Track_fields(track_id).peak(j,2) = FR;
%             end
%             
%             
%             Track_fields(track_id).common_good_cells_ID = place_cell_index;
%             Track_fields(track_id).good_cells_ID = place_fields_BAYESIAN.good_place_cells;           
%         end
%         
%         % Rate fixed version
%         for track_id = 1:2
%             for j = 1 : length(place_fields_BAYESIAN.track(1).peak)
%                 % Average of the peak rates on both tracks (For one bin, it is just the average of the average FR on both tracks)
%                 average_mean_rate = mean([max(Track_fields(1).mean_rate(j)) max(Track_fields(2).mean_rate(j))]);
%                 
%                 % If it's a good cell, calculate scaling factor
%                 if isempty(find(place_fields_BAYESIAN.good_place_cells==j)) % if it's not a good cell, then don't change
%                     scaling_factor = 1;
%                 else
%                     scaling_factor = average_mean_rate./max(Track_fields(track_id).mean_rate(j));
%                 end
%                 
%                 Track_fields(track_id).rate_fixed_scaling_factor(:,j)= scaling_factor;
%                 
%                 % Multiply the place field by its scaling factor
%                 Track_fields(track_id).rate_fixed_mean_rate(:,j) = scaling_factor*Track_fields(track_id).mean_rate(:,j)
%             end
%         end
%         
%         
%         % Rate remapped version
%         for event = 1:length(a.replay_events)
%             for track_id = 1:2
%                 % Sort good cells by raw peak firing rate
%                 [raw_peak_distribution,~] = sort(Track_fields(track_id).mean_rate(:,(place_fields_BAYESIAN.track(track_id).good_cells)));
%                 fraction_of_cells = (0:(length(place_fields_BAYESIAN.track(track_id).good_cells)-1))/(length(place_fields_BAYESIAN.track(track_id).good_cells)-1);
%                 
%                 for j = 1 : length(place_fields_BAYESIAN.track(1).sorted)
%                     % Picks a random value within the real distribution of raw peak FR
%                     random_peak_rate = interp1(fraction_of_cells,raw_peak_distribution,rand(1),'linear');
%                     
%                     % If it's a good cell, calculate scaling factor
%                     if isempty(find(place_fields_BAYESIAN.track(track_id).good_cells==j)) % if it's not a good cell, then don't change
%                         scaling_factor = 1;
%                     else
%                         scaling_factor = random_peak_rate./max(Track_fields(track_id).mean_rate(j));
%                     end
%                     
%                     Track_fields(track_id).rate_remapped_scaling_factor{event}(:,j)= scaling_factor;
%                     
%                     % Multiply each spike and the raw peak FR by it's scaling factor
%                     Track_fields(track_id).rate_remapped_mean_rate{event}(:,j) = scaling_factor*Track_fields(track_id).mean_rate(:,j);
%                 end
%             end
%         end
%          
        
%         save Track_fields Track_fields
%     end
    
    
    %% bayesian decoding

    % Saving replay spike count
    replayEvents_common_good_cell_bayesian_spike_count = a;
    save replayEvents_common_good_cell_bayesian_spike_count replayEvents_common_good_cell_bayesian_spike_count
    
% 
      if timebin < 1 % If 20ms (or other ms timebin)
%         %% Rate fix
%         for track_id = 1:2
%             for j = 1 : length(place_fields_BAYESIAN.track(1).peak)
%                 % Average of the peak rates on both tracks
%                 average_mean_rate = mean([max(place_fields_BAYESIAN.track(1).raw{j}) max(place_fields_BAYESIAN.track(2).raw{j})]);
%                 
%                 % If it's a good cell, calculate scaling factor
%                 if isempty(find(place_fields_BAYESIAN.good_place_cells==j)) % if it's not a good cell, then don't change
%                     scaling_factor = 1;
%                 else
%                     scaling_factor = average_mean_rate./max(place_fields_BAYESIAN.track(track_id).raw{j});
%                 end
%                 
%                 place_fields_BAYESIAN.track(track_id).rate_fixed_scaling_factor(j) = scaling_factor;
%                 
%                 % Multiply the Place field by the scaling factor
%                 place_fields = place_fields_BAYESIAN.track(track_id).raw{j}
%                 %                         place_fields(find(place_fields<parameters.bayesian_threshold)) = parameters.bayesian_threshold;
%                 place_fields_BAYESIAN.track(track_id).rate_fixed_raw{j} = scaling_factor*place_fields;
%             end
%         end
%         save extracted_place_fields_BAYESIAN place_fields_BAYESIAN
%         
%         % Rate fixed spike count (replay spike count scaled by average replay Firing Rate and place field rate fixed)
        rate_fixed_spike_count('replayEvents_common_good_cell_bayesian_spike_count',timebin,'')
        bayesian_spike_count = 'replayEvents_FIXED_spike_count_20ms';
        
%         % Calculate the probability ratio
         estimated_position_rate_fixed_20ms = log_odd_bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count,place_cell_index,timebin,[],'','rate fixed','N');
%         
%         % Save the probability ratio
%         for event = 1:length(estimated_position_rate_fixed_20ms(1).replay_events)
%             rate_fixed_20ms_probability_ratio{1}(event) = estimated_position_rate_fixed_20ms(1).replay_events(event).probability_ratio;
%             disp('Rate Fix (20ms)...')
%         end
%         
%         % Calculate the T1-T2 ratemap (turning curve) shuffle 1000 times
%         % such that there are 1000 shuffles for each event.
%         parfor nshuffle = 1:1000
%             estimated_position_ratemap_shuffled = [];
%             estimated_position_ratemap_shuffled = log_odd_bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count,place_cell_index,timebin,[],'ratemap shuffle','rate fixed','N');
%             for event = 1:length(estimated_position_ratemap_shuffled(1).replay_events)
%                 ratemap_shuffled_probability_ratio{nshuffle}(event) = estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio;
%                 disp('Ratemap shuffling (rate fixed 20ms)...');
%             end
%         end
%         
%         % Save the shuffled probability ratio
%         rate_fixed_20ms_probability_ratio{2} = ratemap_shuffled_probability_ratio;
%         
%         if ~isfolder({'rate_fixed'})
%             mkdir('rate_fixed')
%         end
%         
         cd 'rate_fixed'
         save estimated_position_rate_fixed_20ms estimated_position_rate_fixed_20ms
%         save rate_fixed_20ms_probability_ratio rate_fixed_20ms_probability_ratio
         cd ..
%         
         clear estimated_position_rate_fixed_20ms rate_fixed_20ms_probability_ratio ratemap_shuffled_probability_ratio
%         
        
%         %%    Rate fixed global remapping (per event)
        bayesian_spike_count = 'replayEvents_FIXED_spike_count_20ms';
%         for event = 1:length(a.replay_events)
%             % global remapping by swapping the Cell ID for good cells
%             for track_id = 1:2
%                 global_remapped_place_fields{track_id} = place_fields_BAYESIAN.track(track_id).raw;
%                 
%                 random_cell_index=randperm(length(place_fields_BAYESIAN.good_place_cells));
%                 random_cell=place_fields_BAYESIAN.good_place_cells(random_cell_index);
%                 %                         place_fields_BAYESIAN.track(track_id).random_cell = random_cell;
%                 original_cell=place_fields_BAYESIAN.good_place_cells;
%                 
%                 for j=1:length(random_cell) %only swap good cells
%                     global_remapped_place_fields{track_id}{original_cell(j)}=place_fields_BAYESIAN.track(track_id).raw{random_cell(j)};
%                 end
%                 
%             end
%             
%             % After Cell ID shuffle, rate fix the place fields
%             for track_id = 1:2
%                 for j = 1 : length(place_fields_BAYESIAN.track(1).peak)
%                     % Average of the peak rates on both tracks
%                     average_mean_rate = mean([max(global_remapped_place_fields{1}{j})...
%                         max(global_remapped_place_fields{2}{j})]);
%                     
%                     % If it's a good cell, calculate scaling factor
%                     if isempty(find(place_fields_BAYESIAN.good_place_cells==j)) % if it's not a good cell, then don't change
%                         scaling_factor = 1;
%                     else
%                         scaling_factor = average_mean_rate./max(global_remapped_place_fields{track_id}{j});
%                     end
%                     
%                     % Multiply the Place field by the scaling factor
%                     place_fields = global_remapped_place_fields{track_id}{j};
%                     place_fields_BAYESIAN.track(track_id).global_remapped_raw{event}{j} = scaling_factor*place_fields;
%                 end
%             end
%             
%         end
%         
%         save extracted_place_fields_BAYESIAN place_fields_BAYESIAN
%         
%         % Calculate probability ratio
         estimated_position_global_remapped_20ms = log_odd_bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count,place_cell_index,timebin,[],'','global remapping per event','N');
%         
%         for event = 1:length(estimated_position_global_remapped_20ms(1).replay_events)
%             global_remapped_20ms_probability_ratio{1}(event) = estimated_position_global_remapped_20ms(1).replay_events(event).probability_ratio;
%             disp('Rate Fixed Global Remapped (20ms)...')
%         end
%         
%         % Calculate the T1-T2 ratemap (turning curve) shuffle 1000 times
%         % such that there are 1000 shuffles for each event.
%         parfor nshuffle = 1:1000
%             estimated_position_ratemap_shuffled = [];
%             estimated_position_ratemap_shuffled = log_odd_bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count,place_cell_index,timebin,[],'ratemap shuffle','global remapping per event','N');
%             for event = 1:length(estimated_position_ratemap_shuffled(1).replay_events)
%                 ratemap_shuffled_probability_ratio{nshuffle}(event) = estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio;
%                 disp('Ratemap shuffling (global remapped)...');
%             end
%         end
%         
%         global_remapped_20ms_probability_ratio{2} = ratemap_shuffled_probability_ratio;
%         
         cd 'rate_fixed'
         save estimated_position_global_remapped_20ms estimated_position_global_remapped_20ms
%         save global_remapped_20ms_probability_ratio global_remapped_20ms_probability_ratio
         cd ..
%         
         clear estimated_position_global_remapped_20ms global_remapped_20ms_probability_ratio ratemap_shuffled_probability_ratio
%         
        
%         %% Normal 20ms
        bayesian_spike_count = 'replayEvents_common_good_cell_bayesian_spike_count';
         estimated_position_20ms = log_odd_bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count,place_cell_index,timebin,[],'','','N');
%         
%         for event = 1:length(estimated_position_20ms(1).replay_events)
%             probability_ratio_20ms{1}(event) = estimated_position_20ms(1).replay_events(event).probability_ratio;
%             disp('20ms decoding...')
%         end
%         
%         parfor nshuffle = 1:1000
%             estimated_position_ratemap_shuffled = [];
%             estimated_position_ratemap_shuffled = log_odd_bayesian_decoding(place_fields_BAYESIAN,bayesian_spike_count,place_cell_index,timebin,[],'ratemap shuffle','','N');
%             for event = 1:length(estimated_position_ratemap_shuffled(1).replay_events)
%                 ratemap_shuffled_probability_ratio{nshuffle}(event) = estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio;
%             end
%             disp('Ratemap shuffling (20ms)...')
%         end
%         
%         probability_ratio_20ms{2} = ratemap_shuffled_probability_ratio
%         
%         
%         save probability_ratio_20ms probability_ratio_20ms
         save estimated_position_20ms estimated_position_20ms
%         
         clear estimated_position_20ms probability_ratio_20ms
%         
      elseif timebin >= 1 % if one event time bin (or more than one time bin per event)
%         %% Rate Remapped (rate randomised) 1 time bin (place info removed)
        bayesian_spike_count = 'replayEvents_common_good_cell_bayesian_spike_count';   
         estimated_position_rate_remapped_one_bin = log_odd_bayesian_decoding([],bayesian_spike_count,place_cell_index,timebin,posbin,'','rate remapping','N');
%         
%         for event = 1:length(estimated_position_rate_remapped_one_bin(1).replay_events)
%             probability_ratio_rate_remapped_one_bin{1}(event) = estimated_position_rate_remapped_one_bin(1).replay_events(event).probability_ratio;
%             disp('one bin decoding...')
%         end
%         
%         parfor nshuffle = 1:1000
%             estimated_position_ratemap_shuffled = [];
%             estimated_position_ratemap_shuffled = log_odd_bayesian_decoding([],bayesian_spike_count,place_cell_index,timebin,posbin,'ratemap shuffle','rate remapping','N');
%             for event = 1:length(estimated_position_ratemap_shuffled(1).replay_events)
%                 ratemap_shuffled_probability_ratio{nshuffle}(event) = estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio;
%                 disp('Ratemap shuffling (one bin)...');
%             end
%         end
%         
%         probability_ratio_rate_remapped_one_bin{2} = ratemap_shuffled_probability_ratio;
%         
%         if ~isfolder({'rate remapping'})
%             mkdir('rate remapping')
%         end
%         
         cd 'rate remapping'
%         save probability_ratio_rate_remapped_one_bin probability_ratio_rate_remapped_one_bin
         save estimated_position_rate_remapped_one_bin estimated_position_rate_remapped_one_bin
         cd ..
%         
         clear probability_ratio_rate_remapped_one_bin estimated_position_rate_remapped_one_bin ratemap_shuffled_probability_ratio
%         
        
%         %% 1 time bin (place info removed)
        bayesian_spike_count = 'replayEvents_common_good_cell_bayesian_spike_count';      
         estimated_position_one_bin = log_odd_bayesian_decoding([],bayesian_spike_count,place_cell_index,timebin,posbin,'','','N');
%         
%         for event = 1:length(estimated_position_one_bin(1).replay_events)
%             probability_ratio_one_bin{1}(event) = estimated_position_one_bin(1).replay_events(event).probability_ratio;
%             disp('one bin decoding...')
%         end
%         
%         parfor nshuffle = 1:1000
%             estimated_position_ratemap_shuffled = [];
%             estimated_position_ratemap_shuffled = log_odd_bayesian_decoding([],bayesian_spike_count,place_cell_index,timebin,posbin,'ratemap shuffle','','N');
%             for event = 1:length(estimated_position_ratemap_shuffled(1).replay_events)
%                 ratemap_shuffled_probability_ratio{nshuffle}(event) = estimated_position_ratemap_shuffled(1).replay_events(event).probability_ratio;
%                 disp('Ratemap shuffling (one bin)...');
%             end
%         end
%         
%         probability_ratio_one_bin{2} = ratemap_shuffled_probability_ratio;
%         
%         
%         save probability_ratio_one_bin probability_ratio_one_bin
         save estimated_position_one_bin estimated_position_one_bin
         clear estimated_position_one_bin probability_ratio_one_bin ratemap_shuffled_probability_ratio
        
end

    
    
    %% Ground Truth Information (Bayesian bias, p value, best score, probability ratio)
    
    load scored_replay
    load scored_replay_segments
    method = 'wcorr';
    % gather p value and replay score information 
    [p_values, replay_scores] = extract_score_and_pvalue(scored_replay, scored_replay1, scored_replay2, method, []);
    
    % gather time information for sorting the event as PRE, POST or RUN
    [~,time_range]=sort_replay_events([],'wcorr');
    
   
    
    if timebin >= 1
        load probability_ratio_one_bin
        
        if ~isfolder({'rate remapping'})
            mkdir('rate remapping')
        end
        
%         cd 'rate_fixed'
%         load rate_fixed_one_bin_probability_ratio
%         cd ..
        
        cd 'rate remapping'
        load probability_ratio_rate_remapped_one_bin
        cd ..
        
    elseif timebin < 1
        load probability_ratio_20ms
        
        if ~isfolder({'rate_fixed'})
            mkdir('rate_fixed')
        end
        
        cd 'rate_fixed'
        load rate_fixed_20ms_probability_ratio
        load global_remapped_20ms_probability_ratio
        cd ..
    end
    
    for track=1:2
        for i=1:length(significant_replay_events.track(track).p_value)
            if ~ismember(significant_replay_events.track(track).ref_index(i),significant_replay_events.multi_tracks_index) % if not multi-track
                % save replay event index
                log_odd.index(c) = significant_replay_events.track(track).ref_index(i);
                
                if timebin >= 1
                    % One bin
                    log_odd.one_bin.probability_ratio(c) = probability_ratio_one_bin{1}(log_odd.index(c));
                    for nshuffle = 1:1000
                        log_odd.one_bin.probability_ratio_shuffled{c}(nshuffle) = probability_ratio_one_bin{2}{nshuffle}(log_odd.index(c));
                    end
                    
                     % Calculate and save zscored log odd
                    data = log(log_odd.one_bin.probability_ratio(c));
                    shuffled_data = log(log_odd.one_bin.probability_ratio_shuffled{c});
                     
                    if data-mean(shuffled_data) == 0
                        tempt = 0;
                    else
                        tempt = (data-mean(shuffled_data))/std(shuffled_data);
                    end
                    
                    log_odd.one_bin_zscored.original(c) = tempt;
                   
                    
                    % One bin rate remapped
                    log_odd.one_bin.rate_remapped_probability_ratio(c) = probability_ratio_rate_remapped_one_bin{1}(log_odd.index(c));
                    for nshuffle = 1:1000
                        log_odd.one_bin.rate_remapped_probability_ratio_shuffled{c}(nshuffle) = probability_ratio_rate_remapped_one_bin{2}{nshuffle}(log_odd.index(c));
                    end
                    
                     % Calculate and save zscored log odd
                    data = log(log_odd.one_bin.rate_remapped_probability_ratio(c));
                    shuffled_data = log(log_odd.one_bin.rate_remapped_probability_ratio_shuffled{c});
                     
                    if data-mean(shuffled_data) == 0
                        tempt = 0;
                    else
                        tempt = (data-mean(shuffled_data))/std(shuffled_data);
                    end
                    
                    log_odd.one_bin_zscored.rate_remapping(c) = tempt;                    
                    
                    
                elseif timebin < 1
                    
                    % Original 20ms                    
                    log_odd.normal.probability_ratio(c) = probability_ratio_20ms{1}(log_odd.index(c));
                    
                    for nshuffle = 1:1000
                        log_odd.normal.probability_ratio_shuffled{c}(nshuffle) = probability_ratio_20ms{2}{nshuffle}(log_odd.index(c));
                    end
                    
                    % Calculate and save zscored log odd
                    data = log(log_odd.normal.probability_ratio(c));
                    shuffled_data = log(log_odd.normal.probability_ratio_shuffled{c});
                    tempt = (data-mean(shuffled_data))/std(shuffled_data);
                    log_odd.normal_zscored.original(c) = tempt;
                    
                    
                    % 20ms rate fixed
                    log_odd.normal.rate_fixed_20ms_probability_ratio(c) = rate_fixed_20ms_probability_ratio{1}(log_odd.index(c));
                    for nshuffle = 1:1000
                        log_odd.normal.rate_fixed_20ms_probability_ratio_shuffled{c}(nshuffle) =...
                            rate_fixed_20ms_probability_ratio{2}{nshuffle}(log_odd.index(c));
                    end
                    
                    % Calculate and save zscored log odd
                    data = log(log_odd.normal.rate_fixed_20ms_probability_ratio(c));
                    shuffled_data = log(log_odd.normal.rate_fixed_20ms_probability_ratio_shuffled{c});
                    tempt = (data-mean(shuffled_data))/std(shuffled_data);
                    log_odd.normal_zscored.rate_fixed(c) = tempt;
                    
                                      
                    % 20ms rate fixed global remapped
                    log_odd.normal.rate_fixed_global_remapping_probability_ratio(c) = global_remapped_20ms_probability_ratio{1}(log_odd.index(c));
                    for nshuffle = 1:1000
                        log_odd.normal.rate_fixed_global_remapping_probability_ratio_shuffled{c}(nshuffle) =...
                            global_remapped_20ms_probability_ratio{2}{nshuffle}(log_odd.index(c));
                    end
                    % Calculate and save zscored log odd
                    data = log(log_odd.normal.rate_fixed_global_remapping_probability_ratio(c));
                    shuffled_data = log(log_odd.normal.rate_fixed_global_remapping_probability_ratio_shuffled{c});
                    tempt = (data-mean(shuffled_data))/std(shuffled_data);
                    log_odd.normal_zscored.rate_fixed_global_remapping(c) = tempt;
                    
                    
                end
                
                
                log_odd.pvalue(c)=significant_replay_events.track(track).p_value(i);
                log_odd.best_score(c)=significant_replay_events.track(track).replay_score(i);
                log_odd.track1_best_score(c)=replay_scores.WHOLE(1,log_odd.index(c));
%                 log_odd.active_cells(c)= length(decoded_track_replay(1).replay_active_cells{log_odd.index(c)});
                log_odd.event_duration(c) = significant_replay_events.track(track).event_duration(i);
                
                % The event time of the current event
                event_time = significant_replay_events.track(track).event_times(i);
                behavioural_state = [];
                
                
                if event_time < time_range.pre(2) & event_time > time_range.pre(1) %If PRE
                    
                    log_odd.behavioural_state(c)= -1;
                    
                elseif event_time>time_range.post(1) & event_time < time_range.post(2) %If POST
                    
                    log_odd.behavioural_state(c)= 0;
                   
                    
                elseif event_time>time_range.track(1).behaviour(1) & event_time<time_range.track(1).behaviour(2)
                    log_odd.behavioural_state(c)=1;  % If Run T1
                elseif event_time>time_range.track(2).behaviour(1) & event_time<time_range.track(2).behaviour(2)
                    log_odd.behavioural_state(c)=2;  % If Run T2
                    
                    
                else
                    log_odd.behavioural_state(i)=NaN;  %boundary between behavioural states
                end
                
                % The 'ground truth' track ID label
                log_odd.track(c)=track;
                
                % The session ID
                log_odd.experiment(c)=f;
                c=c+1;
            end
        end
    end
    
    clear global_remapped_probability_ratio
    clear rate_remapped_probability_ratio
    
    clear probability_ratio
    clear estimated_position
    
    clear global_remapped_20ms_probability_ratio
    clear probability_ratio_20ms
    clear probability_ratio_one_bin
%     clear rate_fixed_one_bin_probability_ratio
    clear probability_ratio_rate_remapped_one_bin
    clear rate_fixed_20ms_probability_ratio
    
    cd(current_directory)
end

end



